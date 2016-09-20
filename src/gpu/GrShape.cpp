/*
 * Copyright 2016 Google Inc.
 *
 * Use of this source code is governed by a BSD-style license that can be
 * found in the LICENSE file.
 */

#include "GrShape.h"

GrShape& GrShape::operator=(const GrShape& that) {
    fStyle = that.fStyle;
    this->changeType(that.fType, Type::kPath == that.fType ? &that.path() : nullptr);
    switch (fType) {
        case Type::kEmpty:
            break;
        case Type::kRRect:
            fRRectData = that.fRRectData;
            break;
        case Type::kLine:
            fLineData = that.fLineData;
            break;
        case Type::kPath:
            fPathData.fGenID = that.fPathData.fGenID;
            break;
    }
    fInheritedKey.reset(that.fInheritedKey.count());
    sk_careful_memcpy(fInheritedKey.get(), that.fInheritedKey.get(),
                      sizeof(uint32_t) * fInheritedKey.count());
    return *this;
}

SkRect GrShape::bounds() const {
    // Bounds where left == bottom or top == right can indicate a line or point shape. We return
    // inverted bounds for a truly empty shape.
    static constexpr SkRect kInverted = SkRect::MakeLTRB(1, 1, -1, -1);
    switch (fType) {
        case Type::kEmpty:
            return kInverted;
        case Type::kLine: {
            SkRect bounds;
            if (fLineData.fPts[0].fX < fLineData.fPts[1].fX) {
                bounds.fLeft = fLineData.fPts[0].fX;
                bounds.fRight = fLineData.fPts[1].fX;
            } else {
                bounds.fLeft = fLineData.fPts[1].fX;
                bounds.fRight = fLineData.fPts[0].fX;
            }
            if (fLineData.fPts[0].fY < fLineData.fPts[1].fY) {
                bounds.fTop = fLineData.fPts[0].fY;
                bounds.fBottom = fLineData.fPts[1].fY;
            } else {
                bounds.fTop = fLineData.fPts[1].fY;
                bounds.fBottom = fLineData.fPts[0].fY;
            }
            return bounds;
        }
        case Type::kRRect:
            return fRRectData.fRRect.getBounds();
        case Type::kPath:
            return this->path().getBounds();
    }
    SkFAIL("Unknown shape type");
    return kInverted;
}

SkRect GrShape::styledBounds() const {
    if (Type::kEmpty == fType && !fStyle.hasNonDashPathEffect()) {
        return SkRect::MakeEmpty();
    }
    SkRect bounds;
    fStyle.adjustBounds(&bounds, this->bounds());
    return bounds;
}

int GrShape::unstyledKeySize() const {
    if (fInheritedKey.count()) {
        return fInheritedKey.count();
    }
    switch (fType) {
        case Type::kEmpty:
            return 1;
        case Type::kRRect:
            SkASSERT(!fInheritedKey.count());
            SkASSERT(0 == SkRRect::kSizeInMemory % sizeof(uint32_t));
            // + 1 for the direction, start index, and inverseness.
            return SkRRect::kSizeInMemory / sizeof(uint32_t) + 1;
        case Type::kLine:
            GR_STATIC_ASSERT(2 * sizeof(uint32_t) == sizeof(SkPoint));
            // 4 for the end points and 1 for the inverseness
            return 5;
        case Type::kPath:
            if (0 == fPathData.fGenID) {
                return -1;
            } else {
                // The key is the path ID and fill type.
                return 2;
            }
    }
    SkFAIL("Should never get here.");
    return 0;
}

void GrShape::writeUnstyledKey(uint32_t* key) const {
    SkASSERT(this->unstyledKeySize());
    SkDEBUGCODE(uint32_t* origKey = key;)
    if (fInheritedKey.count()) {
        memcpy(key, fInheritedKey.get(), sizeof(uint32_t) * fInheritedKey.count());
        SkDEBUGCODE(key += fInheritedKey.count();)
    } else {
        switch (fType) {
            case Type::kEmpty:
                *key++ = 1;
                break;
            case Type::kRRect:
                fRRectData.fRRect.writeToMemory(key);
                key += SkRRect::kSizeInMemory / sizeof(uint32_t);
                *key = (fRRectData.fDir == SkPath::kCCW_Direction) ? (1 << 31) : 0;
                *key |= fRRectData.fInverted ? (1 << 30) : 0;
                *key++ |= fRRectData.fStart;
                SkASSERT(fRRectData.fStart < 8);
                break;
            case Type::kLine:
                memcpy(key, fLineData.fPts, 2 * sizeof(SkPoint));
                key += 4;
                *key++ = fLineData.fInverted ? 1 : 0;
                break;
            case Type::kPath:
                SkASSERT(fPathData.fGenID);
                *key++ = fPathData.fGenID;
                // We could canonicalize the fill rule for paths that don't differentiate between
                // even/odd or winding fill (e.g. convex).
                *key++ = this->path().getFillType();
                break;
        }
    }
    SkASSERT(key - origKey == this->unstyledKeySize());
}

void GrShape::setInheritedKey(const GrShape &parent, GrStyle::Apply apply, SkScalar scale) {
    SkASSERT(!fInheritedKey.count());
    // If the output shape turns out to be simple, then we will just use its geometric key
    if (Type::kPath == fType) {
        // We want ApplyFullStyle(ApplyPathEffect(shape)) to have the same key as
        // ApplyFullStyle(shape).
        // The full key is structured as (geo,path_effect,stroke).
        // If we do ApplyPathEffect we get get,path_effect as the inherited key. If we then
        // do ApplyFullStyle we'll memcpy geo,path_effect into the new inherited key
        // and then append the style key (which should now be stroke only) at the end.
        int parentCnt = parent.fInheritedKey.count();
        bool useParentGeoKey = !parentCnt;
        if (useParentGeoKey) {
            parentCnt = parent.unstyledKeySize();
            if (parentCnt < 0) {
                // The parent's geometry has no key so we will have no key.
                fPathData.fGenID = 0;
                return;
            }
        }
        uint32_t styleKeyFlags = 0;
        if (parent.knownToBeClosed()) {
            styleKeyFlags |= GrStyle::kClosed_KeyFlag;
        }
        if (parent.asLine(nullptr, nullptr)) {
            styleKeyFlags |= GrStyle::kNoJoins_KeyFlag;
        }
        int styleCnt = GrStyle::KeySize(parent.fStyle, apply, styleKeyFlags);
        if (styleCnt < 0) {
            // The style doesn't allow a key, set the path gen ID to 0 so that we fail when
            // we try to get a key for the shape.
            fPathData.fGenID = 0;
            return;
        }
        fInheritedKey.reset(parentCnt + styleCnt);
        if (useParentGeoKey) {
            // This will be the geo key.
            parent.writeUnstyledKey(fInheritedKey.get());
        } else {
            // This should be (geo,path_effect).
            memcpy(fInheritedKey.get(), parent.fInheritedKey.get(),
                   parentCnt * sizeof(uint32_t));
        }
        // Now turn (geo,path_effect) or (geo) into (geo,path_effect,stroke)
        GrStyle::WriteKey(fInheritedKey.get() + parentCnt, parent.fStyle, apply, scale,
                          styleKeyFlags);
    }
}

GrShape::GrShape(const GrShape& that) : fStyle(that.fStyle) {
    const SkPath* thatPath = Type::kPath == that.fType ? &that.fPathData.fPath : nullptr;
    this->initType(that.fType, thatPath);
    switch (fType) {
        case Type::kEmpty:
            break;
        case Type::kRRect:
            fRRectData = that.fRRectData;
            break;
        case Type::kLine:
            fLineData = that.fLineData;
            break;
        case Type::kPath:
            fPathData.fGenID = that.fPathData.fGenID;
            break;
    }
    fInheritedKey.reset(that.fInheritedKey.count());
    sk_careful_memcpy(fInheritedKey.get(), that.fInheritedKey.get(),
                      sizeof(uint32_t) * fInheritedKey.count());
}

GrShape::GrShape(const GrShape& parent, GrStyle::Apply apply, SkScalar scale) {
    // TODO: Add some quantization of scale for better cache performance here or leave that up
    // to caller?
    // TODO: For certain shapes and stroke params we could ignore the scale. (e.g. miter or bevel
    // stroke of a rect).
    if (!parent.style().applies() ||
        (GrStyle::Apply::kPathEffectOnly == apply && !parent.style().pathEffect())) {
        this->initType(Type::kEmpty);
        *this = parent;
        return;
    }

    SkPathEffect* pe = parent.fStyle.pathEffect();
    SkTLazy<SkPath> tmpPath;
    const GrShape* parentForKey = &parent;
    SkTLazy<GrShape> tmpParent;
    this->initType(Type::kPath);
    fPathData.fGenID = 0;
    if (pe) {
        const SkPath* srcForPathEffect;
        if (parent.fType == Type::kPath) {
            srcForPathEffect = &parent.path();
        } else {
            srcForPathEffect = tmpPath.init();
            parent.asPath(tmpPath.get());
        }
        // Should we consider bounds? Would have to include in key, but it'd be nice to know
        // if the bounds actually modified anything before including in key.
        SkStrokeRec strokeRec = parent.fStyle.strokeRec();
        if (!parent.fStyle.applyPathEffectToPath(&this->path(), &strokeRec, *srcForPathEffect,
                                                 scale)) {
            tmpParent.init(*srcForPathEffect, GrStyle(strokeRec, nullptr));
            *this = tmpParent.get()->applyStyle(apply, scale);
            return;
        }
        // A path effect has access to change the res scale but we aren't expecting it to and it
        // would mess up our key computation.
        SkASSERT(scale == strokeRec.getResScale());
        if (GrStyle::Apply::kPathEffectAndStrokeRec == apply && strokeRec.needToApply()) {
            // The intermediate shape may not be a general path. If we we're just applying
            // the path effect then attemptToReduceFromPath would catch it. This means that
            // when we subsequently applied the remaining strokeRec we would have a non-path
            // parent shape that would be used to determine the the stroked path's key.
            // We detect that case here and change parentForKey to a temporary that represents
            // the simpler shape so that applying both path effect and the strokerec all at
            // once produces the same key.
            tmpParent.init(this->path(), GrStyle(strokeRec, nullptr));
            tmpParent.get()->setInheritedKey(parent, GrStyle::Apply::kPathEffectOnly, scale);
            if (!tmpPath.isValid()) {
                tmpPath.init();
            }
            tmpParent.get()->asPath(tmpPath.get());
            SkStrokeRec::InitStyle fillOrHairline;
            // The parent shape may have simplified away the strokeRec, check for that here.
            if (tmpParent.get()->style().applies()) {
                SkAssertResult(tmpParent.get()->style().applyToPath(&this->path(), &fillOrHairline,
                                                                    *tmpPath.get(), scale));
            } else if (tmpParent.get()->style().isSimpleFill()) {
                fillOrHairline = SkStrokeRec::kFill_InitStyle;
            } else {
                SkASSERT(tmpParent.get()->style().isSimpleHairline());
                fillOrHairline = SkStrokeRec::kHairline_InitStyle;
            }
            fStyle.resetToInitStyle(fillOrHairline);
            parentForKey = tmpParent.get();
        } else {
            fStyle = GrStyle(strokeRec, nullptr);
        }
    } else {
        const SkPath* srcForParentStyle;
        if (parent.fType == Type::kPath) {
            srcForParentStyle = &parent.path();
        } else {
            srcForParentStyle = tmpPath.init();
            parent.asPath(tmpPath.get());
        }
        SkStrokeRec::InitStyle fillOrHairline;
        SkASSERT(parent.fStyle.applies());
        SkASSERT(!parent.fStyle.pathEffect());
        SkAssertResult(parent.fStyle.applyToPath(&this->path(), &fillOrHairline, *srcForParentStyle,
                                                 scale));
        fStyle.resetToInitStyle(fillOrHairline);
    }
    this->attemptToSimplifyPath();
    this->setInheritedKey(*parentForKey, apply, scale);
}

void GrShape::attemptToSimplifyPath() {
    SkRect rect;
    SkRRect rrect;
    SkPath::Direction rrectDir;
    unsigned rrectStart;
    bool inverted = this->path().isInverseFillType();
    SkPoint pts[2];
    if (this->path().isEmpty()) {
        this->changeType(Type::kEmpty);
    } else if (this->path().isLine(pts)) {
        this->changeType(Type::kLine);
        fLineData.fPts[0] = pts[0];
        fLineData.fPts[1] = pts[1];
        fLineData.fInverted = inverted;
    } else if (this->path().isRRect(&rrect, &rrectDir, &rrectStart)) {
        this->changeType(Type::kRRect);
        fRRectData.fRRect = rrect;
        fRRectData.fDir = rrectDir;
        fRRectData.fStart = rrectStart;
        fRRectData.fInverted = inverted;
        // Currently SkPath does not acknowledge that empty, rect, or oval subtypes as rrects.
        SkASSERT(!fRRectData.fRRect.isEmpty());
        SkASSERT(fRRectData.fRRect.getType() != SkRRect::kRect_Type);
        SkASSERT(fRRectData.fRRect.getType() != SkRRect::kOval_Type);
    } else if (this->path().isOval(&rect, &rrectDir, &rrectStart)) {
        this->changeType(Type::kRRect);
        fRRectData.fRRect.setOval(rect);
        fRRectData.fDir = rrectDir;
        fRRectData.fInverted = inverted;
        // convert from oval indexing to rrect indexiing.
        fRRectData.fStart = 2 * rrectStart;
    } else if (SkPathPriv::IsSimpleClosedRect(this->path(), &rect, &rrectDir, &rrectStart)) {
        this->changeType(Type::kRRect);
        // When there is a path effect we restrict rect detection to the narrower API that
        // gives us the starting position. Otherwise, we will retry with the more aggressive
        // isRect().
        fRRectData.fRRect.setRect(rect);
        fRRectData.fInverted = inverted;
        fRRectData.fDir = rrectDir;
        // convert from rect indexing to rrect indexiing.
        fRRectData.fStart = 2 * rrectStart;
    } else if (!this->style().hasPathEffect()) {
        bool closed;
        if (this->path().isRect(&rect, &closed, nullptr)) {
            if (closed || this->style().isSimpleFill()) {
                this->changeType(Type::kRRect);
                fRRectData.fRRect.setRect(rect);
                // Since there is no path effect the dir and start index is immaterial.
                fRRectData.fDir = kDefaultRRectDir;
                fRRectData.fStart = kDefaultRRectStart;
                // There isn't dashing so we will have to preserver inverseness.
                fRRectData.fInverted = inverted;
            }
        }
    }
    if (Type::kPath != fType) {
        fInheritedKey.reset(0);
        if (Type::kRRect == fType) {
            this->attemptToSimplifyRRect();
        } else if (Type::kLine == fType) {
            this->attemptToSimplifyLine();
        }
    } else {
        if (fInheritedKey.count() || this->path().isVolatile()) {
            fPathData.fGenID = 0;
        } else {
            fPathData.fGenID = this->path().getGenerationID();
        }
        if (!this->style().hasNonDashPathEffect()) {
            if (this->style().strokeRec().getStyle() == SkStrokeRec::kStroke_Style ||
                this->style().strokeRec().getStyle() == SkStrokeRec::kHairline_Style) {
                // Stroke styles don't differentiate between winding and even/odd.
                // Moreover, dashing ignores inverseness (skbug.com/5421)
                bool inverse = !this->style().isDashed() && this->path().isInverseFillType();
                if (inverse) {
                    this->path().setFillType(kDefaultPathInverseFillType);
                } else {
                    this->path().setFillType(kDefaultPathFillType);
                }
            } else if (this->path().isConvex()) {
                // There is no distinction between even/odd and non-zero winding count for convex
                // paths.
                if (this->path().isInverseFillType()) {
                    this->path().setFillType(kDefaultPathInverseFillType);
                } else {
                    this->path().setFillType(kDefaultPathFillType);
                }
            }
        }
    }
}

void GrShape::attemptToSimplifyRRect() {
    SkASSERT(Type::kRRect == fType);
    SkASSERT(!fInheritedKey.count());
    if (fRRectData.fRRect.isEmpty()) {
        fType = Type::kEmpty;
        return;
    }
    if (!this->style().hasPathEffect()) {
        fRRectData.fDir = kDefaultRRectDir;
        fRRectData.fStart = kDefaultRRectStart;
    } else if (fStyle.isDashed()) {
        // Dashing ignores the inverseness (currently). skbug.com/5421
        fRRectData.fInverted = false;
    }
    // Turn a stroke-and-filled miter rect into a filled rect. TODO: more rrect stroke shortcuts.
    if (!fStyle.hasPathEffect() &&
        fStyle.strokeRec().getStyle() == SkStrokeRec::kStrokeAndFill_Style &&
        fStyle.strokeRec().getJoin() == SkPaint::kMiter_Join &&
        fStyle.strokeRec().getMiter() >= SK_ScalarSqrt2 &&
        fRRectData.fRRect.isRect()) {
        SkScalar r = fStyle.strokeRec().getWidth() / 2;
        fRRectData.fRRect = SkRRect::MakeRect(fRRectData.fRRect.rect().makeOutset(r, r));
        fStyle = GrStyle::SimpleFill();
    }
}

void GrShape::attemptToSimplifyLine() {
    SkASSERT(Type::kLine == fType);
    SkASSERT(!fInheritedKey.count());
    if (fStyle.isDashed()) {
        // Dashing ignores inverseness.
        fLineData.fInverted = false;
        return;
    } else if (fStyle.hasPathEffect()) {
        return;
    }
    if (fStyle.strokeRec().getStyle() == SkStrokeRec::kStrokeAndFill_Style) {
        // Make stroke + fill be stroke since the fill is empty.
        SkStrokeRec rec = fStyle.strokeRec();
        rec.setStrokeStyle(fStyle.strokeRec().getWidth(), false);
        fStyle = GrStyle(rec, nullptr);
    }
    if (fStyle.isSimpleFill() && !fLineData.fInverted) {
        this->changeType(Type::kEmpty);
        return;
    }
    SkPoint* pts = fLineData.fPts;
    if (fStyle.strokeRec().getStyle() == SkStrokeRec::kStroke_Style) {
        // If it is horizontal or vertical we will turn it into a filled rrect.
        SkRect rect;
        rect.fLeft = SkTMin(pts[0].fX, pts[1].fX);
        rect.fRight = SkTMax(pts[0].fX, pts[1].fX);
        rect.fTop = SkTMin(pts[0].fY, pts[1].fY);
        rect.fBottom = SkTMax(pts[0].fY, pts[1].fY);
        bool eqX = rect.fLeft == rect.fRight;
        bool eqY = rect.fTop == rect.fBottom;
        if (eqX || eqY) {
            SkScalar r = fStyle.strokeRec().getWidth() / 2;
            bool inverted = fLineData.fInverted;
            this->changeType(Type::kRRect);
            switch (fStyle.strokeRec().getCap()) {
                case SkPaint::kButt_Cap:
                    if (eqX && eqY) {
                        this->changeType(Type::kEmpty);
                        return;
                    }
                    if (eqX) {
                        rect.outset(r, 0);
                    } else {
                        rect.outset(0, r);
                    }
                    fRRectData.fRRect = SkRRect::MakeRect(rect);
                    break;
                case SkPaint::kSquare_Cap:
                    rect.outset(r, r);
                    fRRectData.fRRect = SkRRect::MakeRect(rect);
                    break;
                case SkPaint::kRound_Cap:
                    rect.outset(r, r);
                    fRRectData.fRRect = SkRRect::MakeRectXY(rect, r, r);
                    break;
            }
            fRRectData.fInverted = inverted;
            fRRectData.fDir = kDefaultRRectDir;
            fRRectData.fStart = kDefaultRRectStart;
            if (fRRectData.fRRect.isEmpty()) {
                // This can happen when r is very small relative to the rect edges.
                this->changeType(Type::kEmpty);
                return;
            }
            fStyle = GrStyle::SimpleFill();
            return;
        }
    }
    // Only path effects could care about the order of the points. Otherwise canonicalize
    // the point order.
    if (pts[1].fY < pts[0].fY || (pts[1].fY == pts[0].fY && pts[1].fX < pts[0].fX)) {
        SkTSwap(pts[0], pts[1]);
    }
}
