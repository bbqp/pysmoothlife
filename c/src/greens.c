#define CORNERS_MASK 0x0000ff
#define SET_CORNERS_BIT(bits, i) ((bits) |= (((int)1) << (i)))
#define GET_CORNERS_BIT(i) (((bits) >> (i)) & ((int)1))
#define GET_CORNERS_NIBBLE(n) ((n) & CORNERS_MASK)
#define SET_CORNERS_NIBBLE(n, bits) ((n) |= (bits))

#define EDGES_MASK 0x00ff00
#define SET_EDGES_BIT(i) ((bits) |= (((int)1) << ((i) + 4))
#define GET_EDGES_BIT(i) (((bits) >> ((i) + 4)) & ((int)1))
#define GET_EDGES_NIBBLE(n) (((n) & EDGES_MASK) >> 4)
#define SET_EDGES_NIBBLE(n, bits) ((n) |= ((bits) << 4))

#define QUADRANTS_MASK 0xff0000
#define SET_QUADRANTS_BIT(i) ((bits) |= (((int)1) << ((i) + 8))
#define GET_QUADRANTS_BIT(i) (((bits) >> ((i) + 8)) & ((int)1))
#define GET_QUADRANTS_NIBBLE(n) (((n) & QUADRANTS_MASK) >> 8)
#define SET_QUADRANTS_NIBBLE(n, bits) ((n) |= ((bits) << 8))

#define F_PI atan2f(0.0f, -1.0f)

#define NIBBLE_SUM(n, mask) ( \
    (((n & mask) >> 0) & 1) + \
    (((n & mask) >> 1) & 1) + \
    (((n & mask) >> 2) & 1) + \
    (((n & mask) >> 3) & 1) \
)

enum quadrant {
    QUADRANT1 = 1,
    QUADRANT2 = 2,
    QUADRANT3 = 4,
    QUADRANT4 = 8
};

enum axis {
    AXIS1 = 1,
    AXIS2 = 1,
    AXIS3 = 1,
    AXIS4 = 1,
};

enum point_type
{
    PTYPE_INTERIOR = 1,
    PTYPE_BOUNDARY = 2,
    PTYPE_EXTERIOR = 4
};

enum location_type
{
    LTYPE_CORNER = 1,
    LTYPE_FACE = 2
};

enum segment_type
{
    STYPE_NONE = 0,
    STYPE_LINE = 1,
    STYPE_CIRCLE = 2
};

struct point_s
{
    float t;
    float x, y;
    enum point_type ptype;
    enum location_type ltype;
    enum segment_type stype;
    int location_index;

    struct point_s *next, *prev;
};

struct segment_s
{
    struct point_s start, end;
    enum segment_type type;
};

struct point_list_s
{
    int length;
    struct point_s *head, *tail;
};

struct point_s *point_alloc(float t, float x, float y, enum point_type ptype, enum location_type ltype, int lindex)
{
    struct point_s *point = NULL;

    point = calloc(1, sizeof(struct point_s));

    if (point != NULL) {
        point->t = t;
        point->x = x;
        point->y = y;
        point->ptype = ptype;
        point->ltype;
        point->stype = STYPE_NONE;
        point->location_index = lindex;

        point->prev = point->next = NULL;
    }

    return point;
}

// Assumes point->next != NULL
void point_set_segment_type(struct point_s *p1, const struct point_s *p2)
{
    point_type p1_ptype = p1->ptype;
    location_type p1_ltype = p1->ltype;
    int p1_lindex = p1->location_index;

    point_type p2_ptype = p2->ptype;
    location_type p2_ltype = p2->ltype;
    int p2_lindex = p2->location_index;

    if (p1->ptype == p2->ptype) {
        if (p1->ptype == PTYPE_INTERIOR) {
            p1->stype = STYPE_LINE;
        } else {
            // For boundary points, we need to be more careful. At the given
            // and next point, we'll use information about the derivative of
            // the circular countour at the intersection points to determine
            // if the circular contour is inside the square region or outside
            // of it. If inside, we integrate along the circular contour.
            // Otherwise, we integrate across the line segment along a given
            // face.

            if (p1->location_index == p2->location_index) {
                // This means both points are on the same face, so we just need
                // to check if the circle is in or outside of the rectangle. We
                // can do this by checking how the x-coordinates change going
                // from point p1 to point p2.

                switch (p1->location_index) {
                    case 0:
                        if (p1->x < p2->x) {
                            p1->stype = STYPE_LINE;
                        } else {
                            p1->stype = STYPE_CIRCLE;
                        }

                        break;
                    case 1:
                        if (p1->y < p2->y) {
                            p1->stype = STYPE_LINE;
                        } else {
                            p1->stype = STYPE_CIRCLE;
                        }

                        break;
                    case 2:
                        if (p1->x > p2->x) {
                            p1->stype = STYPE_LINE;
                        } else {
                            p1->stype = STYPE_CIRCLE;
                        }

                        break;
                    case 3:
                        if (p1->y > p2->y) {
                            p1->stype = STYPE_LINE;
                        } else {
                            p1->stype = STYPE_CIRCLE;
                        }
                }
            } else {
                // In this case, we're not indexed to the same face, but we
                // could, for example, have intersection points at corners
                // where the x- or y-coordinates are equal, meaning the
                // starting and ending points lie on the same horizontal or
                // vertical surface of the rectangle. Similarly, one point
                // could cross through a face and another through a corner.

                int adjacent_locations = 0;
                switch (p1->location_index) {
                    case 0:
                        switch (p2->location_index) {
                            case 1:
                            case 3:
                                adjacent_locations = 1;
                        }
                        break;
                  case 1:
                        switch (p2->location_index) {
                            case 0:
                            case 2:
                                adjacent_locations = 1;
                        }
                        break;
                  case 2:
                        switch (p2->location_index) {
                            case 1:
                            case 3:
                                adjacent_locations = 1;
                        }
                        break;
                  case 3:
                        switch (p2->location_index) {
                            case 0:
                            case 2:
                                adjacent_locations = 1;
                        }
                }

                if (p1->ltype == LTYPE_FACE && p2->ltype == LTYPE_FACE) {
                    p1->stype = STYPE_CIRCLE;
                } else {
                    int both_corners = p1->ltype == p2->ltype && p1->ltype == LTYPE_CORNER;
                    int ends_at_corner = both_corners || p2->ltype == LTYPE_CORNER;

                    if (adjacent_locations) {
                        if (ends_at_corner) {
                            switch (p1->location_index) {
                                case 0:
                                    if (p1->x < p2->x) {
                                        p1->stype = STYPE_LINE;
                                    } else {
                                        p1->stype = STYPE_CIRCLE;
                                    }

                                    break;
                                case 1:
                                    if (p1->y < p2->y) {
                                        p1->stype = STYPE_LINE;
                                    } else {
                                        p1->stype = STYPE_CIRCLE;
                                    }

                                    break;
                                case 2:
                                    if (p1->x > p2->x) {
                                        p1->stype = STYPE_LINE;
                                    } else {
                                        p1->stype = STYPE_CIRCLE;
                                    }

                                    break;
                                case 3:
                                    if (p1->y > p2->y) {
                                        p1->stype = STYPE_LINE;
                                    } else {
                                        p1->stype = STYPE_CIRCLE;
                                    }
                            }
                        } else {
                            p1->stype = STYPE_CIRCLE;
                        }
                    } else {
                        p1->stype = STYPE_CIRCLE;
                    }  // else point locations are not on adjacent portions of the boundary.
                }  // else one or both points are corners
            }  // else we are not indexed to the same location
        }  // else we gave two boundary intersection points.
    } else {
        p1->stype = STYPE_LINE;
    }
}

void point_free(struct point_s *point)
{
    free(point);
}

void point_list_init(struct point_list_s *plist, const struct point_s *head)
{
    plist->head = head;
    plist->length = 1;
}

int point_list_add(struct point_list_s *plist, const struct point_s *point)
{
    plist->tail->next = point;
    plist->tail = point;
    plist->length++;

    return plist->length;
}

int point_list_remove_next(struct point_list_s *plist, struct point_s *point)
{
    struct point_s *removed = NULL;
    int removal_status = 0;

    if (point->next != NULL) {
        if (point->next != plist->tail) {
            removed = point->next;
            removed->next->prev = point;
            point->next = removed->next;
            removed->next = removed->prev = NULL;            
        } else {
            removed = point->next;
            point->next = removed->next;
            removed->prev = NULL; 
        }

        free(removed);

        removal_status = plist->length--;
    }

    return removal_status;
}

void point_list_clear(struct point_list_s *plist)
{
    struct point_s *curr = plist->head;

    while (curr->next != NULL) {
        plist->head = plist->head->next;
        curr->prev = curr->next = plist->head->prev = NULL;

        free(curr);
        curr = plist->head;
    }

    curr->next = curr->prev = NULL;
    free(curr);

    plist->head = plist->tail = NULL;
    plist->length = 0;
}

void get_quadrant_and_axis(float x, float y, int *quadrant, int *axis)
{
    int quad = 0;
    int axis = -1;

    if (x > 0) {
        if (y > 0) {
            quad = (1 << 0);
        } else if (y < 0) {
            quad = (1 << 3);
        } else {
            quad = (1 << 0) | (1 << 3);
            ax = (1 << 0);
        }
    } else if (x < 0) {
        if (y > 0) {
            quad = (1 << 1);
        } else if (y < 0) {
            quad = (1 << 2);
        } else {
            quad = (1 << 1) | (1 << 2);
            ax = (1 << 2);
        }
    } else {
        if (y > 0) {
            quad = (1 << 0) | (1 << 1);
            ax = (1 << 1);
        } else if (y < 0) {
            quad = (1 << 2) | (1 << 3);
            ax = (1 << 3);
        } else {
            quad = (1 << 0) | (1 << 1) | (1 << 2) | (1 << 3);
            axis = (1 << 0) | (1 << 1) | (1 << 2) | (1 << 3);
        }
    }

    *quadrant = quad;
    *axix = ax;
}

void point_list_build(struct point_list_s *plist, float x0, float y0, float x1, float y1, float xc, float yc, float r)
{
    // Mathematical constant pi.
    const float f_pi = atan2f(0.0f, -1.0f);
    const float f_pi2 = f_pi / 2;
    const float f_2pi = 2 * f_pi;

    // Tolerance for corner proximity.
    const float ctol = 10 * FLT_EPSILON;

    // Intersection points.
    float ti[2];
    float xi[2];
    float yi[2];

    // Shorthands.
    float r2 = r * r;

    // Point (vector) length.
    float mag2;
    float absdiff;

    // Placeholder struct for new point.
    struct point_s *new_point = NULL;

    // Placeholders for point and location types.
    enum point_type ptype;
    enum location_type ltype;
    int location_index;

    // Shift the points to the origin.
    x0 -= xc;
    y0 -= yc;
    x1 -= xc;
    y1 -= yc;

    //------------------------------------------------------------------------
    // Add the points from face 0.
    //------------------------------------------------------------------------

    location_index = 0;
    mag2 = x0 * x0 + y0 * y0;
    mag = sqrtf(mag2);
    absdiff = mag > r ? mag - r : r - mag;

    if (absdiff < ctol) {
        ptype = PTYPE_BOUNDARY;
    } else if (mag2 < r2) {
        ptype = PTYPE_INTERIOR;
    } else if (mag2 > r2) {
        ptype = PTYPE_EXTERIOR;
    }

    ltype = LTYPE_CORNER;

    // Add the corner point.
    new_point = point_alloc(0, x0, y0, ptype, ltype, location_index);
    point_list_add(plist, new_point);

    // Now find the intersection points on face 0.
    ti[0] = asinf(y0 / r);
    if (ti[0] < 0) {
        ti[1] = f_2pi + ti[0];
        ti[0] = f_pi - ti[0];
    } else {
        ti[1] = f_pi - ti[0];
    }

    yi[0] = yi[1] = y0;
    xi[0] = r * cosf(ti[0]);
    xi[1] = r * cosf(ti[1]);

    // For now, just set the point type since it's known by construction.
    // We'll set the location type depending on the proximity to a corner.
    ptype = PTYPE_BOUNDARY;

    for (int k = 0; k < 2; k++) {
        if (x0 <= xi[k] && xi[k] <= x1) {
            new_point = point_alloc(ti[k], xi[k], yi[k], ptype, ltype, location_index);

            // Now check if this point is on a proper face or boundary.
            float dx = xi[k] - x0;
            if (dx <= ctol) {
                new_point->ltype = LTYPE_CORNER;
                new_point->location_index = location_index;
            } else {
                dx = x1 - xi[k];

                if (dx <= ctol) {
                    new_point->ltype = LTYPE_CORNER;
                    new_point->location_index = location_index + 1;
                }
            }

            // Add the point to the list.
            point_list_add(plist, new_point);
        }
    }

    //------------------------------------------------------------------------
    // Add the points from face 1.
    //------------------------------------------------------------------------

    location_index = 1;
    mag2 = x1 * x1 + y0 * y0;
    mag = sqrtf(mag2);
    absdiff = mag > r ? mag - r : r - mag;

    if (absdiff < ctol) {
        ptype = PTYPE_BOUNDARY;
    } else if (mag2 < r2) {
        ptype = PTYPE_INTERIOR;
    } else if (mag2 > r2) {
        ptype = PTYPE_EXTERIOR;
    }

    ltype = LTYPE_CORNER;

    // Add the corner point.
    new_point = point_alloc(0, x1, y0, ptype, ltype, location_index);
    point_list_add(plist, new_point);

    // Now find the intersection points on face 1.
    ti[0] = acosf(x1 / r);
    if (ti[0] < f_pi2) {
        ti[1] = t[0];
        ti[0] = f_2pi - ti[0];
    } else {
        ti[1] = f_2pi - ti[0];
    }

    xi[0] = xi[1] = x1;
    yi[0] = r * sinf(ti[0]);
    yi[1] = r * sinf(ti[1]);

    // For now, just set the point type since it's known by construction.
    // We'll set the location type depending on the proximity to a corner.
    ptype = PTYPE_BOUNDARY;

    for (int k = 0; k < 2; k++) {
        if (y0 <= yi[k] && yi[k] <= y1) {
            new_point = point_alloc(ti[k], xi[k], yi[k], ptype, ltype, location_index);

            // Now check if this point is on a proper face or boundary.
            float dy = yi[k] - y0;
            if (dy <= ctol) {
                new_point->ltype = LTYPE_CORNER;
                new_point->location_index = location_index;
            } else {
                dy = y1 - yi[k];

                if (dy <= ctol) {
                    new_point->ltype = LTYPE_CORNER;
                    new_point->location_index = location_index + 1;
                }
            }

            // Add the point to the list.
            point_list_add(plist, new_point);
        }
    }

    //------------------------------------------------------------------------
    // Add the points from face 2.
    //------------------------------------------------------------------------

    location_index = 2;
    mag2 = x1 * x1 + y1 * y1;
    mag = sqrtf(mag2);
    absdiff = mag > r ? mag - r : r - mag;

    if (absdiff < ctol) {
        ptype = PTYPE_BOUNDARY;
    } else if (mag2 < r2) {
        ptype = PTYPE_INTERIOR;
    } else if (mag2 > r2) {
        ptype = PTYPE_EXTERIOR;
    }

    ltype = LTYPE_CORNER;

    // Add the corner point.
    new_point = point_alloc(0, x1, y1, ptype, ltype, location_index);
    point_list_add(plist, new_point);

    // Now find the intersection points on face 2.
    ti[0] = asinf(y1 / r);
    if (ti[0] < 0) {
        ti[1] = f_2pi + ti[0];
        ti[0] = f_pi - ti[0];
    } else {
        ti[1] = f_pi - ti[0];
    }

    yi[0] = yi[1] = y1;
    xi[0] = r * cosf(ti[0]);
    xi[1] = r * cosf(ti[1]);

    // For now, just set the point type since it's known by construction.
    // We'll set the location type depending on the proximity to a corner.
    ptype = PTYPE_BOUNDARY;

    for (int k = 0; k < 2; k++) {
        if (x0 <= xi[k] && xi[k] <= x1) {
            new_point = point_alloc(ti[k], xi[k], yi[k], ptype, ltype, location_index);

            // Now check if this point is on a proper face or boundary.
            float dx = xi[k] - x0;
            if (dx <= ctol) {
                new_point->ltype = LTYPE_CORNER;
                new_point->location_index = location_index;
            } else {
                dx = x1 - xi[k];

                if (dx <= ctol) {
                    new_point->ltype = LTYPE_CORNER;
                    new_point->location_index = location_index + 1;
                }
            }

            // Add the point to the list.
            point_list_add(plist, new_point);
        }
    }

    //------------------------------------------------------------------------
    // Add the points from face 3.
    //------------------------------------------------------------------------

    location_index = 3;
    mag2 = x0 * x0 + y1 * y1;
    mag = sqrtf(mag2);
    absdiff = mag > r ? mag - r : r - mag;

    if (absdiff < ctol) {
        ptype = PTYPE_BOUNDARY;
    } else if (mag2 < r2) {
        ptype = PTYPE_INTERIOR;
    } else if (mag2 > r2) {
        ptype = PTYPE_EXTERIOR;
    }

    ltype = LTYPE_CORNER;

    // Add the corner point.
    new_point = point_alloc(0, x0, y1, ptype, ltype, location_index);
    point_list_add(plist, new_point);

    // Now find the intersection points on face 1.
    ti[0] = acosf(x0 / r);
    if (ti[0] < f_pi2) {
        ti[1] = t[0];
        ti[0] = f_2pi - ti[0];
    } else {
        ti[1] = f_2pi - ti[0];
    }

    xi[0] = xi[1] = x0;
    yi[0] = r * sinf(ti[0]);
    yi[1] = r * sinf(ti[1]);

    // For now, just set the point type since it's known by construction.
    // We'll set the location type depending on the proximity to a corner.
    ptype = PTYPE_BOUNDARY;

    for (int k = 0; k < 2; k++) {
        if (y0 <= yi[k] && yi[k] <= y1) {
            new_point = point_alloc(ti[k], xi[k], yi[k], ptype, ltype, location_index);

            // Now check if this point is on a proper face or boundary.
            float dy = yi[k] - y0;
            if (dy <= ctol) {
                new_point->ltype = LTYPE_CORNER;
                new_point->location_index = location_index;
            } else {
                dy = y1 - yi[k];

                if (dy <= ctol) {
                    new_point->ltype = LTYPE_CORNER;
                    new_point->location_index = location_index + 1;
                }
            }

            // Add the point to the list.
            point_list_add(plist, new_point);
        }
    }
}

void point_list_prune(struct point_list_s *plist)
{
    struct point_s *curr = plist->head;

    while (curr != plist->tail) {
        if (curr->ptype = PTYPE_EXTERIOR) {
            before = curr->prev;
            after = curr->next;

            before->next = curr->next;
            after->prev = curr->prev;

            curr->prev = curr->next = NULL;
            free(curr);
        }

        curr = curr->next;
    }

    if (curr->ptype == PTYPE_EXTERIOR) {
        before = curr->prev;
        before->next = curr->next;
        plist->tail = before;

        curr->prev = curr->next = NULL;
        free(curr);
    }
}

void point_list_set_segment_types(struct point_list_s *plist)
{
    struct point_s *curr = plist->head;
    struct point_s *next = NULL;

    while (curr != NULL) {
        next = curr->next != NULL ? curr->next : plist->head;
        point_set_segment_type(curr, next);
    }
}

//
// void point_list_compute_greens()
//
// Uses Green's theorem to compute the area of a region comprising the overlap
// of a circle of radius r with the rectangle [x0, x1] x [y0, y1].
void point_list_compute_greens(const struct point_list_s *plist, float r)
{
    float contour_sums = 0.0f;
    struct point_s *curr, *p1, *p2;
    float r2 = r*r;

    int start_face;

    float p1t;
    float p2t;

    float p1x;
    float p2x;

    float p1y;
    float p2y;

    for (curr = plist->head; curr != NULL; curr = curr->next) {
        // Grab the points connected by either a line or circular curve.
        p1 = curr;
        p2 = curr->next != NULL ? curr->next : plist->head;

        // Extract the coordinates, type, and location of the first point.
        p1t = p1->t;
        p1x = p1->x;
        p1y = p1->y;

        // Extract the coordinates, type, and location of the second point.
        p2t = p2->t;
        p2x = p2->x;
        p2y = p2->y;

        if (p1->stype == STYPE_CIRCLE) {
            contour_sums += 0.5 * r2 * (p2t - p1t);
        } else {
            switch (p1->location_index) {
                case 0:
                case 2:
                    contour_sums -= 0.5 * p1y * (p2x - p1x);
                    break;
                case 1:
                case 3:
                    contour_sums += 0.5 * p1x * (p2y - p1y);
            }
        }
    }

    return countour_sums;
}

