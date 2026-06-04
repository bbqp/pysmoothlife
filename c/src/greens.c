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

enum point_type
{
    PTYPE_INTERIOR = -1,
    PTYPE_BOUNDARY,
    PTYPE_EXTERIOR
};

struct point_s
{
    float t;
    float x, y;
    enum point_type type;
};

int overlap_contains_corners(unsigned cbits)
{
    return CORNERS_NIBBLE(cbits);
}

int set_corner_bits(float x0, float x1, float y0, float y1, float xc, float yc, float r)
{
    int bits;
    float x0mc = x0 - xc;
    float x1mc = x1 - xc;
    float y0mc = y0 - yc;
    float y1mc = y1 - yc;

    float d[4] = {
        x0mc*x0mc + y0mc*y0mc,
        x1mc*x1mc + y0mc*y0mc,
        x1mc*x1mc + y1mc*y1mc,
        x0mc*x0mc + y1mc*y1mc
    };
 
    int inside = 0;
    float r2 = r*r;
    unsigned bits = 0;

    // Determine which of the corners lie within the circle.
    for (int i = 0; i < 4; i++) {
        if (d[i] <= r2) {
            bits |= CORNERS_BIT(i);
        }
    }

    return bits;
}

int set_edge_bits(float x0, float x1, float y0, float y1, float xc, float yc, float r)
{
    int bits = 0;
    int intersects_edge;

    // Determine which of the edges overlaps with the circle.
    int curr, next;
    float x2, y2;
    float x02 = x0*x0;
    float x12 = x1*x1;
    float y02 = y0*y0;
    float y12 = y1*y1;

    float d[4] = {
        x0mc*x0mc + y0mc*y0mc,
        x1mc*x1mc + y0mc*y0mc,
        x1mc*x1mc + y1mc*y1mc,
        x0mc*x0mc + y1mc*y1mc
    };

    float r2 = r*r;

    float minl = d[0] < d[1] ? d[0] : d[1];
    float minr = d[2] < d[3] ? d[2] : d[3];
    float mind = minl < minr ? minl : minr;

    if (mind >= r2) {
        return bits;
    }

    for (int i = 0; i < 4; i++) {
        switch (i) {
            case 0:
                x2 = r2 - y0*y0; 
                intersects_edge = x0*x0 < x2 && x2 < x1*x1;
                break;
            case 1:
                y2 = r2 - x1*x1; 
                intersects_edge = y0*y0 < y2 && y2 < y1*y1;
                break;
            case 2:
                x2 = r2 - y1*y1; 
                intersects_edge = x0*x0 < x2 && x2 < x1*x1;
                break;
            case 3:
                y2 = r2 - x0*x0; 
                intersects_edge = y0*y0 < y2 && y2 < y1*y1;
        }

        if (intersects_edge) { 
            bits |= EDGES_BIT(i);
        }
    }

    return bits;
}


int set_quadrant_bits(float x0, float x1, float y0, float y1, float xc, float yc, float r)
{
    int bits = 0;

    float x0mc = x0 - xc;
    float x1mc = x1 - xc;
    float y0mc = y0 - yc;
    float y1mc = y1 - yc;

    // Check the quadrants in which the shape is located.
    if (x0mc >= 0) {
        if (y0mc >= 0) {
            bits |= QUADRANT_BIT(0);
        } else {
            bits |= QUADRANT_BIT(3);
        }

        if (y1mc >= 0) {
            bits |= QUADRANT_BIT(0);
        } else {
            bits |= QUADRANT_BIT(3);
        }
    } else {
        if (y0mc >= 0) {
            bits |= QUADRANT_BIT(1);
        } else {
            bits |= QUADRANT_BIT(2);
        }

        if (y1mc >= 0) {
            bits |= QUADRANT_BIT(1);
        } else {
            bits |= QUADRANT_BIT(2);
        }
    }

    if (x1mc >= 0) {
        if (y0mc >= 0) {
            bits |= QUADRANT_BIT(0);
        } else {
            bits |= QUADRANT_BIT(3);
        }

        if (y1mc >= 0) {
            bits |= QUADRANT_BIT(0);
        } else {
            bits |= QUADRANT_BIT(3);
        }
    } else {
        if (y0mc >= 0) {
            bits |= QUADRANT_BIT(1);
        } else {
            bits |= QUADRANT_BIT(2);
        }

        if (y1mc >= 0) {
            bits |= QUADRANT_BIT(1);
        } else {
            bits |= QUADRANT_BIT(2);
        }
    }

    return bits;
}

int get_quadrant(float x, float y)
{
    int quadrant;

    if (x >= 0) {
        quadrant = y >= 0 ? 0 : 3;
    } else {
        quadrant = y >= 0 ? 1 : 2;
    }

    return quadrant;
}


// Prune the queue of exterior points.
void prune_point_queue(struct point_s *pqueue, int *head, int *tail)
{
    int start = *head;
    int end = *tail;    
    int length = end - start + 1;

    while (pqueue[start].type == PTYPE_EXTERIOR && start <= end) {
        start++;
        length--;
    }

    if (start <= end) {
        // Move everything between start and end to the beginning.
        memmove(pqueue, pqueue + start, length * sizeof(struct pqueue));

        int i, j;
        struct point_s temp;

        for (i = start + 1; i <= end; i++) {
            // S
            if (pqueue[start].type == PTYPE_EXTERIOR) {
                for (j = i + 1; j <= end; j++) {
                    if (pqueue[j].type = PTYPE_EXTERIOR) {
                        temp = pqueue[i];
                        pqueue[i] = pqueue[j];
                        pqueue[j] = temp;

                        length--;
                        break;                       
                    }
                }
            }
        }

        end = start + length - 1;
    }

    *head = start;
    *tail = end;
}

// For now, we'll assume that sqrt((x1-x0)**2 + (y1-y0)*2) <= 2*r
// This means that if all points in a square are exterior points, then
// the entire square has to be outside of the circle.
float fgreens(float x0, float x1, float y0, float y1, float xc, float yc, float r)
{
    float area = 0;
    int cbits = set_corner_bits(x0, x1, y0, y1, xc, yc, r); 
    int ebits = set_edge_bits(x0, x1, y0, y1, xc, yc, r); 
    int qbits = set_quadrant_bits(x0, x1, y0, y1, xc, yc, r);
    int bits = cbits | ebits | qbits;

    int head = 0, tail = 0;
    struct point_s pqueue[12];
    float lsq;
    float r2 = r*r;
    float t[2];

    float x0c = x0 - xc;
    float y0c = y0 - yc;
    float x1c = x1 - xc;
    float y1c = y1 - yc;
    float d[4] = {
        x0c*x0c + y0c*y0c,
        x1c*x1c + y0c*y0c,
        x0c*x0c + y1c*y1c,
        x1c*x1c + y1c*y1c,
    };

    int all_interior = d[0] <= r2 && d[1] <= r2 && d[2] <= r2 && d[3] <= r2;
    int all_exterior = d[0] >= r2 && d[1] >= r2 && d[2] >= r2 && d[3] >= r2;

    if (all_interior) {
        return (x1 - x0) * (y1 - y0);
    }

    if (all_exterior) {
        return 0;
    }

    // Add the points moving counterclockwise, starting from the lower left corner.
    pqueue[tail].x = x0;
    pqueue[tail].y = y0;

    lsq = x0 * x0 + y0 * y0;
    absdiff = lsq < r2 ? r2 - lsq : lsq - r2;

    if (lsq < r2) {
        pqueue[tail].type = PTYPE_INTERIOR;
    } else (lsq > r2) {
        pqueue[tail].type = PTYPE_EXTERIOR;
    } else if (absdiff <= 10 * FLT_EPSILON) {
        pqueue[tail].type = PTYPE_BOUNDARY;
    }

    // Next, find the intersections along the bottom face.
    xa[0] = -sqrt(r2 - y0*y0);
    xa[1] = -xa[0];
    ya[0] = ya[1] = y0;

    for (int i = 0; i < 2; i++) {
        if (x0 <= xa[i] && xa[i] <= x1) {
            tail++;

            pqueue[tail].x = xa[i];
            pqueue[tail].y = ya[i];
            pqueue[tail].type = PTYPE_BOUNDARY;
        }
    }

    // Add the next set of points starting from the lower-right corner and along the right face.
    pqueue[tail].x = x1;
    pqueue[tail].y = y0;

    lsq = x1 * x1 + y0 * y0;
    absdiff = lsq < r2 ? r2 - lsq : lsq - r2;

    if (lsq < r2) {
        pqueue[tail].type = PTYPE_INTERIOR;
    } else (lsq > r2) {
        pqueue[tail].type = PTYPE_EXTERIOR;
    } else if (absdiff <= 10 * FLT_EPSILON) {
        pqueue[tail].type = PTYPE_BOUNDARY;
    }

    // Next, find the intersections along the right face.
    xa[0] = xa[1] = x1;
    ya[0] = -sqrt(r2 - x1*x1);
    ya[1] = -ya[0];

    for (int i = 0; i < 2; i++) {
        if (y0 <= ya[i] && ya[i] <= y1) {
            tail++;

            pqueue[tail].x = xa[i];
            pqueue[tail].y = ya[i];
            pqueue[tail].type = PTYPE_BOUNDARY;
        }
    }

    // Add the next set of points starting from the upper-right corner and along the top face.
    pqueue[tail].x = x1;
    pqueue[tail].y = y1;

    lsq = x1 * x1 + y1 * y1;
    absdiff = lsq < r2 ? r2 - lsq : lsq - r2;

    if (lsq < r2) {
        pqueue[tail].type = PTYPE_INTERIOR;
    } else (lsq > r2) {
        pqueue[tail].type = PTYPE_EXTERIOR;
    } else if (absdiff <= 10 * FLT_EPSILON) {
        pqueue[tail].type = PTYPE_BOUNDARY;
    }

    // Next, find the intersections along the bottom face.
    xa[0] = sqrt(r2 - y1*y1);
    xa[1] = -xa[0];
    ya[0] = ya[1] = y1;

    for (int i = 1; i >= 0; i--) {
        if (x0 <= xa[i] && xa[i] <= x1) {
            tail++;

            pqueue[tail].x = xa[i];
            pqueue[tail].y = ya[i];
            pqueue[tail].type = PTYPE_BOUNDARY;
        }
    }

    // Add the next set of points starting from the lower-right corner and along the right face.
    pqueue[tail].x = x0;
    pqueue[tail].y = y1;

    lsq = x0 * x0 + y0 * y0;
    absdiff = lsq < r2 ? r2 - lsq : lsq - r2;

    if (lsq < r2) {
        pqueue[tail].type = PTYPE_INTERIOR;
    } else (lsq > r2) {
        pqueue[tail].type = PTYPE_EXTERIOR;
    } else if (absdiff <= 10 * FLT_EPSILON) {
        pqueue[tail].type = PTYPE_BOUNDARY;
    }

    // Next, find the intersections along the right face.
    xa[0] = xa[1] = x0;
    ya[0] = sqrt(r2 - x0*x0);
    ya[1] = -ya[0];

    for (int i = 1; i >= 0; i--) {
        if (y0 <= ya[i] && ya[i] <= y1) {
            tail++;

            pqueue[tail].x = xa[i];
            pqueue[tail].y = ya[i];
            pqueue[tail].type = PTYPE_BOUNDARY;
        }
    }

    // Add the first point to the end.
    pqueue[++tail] = (struct point_s)pqueue[0];

    // Prune the queue.
    prune_point_queue(pqueue, &head, &tail);

    // Compute the contributions to the sum using green's theorem.
    struct point_s curr, next;
    float curve_sum = 0;

    while (head <= tail) {
        // Set the boundary integral for the curve/segment to be zero.
        curve_sum = 0;

        // Get the current and next point.
        curr = point[head];
        next = point[head + 1];

        if (curr.type == PTYPE_INTERIOR) {
            if (curr.y == next.y) {
                curve_sum = -0.5 * curr.y * (next.x - curr.x);
            } else if (curr.x == next.x) {
                curve_sum = 0.5 * curr.x * (next.y - curr.y);
            }               
        } else if (curr.type == PTYPE_BOUNDARY) {
            if (next.type == PTYPE_INTERIOR) {
                if (curr.y == next.y) {
                    curve_sum = -0.5 * curr.y * (next.x - curr.x);
                } else if (curr.x == next.x) {
                    curve_sum = 0.5 * curr.x * (next.y - curr.y);
                }                
            } else if (next.type == PTYPE_BOUNDARY) {
                curve_sum = 0.5 * r * r * (next.t - curr.t);
            }
        }

        // Add the contribution to the area.
        area += curve_sum;

        head++;
    }

    return area;
}
