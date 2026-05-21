#define CORNERS_MASK 0x0000ff
#define CORNERS_BIT(i) (((int)1) << i)
#define CORNERS_NIBBLE(n) ((n) & CORNERS_MASK)

#define EDGES_MASK 0x00ff00
#define EDGES_BIT(i) (((int)1) << (i + 4))
#define EDGES_NIBBLE(n) ((n) & EDGES_MASK)

#define QUADRANTS_MASK 0xff0000
#define QUADRANTS_BIT(i) (((int)1) << (i + 8))
#define QUADRANTS_NIBBLE(n) ((n) & QUADRANTS_MASK)

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
    // Determine which of the edges overlaps with the circle.
    int curr, next;
    float xa, ya;

    float d[4] = {
        x0mc*x0mc + y0mc*y0mc,
        x1mc*x1mc + y0mc*y0mc,
        x1mc*x1mc + y1mc*y1mc,
        x0mc*x0mc + y1mc*y1mc
    };
 
    for (int i = 0; i < 4; i++) {
        curr = i;
        next = (i + 1) % 4;

        if (d[curr] >= r2) {
            
        }

        if (d[curr] <= r2 && d[next] >= r2) {
            switch (i) {
                case 0:
                    xa = r2 - y0*y0; 

                    if (xa > x0*x0 && xa <= x1*x1) {
                        cbits |= EDGES_BIT(i);
                    }
                    break;
                case 1:
                    ya = r2 - x1*x1; 

                    if (ya > y0*y0 && ya <= y1*y1) {
                        cbits |= EDGES_BIT(i);
                    }
                    break;
                case 2:
                    xa = r2 - y1*y1; 

                    if (xa >= x0*x0 && xa < x1*x1) {
                        cbits |= EDGES_BIT(i);
                    }
                    break;
                case 3:
                    ya = r2 - x1*x1; 

                    if (ya >= y0*y0 && ya < y1*y1) {
                        cbits |= EDGES_BIT(i);
                    }
                    break;
            }
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

float fgreens(float x0, float x1, float y0, float y1, float xc, float yc, float r)
{
    float area = 0;
    
    if (cbits) {
        int quadrant = QUADRANTS_NIBBLE(cbits, QUADRANTS_MASK);

        switch (quadrant) {
            case QUADRANT1:
                float xa, xb;
                float ya, yb;

                float ta = acosf(y0 / r);
                float tb = asinf(x0 / r);

                xa = r * cosf(ta);
                ya = r * sinf(ta);
                xb = r * cosf(tb);
                yb = r * sinf(tb);

                area = 0.5 * r * r * (tb - ta)
                     + 0.5 * x0 * (yb - y0)
                     - 0.5 * y0 * (xb - x0);
            break;

            case QUADRANT2:
            break;

            case QUADRANT3:
            break;

            case QUADRANT4:
            break;
        }
    }

    return area;
}
