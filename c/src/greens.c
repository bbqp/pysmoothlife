#define CORNERS_MASK 0x0000ff
#define CORNERS_BIT(i) (((int)1) << i)
#define GET_CORNERS_NIBBLE(n) ((n) & CORNERS_MASK)
#define SET_CORNERS_NIBBLE(n, bits) ((n) |= bits)

#define EDGES_MASK 0x00ff00
#define EDGES_BIT(i) (((int)1) << (i + 4))
#define GET_EDGES_NIBBLE(n) (((n) & EDGES_MASK) >> 4)
#define SET_EDGES_NIBBLE(n, bits) ((n) |= ((bits) << 4))

#define QUADRANTS_MASK 0xff0000
#define QUADRANTS_BIT(i) (((int)1) << (i + 8))
#define GET_QUADRANTS_NIBBLE(n) (((n) & QUADRANTS_MASK) >> 8)
#define SET_QUADRANTS_NIBBLE(n, bits) ((n) |= ((bits) << 8))

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

float fgreens(float x0, float x1, float y0, float y1, float xc, float yc, float r)
{
    float area = 0;

    int cbits = 
    
    return area;
}
