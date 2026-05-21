#include <math.h>


static float circle(float t, const float *radius, int n)
{
    return 
}

struct curve_s
{
    float t0, tn;
    float (*curve)(float, const float *, int, float *, float *);
}

struct region_s
{
    float *endpoints;
    int num_endpoints;
    float min_x, min_y;
    float max_x, max_y;
}


float intersection(float x0, float x1, float y0, float y1, float cx, float cy, float r)
{
    float xstack[12], ystack[12];
    int head, tail;
}


int is_point_inside(const struct region_s *region, float x, float y)
{
    int inside = 0;

    float min_x = region->min_x;
    float max_x = region->max_x;
    float min_y = region->min_y;
    float max_y = region->max_y;

    if (min_x <= x && x <= max_x && min_y <= y && y <= max_y) {
        
    }

    return inside;
}

float solve_on_boundary(const struct region_s *region, float x, float y)
{
    float tvalue = -1;
    int boundary_region = -1;

    // Determine if we are inside or outside.
    int inside = 0, outside = 0;


    
    return tvalue;
}

int inside(const struct region_s *region, float, float)
{
    // Find the value(s) of t for which (x,y) is on the boundary.
    float t = solve_on_boundary(boundary, x, y);
}
