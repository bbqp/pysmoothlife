#ifndef GRID_H
#define GRID_H

#include "transitionfunction.h"

#define CIRCLE_AREA(r) (M_PI * ((r) * (r)))
#define ANNULUS_AREA(ri, ro) (M_PI * ((ro) * (ro) - (ri) * (ri)))

enum grid_region {
	GR_OUTSIDE = -1,
	GR_CIRCLE,
	GR_INNER_ALIAS,
	GR_ANNULUS,
	GR_OUTER_ALIAS
};

struct grid_s
{
	float *data;
	float *distsq_closest;
	float *distsq_furthest;
	float *dist_midpoint;
	float *weights;
	int *domain_indices;
	
	struct tfunc_state_s tfs;

	float x0, xn;
	float y0, yn;
	int numx, numy;
	float dx, dy, dA;
	float ri, ri_ratio, ro, b;
	
	int cindexi_start, cindexi_end;
	int cindexj_start, cindexj_end;
	int cindexk_start, cindexk_end;
	
	int aindexi_start, aindexi_end;
	int aindexj_start, aindexj_end;
	int aindexk_start, aindexk_end;

	int cspani, cspano;
	int cdiami, cdiamo;
	
	int cstart, cend;
	int astart, aend;
	int fstart, fend;
	
	float circle_area, annulus_area;
};

int grid_init(struct grid_s *, float, float, int, float, float, int, float, float);
void grid_set_indices(struct grid_s *, int, int);
void grid_compute_grid_index_modulus(struct grid_s *);
void grid_compute_midpoint_distances(struct grid_s *, float, float);
void grid_get_extreme_points_on_line(struct grid_s *, float, float, char, int, int, float *, float *, float *);
void grid_get_extreme_indices(float *, int, int *, int *);
void grid_get_extreme_points(struct grid_s *);
void grid_compute_weights(struct grid_s *);

int grid_is_cell_inside_circle(const struct grid_s *, int, float, float);
int grid_is_cell_inside_circle_aliasing_zone(const struct grid_s *, int, float, float);
int grid_is_cell_in_annulus(const struct grid_s *, int, float, float, float);
int grid_is_shape_in_inner_annulus_aliasing_zone(const struct grid_s *, int, float, float, float);
int grid_is_shape_in_outer_annulus_aliasing_zone(const struct grid_s *, int, float, float, float);

void grid_clear(struct grid_s *);
void grid_set_random(struct grid_s *);
void grid_set_value(struct grid_s *, float);
void grid_set_values(struct grid_s *, const float *, int);
void grid_set_value_triangular(struct grid_s *, float, char);
int grid_set_value_fromfile(struct grid_s *, const char *);
int grid_write_state_tofile(struct grid_s *, const char *);

void grid_compute_integral(struct grid_s *, char, float *, float *);
void grid_copy_index_data(struct grid_s *);
void grid_copy_state_data(struct grid_s *);
void grid_center_state_data_at(struct grid_s *, int, int);
void grid_update_state_function(struct grid_s *);

void grid_print(const struct grid_s *);
void grid_print_domain(const struct grid_s *, char);
void grid_print_weights(const struct grid_s *, char, const char *);
void grid_print_state(const struct grid_s *, char, const char *);
void grid_print_indices(const struct grid_s *);

#endif