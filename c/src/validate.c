#include "grid.h"
#include "transitionfunction.h"

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>


void validate_integrals(struct grid_s *grid)
{
	int len = (grid->numx - 1) * (grid->numy - 1);
	float *integrals = calloc(2 * len, sizeof(float));
	float *circle, *annulus;
	int i, j, jind;
	float cint, cmean;
	float aint, amean;
	int minidx, maxidx;
	float minrelerr, maxrelerr;
	int icenter, jcenter;
	
	icenter = grid->numx / 2;
	jcenter = grid->numy / 2;

	if (integrals != NULL) {
		circle = integrals;
		annulus = integrals + (grid->numx - 1) * (grid->numy - 1);

		printf("INFO: Starting Loop\n");
		for (j = 0; j < grid->numy - 1; j++) {
			jind = j * (grid->numx - 1);

			for (i = 0; i < grid->numx - 1; i++) {
				// printf("INFO: Centering data at (%d, %d)\n", i, j);
				grid_center_state_data_at(grid, i, j);
				
				//printf("INFO: Computing circular integral\n");
				grid_compute_integral(grid, 'i', &cint, &cmean);
				//printf("INFO: Computing annular integral\n");
				grid_compute_integral(grid, 'o', &aint, &amean);
				
				circle[jind + i] = cint;
				annulus[jind + i] = aint;
			}
		}
		
		grid_get_extreme_indices(circle, len, &minidx, &maxidx);
		printf("Circular Integrals:\n");
		printf("\tExact Area: %f\n", grid->circle_area);
		printf("\tMinimum Value: circle[%d] = %f\n", minidx, circle[minidx]);
		printf("\tMaximum Value: circle[%d] = %f\n", maxidx, circle[maxidx]);
		printf("\tAbsolute Difference: %e\n", circle[maxidx] - circle[minidx]);
		
		minrelerr = (circle[minidx] - grid->circle_area) / grid->circle_area * 100;
		maxrelerr = (circle[maxidx] - grid->circle_area) / grid->circle_area * 100;
		printf("\tMinimum Percentage Error: %+f%%\n", minrelerr);
		printf("\tMaximum Percentage Error: %+f%%\n\n", maxrelerr);
		
		grid_get_extreme_indices(annulus, len, &minidx, &maxidx);
		printf("Annular Integrals:\n");
		printf("\tExact Area: %f\n", grid->annulus_area);
		printf("\tMinimum Value: annulus[%d] = %f\n", minidx, annulus[minidx]);
		printf("\tMaximum Value: annulus[%d] = %f\n", maxidx, annulus[maxidx]);
		printf("\tAbsolute Difference: %e\n", annulus[maxidx] - annulus[minidx]);

		minrelerr = (annulus[minidx] - grid->annulus_area) / grid->annulus_area * 100;
		maxrelerr = (annulus[maxidx] - grid->annulus_area) / grid->annulus_area * 100;
		printf("\tMinimum Percentage Error: %+f%%\n", minrelerr);
		printf("\tMaximum Percentage Error: %+f%%\n\n", maxrelerr);
		
		free(integrals);
	} else {
		printf("ERROR: Insufficent memory to validate areas.\n");
	}
}

void validate(float x0, float xn, int numx, float y0, float yn, int numy, float ri_ratio, float b)
{
	struct grid_s grid;
	int init_status;
	float cint, cmean;
	float aint, amean;
	float percerr;
	int i, j;
	int icenter, jcenter;
	
	init_status = grid_init(&grid, x0, xn, numx, y0, yn, numy, ri_ratio, b);
	
	if (init_status) {
		printf("INFO: Printing grid object.\n");
		grid_print(&grid);
		
		printf("INFO: Setting indices centered at (i,j) = (0,0).\n");
		grid_set_indices(&grid, 0, 0);
		grid_print_indices(&grid);

		printf("INFO: Setting the state function value to 1.\n");
		grid_set_value(&grid, 1);
		grid_print_weights(&grid, 'i', NULL);
		grid_print_state(&grid, 'i', "%.2e ");
		
		printf("INFO: Printing inner domain.\n");
		grid_print_domain(&grid, 'i');
		
		printf("INFO: Printing outer domain.\n");
		grid_print_domain(&grid, 'o');
		
		icenter = numx / 2;
		jcenter = numy / 2;
		printf("INFO: Centering data at (i,j) = (%3d,%3d)\n", icenter, jcenter);
		grid_center_state_data_at(&grid, icenter, jcenter);
		grid_print_indices(&grid);
		grid_print_state(&grid, 'i', "%.2e ");
		
		printf("INFO: Validating Integrals\n");
		validate_integrals(&grid);
		
		printf("INFO: Freeing allocated memory.\n");
		grid_clear(&grid);
	} else {
		printf("ERROR: Unable to allocate internal memory in struct grid_s object.\n");
	}
	
}


int main(int argc, char *argv[])
{
	float x0 = -1;
	float xn = 1;
	int numx = 129;
	float y0 = -1;
	float yn = 1;
	int numy = 129;
	float ri_ratio = 0.025;
	float ri = (xn - x0) * ri_ratio;
	float ro = 3 * ri;
	float b = 0.25 * (ro - ri);
	
	validate(x0, xn, numx, y0, yn, numy, ri_ratio, b);

	return 0;
}

