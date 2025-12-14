
#ifdef USE_INTRINSICS
void set_loadmask(__m256i *loadmask, int cutoff)
{
	int allbits = 0xFFFFFFFF;
	
	switch (cutoff) {
		case 1:
			*loadmask = _mm256_setr_epi32(allbits, 0, 0, 0, 0, 0, 0, 0);
			break;
		case 2:
			*loadmask = _mm256_setr_epi32(allbits, allbits, 0, 0, 0, 0, 0, 0);
			break;
		case 3:
			*loadmask = _mm256_setr_epi32(allbits, allbits, allbits, 0, 0, 0, 0, 0);
			break;
		case 4:
			*loadmask = _mm256_setr_epi32(allbits, allbits, allbits, allbits, 0, 0, 0, 0);
			break;
		case 5:
			*loadmask = _mm256_setr_epi32(allbits, allbits, allbits, allbits, allbits, 0, 0, 0);
			break;
		case 6:
			*loadmask = _mm256_setr_epi32(allbits, allbits, allbits, allbits, allbits, allbits, 0, 0);
			break;
		case 7:
			*loadmask = _mm256_setr_epi32(allbits, allbits, allbits, allbits, allbits, allbits, allbits, 0);
	}
}

float compute_register_sum(__m256 a)
{
	float regsum = 0;
	__m256 sreg;
	__m256i idx = _mm256_setr_epi32(0, 4, 2, 6, 1, 5, 3, 7);
	
	// Compute the sum
	//     s[0] + s[1],
	//     s[2] + s[3],
	//     s[0] + s[1],
	//     s[2] + s[3],
	//     s[4] + s[5],
	//     s[6] + s[7],
	//     s[4] + s[5],
	//     s[6] + s[7]
	sreg = _mm256_hadd_ps(a, a);

	// Compute the sum
	//     s[0] + s[1] + s[2] + s[3],
	//     s[0] + s[1] + s[2] + s[3],
	//     s[0] + s[1] + s[2] + s[3],
	//     s[0] + s[1] + s[2] + s[3],
	//     s[4] + s[5] +  s[6] + s[7],
	//     s[4] + s[5] +  s[6] + s[7],
	//     s[4] + s[5] +  s[6] + s[7],
	//     s[4] + s[5] +  s[6] + s[7]
	sreg = _mm256_hadd_ps(sreg, sreg);
	
	// Swap slots 1 and 4 followed by slots 3 and 6 to get
	//     s[0] + s[1] + s[2] + s[3],
	//     s[4] + s[5] +  s[6] + s[7],
	//     s[0] + s[1] + s[2] + s[3],
	//     s[4] + s[5] +  s[6] + s[7],
	//     s[0] + s[1] + s[2] + s[3],
	//     s[4] + s[5] +  s[6] + s[7],
	//     s[0] + s[1] + s[2] + s[3],
	//     s[4] + s[5] +  s[6] + s[7]
	sreg = _mm256_permutevar8x32_ps(sreg, idx);
	
	// Compute the sum
	//     s[0] + s[1] + s[2] + s[3] + s[4] + s[5] +  s[6] + s[7],
	//     s[0] + s[1] + s[2] + s[3] + s[4] + s[5] +  s[6] + s[7],
	//     s[0] + s[1] + s[2] + s[3] + s[4] + s[5] +  s[6] + s[7],
	//     s[0] + s[1] + s[2] + s[3] + s[4] + s[5] +  s[6] + s[7],
	//     s[0] + s[1] + s[2] + s[3] + s[4] + s[5] +  s[6] + s[7],
	//     s[0] + s[1] + s[2] + s[3] + s[4] + s[5] +  s[6] + s[7],
	//     s[0] + s[1] + s[2] + s[3] + s[4] + s[5] +  s[6] + s[7],
	//     s[0] + s[1] + s[2] + s[3] + s[4] + s[5] +  s[6] + s[7]
	sreg = _mm256_hadd_ps(a, a);
	
	// Grab the sum from the first slot in the register.
	return _mm256_cvtss_f32(sreg);
}
#endif

float grid_compute_integrals(struct grid_s *grid, char which, int i, int j, float *circular_integral, float *annular_integral)
{
	int k;
	float circular_sum = 0, annular_sum = 0;
	float *data = (float *)grid->data;
	float *weights = (float *)grid->weights;
	float diffarea = grid->dx * grid->dy;
	
	grid_center_state_data_at(grid, i, j);
	
#ifdef USE_INTRINSICS

	__m256 freg;  // State function registers.
	__m256 wreg;  // Integral weight registers.
	__m256 preg;  // Product registers.
	__m256 sreg = _mm256_set1_ps(0);
	__m256i cloadmask, aloadmask
	int circular_cutoff = (grid->cend - grid->cstart) % FLOATS_PER_YMM_REG;
	int annular_cutoff = (grid->aend - grid->astart) % FLOATS_PER_YMM_REG;
	
	// Compute the circular sum.
	if (circular_cutoff > 0) {
		set_loadmask(&loadmask, circular_cutoff);
		
		freg = _mm256_maskload_ps(data + grid->cstart, loadmask);
		wreg = _mm256_maskload_ps(weights + grid->cstart, loadmask);
		preg = _mm256_mul_ps(freg, wreg);
		sreg = _mm256_add_ps(sreg, preg);
	}
	
	for (k = grid->cstart + circular_cutoff; k < grid->cend; k += FLOATS_PER_YMM_REG) {
		freg = _mm256_load_ps(data + k);
		wreg = _mm256_load_ps(weights + k);
		preg = _mm256_mul_ps(freg, wreg);
		sreg = _mm256_add_ps(sreg, preg);
	}
	
	circular_sum = compute_register_sum(sreg);
	
	// Compute the annular sum.
	sreg = _mm256_set1_ps(0);

	if (annular_cutoff > 0) {
		set_loadmask(&loadmask, annular_cutoff);
		
		freg = _mm256_maskload_ps(data + grid->astart, loadmask);
		wreg = _mm256_maskload_ps(weights + grid->astart, loadmask);
		preg = _mm256_mul_ps(freg, wreg);
		sreg = _mm256_add_ps(sreg, preg);
	}
	
	for (k = grid->astart + annular_cutoff; k < grid->aend; k += FLOATS_PER_YMM_REG) {
		freg = _mm256_load_ps(data + k);
		wreg = _mm256_load_ps(weights + k);
		preg = _mm256_mul_ps(freg, wreg);
		sreg = _mm256_add_ps(sreg, preg);
	}
	
	annular_sum = compute_register_sum(sreg);

#else
	
	for (k = grid->cstart; k < grid->cend; k++) {
		circular_sum += data[k] * weights[k];
	}
	
	for (k = grid->astart; k < grid->aend; k++) {
		annular_sum += data[k] * weights[k];
	}
	
#endif

	*circular_integral = circular_sum * diffarea / CIRCLE_AREA(grid->ri);
	*annular_integral = annular_sum * diffarea / ANNULUS_AREA(grid->ri, grid->ro);
}