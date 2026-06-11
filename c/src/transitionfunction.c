#include "transitionfunction.h"
#include <math.h>
#include <stdio.h>

void tfunc_state_init(struct tfunc_state_s *tfs)
{
	tfs->alpham = 0.147;
	tfs->alphan = 0.028;
	
	tfs->bint[0] = 0.278;
	tfs->bint[1] = 0.365;
	
	tfs->dint[0] = 0.267;
	tfs->dint[1] = 0.445;
}

float tfunc_sig1(const struct tfunc_state_s *tfs, float x, float a, float alpha)
{
    float sig1_val = 1 / (1 + expf(-1 * (x - a) / (4*alpha)));
    //printf("-D- sig1: sig1_val = %f\n", sig1_val);
	return sig1_val;
}

float tfunc_sig2(const struct tfunc_state_s *tfs, float x, float a, float b, float alpha)
{
    float sig1a = tfunc_sig1(tfs, x, a, alpha);
    float sig1b = tfunc_sig1(tfs, x, b, alpha);
    float sig2_val = sig1a * (1 - sig1b);
    //printf("-D- sig2: sig1a = %f, sig1b = %f, sig2_val = %f\n", sig1a, sig1b, sig2_val);

	return sig1a * (1 - sig1b);
}

float tfunc_sig3(const struct tfunc_state_s *tfs, float x, float y, float m, float alpha)
{
	float sig1m = tfunc_sig1(tfs, m, 0.5, alpha);
	float sig3_val = x * (1 - sig1m) + y * sig1m;
    //printf("-D- sig3: sig3_val = %f, x = %f, y = %f\n", sig3_val, x, y);

	return sig3_val;
}

float tfunc_apply(const struct tfunc_state_s *tfs, float n, float m)
{
	float b0 = tfs->bint[0];
	float b1 = tfs->bint[1];
	float d0 = tfs->dint[0];
	float d1 = tfs->dint[1];

	float a = tfunc_sig3(tfs, b0, d0, m, tfs->alpham);
	float b = tfunc_sig3(tfs, b1, d1, m, tfs->alpham);
    float sig2_val = tfunc_sig2(tfs, n, a, b, tfs->alphan);

    //printf("-D- tfunc_apply: sig2_val = %f, a = %f, b = %f, m = %f, n = %f\n", sig2_val, a, b, m, n);

	return sig2_val;
}

void tfunc_set_alpham(struct tfunc_state_s *tfs, float alpham)
{
	tfs->alpham = alpham;
}

void tfunc_set_alphan(struct tfunc_state_s *tfs, float alphan)
{
	tfs->alphan = alphan;
}

void tfunc_set_birth_interval(struct tfunc_state_s *tfs, float b0, float b1)
{
	tfs->bint[0] = b0;
	tfs->bint[1] = b1;
}

void tfunc_set_death_interval(struct tfunc_state_s *tfs, float d0, float d1)
{
	tfs->dint[0] = d0;
	tfs->dint[1] = d1;
}
