#include "transitionfunction.h"
#include <math.h>

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
	return 1 / (1 + expf(-4 *(x - a) / alpha));
}

float tfunc_sig2(const struct tfunc_state_s *tfs, float x, float a, float b, float alpha)
{
	return tfunc_sig1(tfs, x, a, alpha) * (1 - tfunc_sig1(tfs, x, b, alpha));
}

float tfunc_sig3(const struct tfunc_state_s *tfs, float x, float y, float m, float alpha)
{
	float sig1m = tfunc_sig1(tfs, m, 0.5, alpha);

	return x * (1 - sig1m) + y * sig1m;
}

float tfunc_apply(const struct tfunc_state_s *tfs, float n, float m)
{
	float b0 = tfs->bint[0];
	float b1 = tfs->bint[1];
	float d0 = tfs->dint[0];
	float d1 = tfs->dint[1];

	float a = tfunc_sig3(tfs, b0, d0, m, tfs->alpham);
	float b = tfunc_sig3(tfs, b1, d1, m, tfs->alpham);

	return tfunc_sig2(tfs, n, a, b, tfs->alphan);
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