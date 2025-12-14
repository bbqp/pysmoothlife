#ifndef TRANSITIONFUNCTION_H
#define TRANSITIONFUNCTION_H

struct tfunc_state_s
{
	float alpham, alphan;
	float bint[2], dint[2];
};

void tfunc_state_init(struct tfunc_state_s *);

float tfunc_sig1(const struct tfunc_state_s *, float, float, float);
float tfunc_sig2(const struct tfunc_state_s *, float, float, float, float);
float tfunc_sig3(const struct tfunc_state_s *, float, float, float, float);
float tfunc_apply(const struct tfunc_state_s *, float, float);

void tfunc_set_alpham(struct tfunc_state_s *, float);
void tfunc_set_alphan(struct tfunc_state_s *, float);
void tfunc_set_birth_interval(struct tfunc_state_s *, float, float);
void tfunc_set_death_interval(struct tfunc_state_s *, float, float);

#endif