#include <math.h>

#include "polymer.h"
#include "energy.h"
#include "vec3.h"

#define R_CUTOFF 1.122462048309372981433533

// A Finite Extendible Non-linear Elastic potential
static double U_fene(double r2, double k, double a){
	if( r2<a*a ){
		return -0.5*k*a*a*log(1.0-r2/a/a);
	}
	return HUGE_VAL;
}

// Truncated Lennard-Jones potential
static double U_tlj( double r2, double ep, double sig){
	double inv_r2, inv_r6;
	r2 /= sig*sig;
	if(R_CUTOFF){
		inv_r2 = 1.0/r2;
		inv_r6 = inv_r2*inv_r2*inv_r2;
		return 4.0*ep*(inv_r6*inv_r6-inv_r6)+1.0;
	}
	return 0.0;
}

double r_max_calc( const void *params){
	en_params *p = (en_params *) params;
	return p->a;
}

double r_min_calc( const void *params){
	en_params *p = (en_params *) params;
	return 0.8*p->sig;
}

double r_co_calc( const void *params){
	en_params *p = (en_params *) params;
	return R_CUTOFF*p->sig;
}

double E_bond( const double *x1, const double *x2, 
	           const void *params){
	double k, a, r2;
	en_params *p = (en_params *) params;
	k = p->k;
	a = p->a;
	r2 = dist2(x1, x2);
	return U_fene(r2, k, a);
}

double E_pair( const double *x1, const double *x2, 
			  const void *params){
	double ep, sig, r2;
	en_params *p = (en_params *) params;
	ep = p->ep;
	sig = p->sig;

	r2 = dist2(x1, x2);
	return U_tlj(r2, ep, sig);
}