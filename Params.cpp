#include "Params.h"

/*
default values
*/
UL Params::gamma(0); 
UL Params::eta(0);
UL Params::ro(0);
UL Params::tau(0);
UL Params::ro_prim(0);
ZZ Params::tau_LS(0);


void Params::set_params(UL g, UL e, UL r, UL t, UL rp)
{
	gamma = g;
	eta = e;
	ro = r;
	tau = t;
	ro_prim = rp;
}
