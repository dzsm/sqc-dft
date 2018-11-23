#include "sto_params.h"

#include <cmath>

int factorial(int n)
{
    return (n == 1 || n == 0) ? 1 : factorial(n - 1) * n;
}

double radial_norm(double zeta,int n)
{
    return std::pow(2.0*zeta,0.5+(double)(n))/std::sqrt((double)(factorial(2*n)));
}

sto_params::sto_params(int n,int l,int m,double zeta) : n(n),l(l),m(m),zeta(zeta)
{
    norm = radial_norm(zeta,n);
}


