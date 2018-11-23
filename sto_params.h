#ifndef STO_PARAMS_H
#define STO_PARAMS_H

struct sto_params {
    int n,l,m;
    double zeta;
    double norm;

    sto_params(int n,int l,int m,double zeta);
};

#endif // STO_PARAMS_H
