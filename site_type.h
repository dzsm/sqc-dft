#ifndef SITE_TYPE_H
#define SITE_TYPE_H

#include <Eigen/Eigen>
#include "sto_params.h"

struct site_type {
    std::string label;    // this is to what we associated from the xyz file

    std::string symbol;   // this is in the  file
    std::string filename;

    std::vector<sto_params> corebasis;
    std::vector<double> corebasis_coeff;

    std::vector<sto_params> basis;

    std::vector<sto_params> qbasis;
    std::vector<double> qbasis_coeff;

    double zindex;
    double nelectron;
    double radius;
    double mass;

    int    basis_maxl;
    int    basis_maxn;

    int    qbasis_maxl;
    int    qbasis_maxn;

    std::vector<Eigen::Vector3d> grid;
    std::vector<double>          weight;

    std::vector<int> grid_order;
    double grid_scale;
    void generate_grid(const std::vector<int> & lords,double scale);

};

#endif // SITE_TYPE_H

