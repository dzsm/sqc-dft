#ifndef SITE_H
#define SITE_H

#include <Eigen/Eigen>
#include "site_type.h"
#include "json/json.h"

struct Site : site_type {
    Eigen::Vector3d     R;
    std::vector<double> projector;
    std::vector<int>    gi; // globail grid idx;
    std::vector<int>    bai; // globail bi idx;
    std::vector<int>    qbai; // globail qb idx;
    //mutable std::vector<double> rho_grid;
    //mutable std::vector<double> hartree_grid;

};

//struct SiteGridVariable : site_type {
 //   std::vector<double> rho_grid;
 //   std::vector<double> hartree_grid;
//};

struct MultipleSites {

    std::vector<sto_params*> basis_ref;
    std::vector<sto_params*> qbasis_ref;
    std::vector<Site> sites;
    //std::vector<SiteGridVariable> sitegridvariables;

    std::vector<std::vector<double>> rho_grid;
    std::vector<std::vector<double>> hartree_grid;


    std::vector<site_type>  sitetypes;

    double total_electron;

    void init_from_json(const Json::Value & root);



};

#endif // SITE_H
