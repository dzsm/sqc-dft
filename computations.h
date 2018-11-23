#ifndef COMPUTATIONS_H
#define COMPUTATIONS_H

#include <Eigen/Eigen>
#include "site.h"
#include "json/json.h"
#include <iostream>


size_t grid_size(const std::vector<Site> & sites);

size_t qbasis_size(const std::vector<Site> & sites);

size_t basis_size(const std::vector<Site> & sites);

void compute_K_Un_S(const std::vector<Site> & sites,Eigen::MatrixXd & K,Eigen::MatrixXd & Un,Eigen::MatrixXd & S);

void compute_QS_scoeff(const std::vector<Site> & sites,Eigen::MatrixXd & QS,Eigen::VectorXd & scoeff);
void solve_SE(const Eigen::MatrixXd & H,const Eigen::MatrixXd & S,Eigen::VectorXd & eE,Eigen::MatrixXd & eC);
void compute_occcupancy_with_spindegeneracy(const Eigen::VectorXd & eE,double nelectrons, double smearing, Eigen::VectorXd & occupancy, int & ihomo, int & ilumo);
void compute_Eorb(const Eigen::VectorXd & eE,const Eigen::VectorXd & occupancy, double & Eorb);
void compute_rho_grid(const std::vector<Site> & sites, const Eigen::VectorXd & occupancy, const Eigen::MatrixXd & eC, Eigen::VectorXd & rho_grid);
void compute_qcoeff_from_init(const std::vector<Site> & sites,Eigen::VectorXd & qcoeff);
void compute_rho_grid_from_qcoeff(const std::vector<Site> & sites,const Eigen::VectorXd & qcoeff,Eigen::VectorXd & rho_grid);
void compute_tcoeff(const std::vector<Site> & sites,const Eigen::VectorXd & rho_grid,Eigen::VectorXd & tcoeff);
void compute_hartree_grid_from_qcoeff(const std::vector<Site> & sites,const Eigen::VectorXd & qcoeff, Eigen::VectorXd & hartree_grid);
void compute_Uh_Uxc_Ecoul_Exc_Eharris(const std::vector<Site> & sites,const Eigen::VectorXd & rho_grid, const Eigen::VectorXd & hartree_grid,Eigen::MatrixXd & Uh,Eigen::MatrixXd & Uxc,double & Ecoul,double & Exc, double & Eharris);
void compute_volume_integral(const std::vector<Site> & sites,const Eigen::VectorXd & grid,double & total);
void compute_D(const Eigen::VectorXd & occupancy,const Eigen::MatrixXd & eC,Eigen::MatrixXd & D);
void compute_rho_grid_from_D(const std::vector<Site> & sites, const Eigen::MatrixXd & D, Eigen::VectorXd & rho_grid);
void solve_poisson(const std::vector<Site> & sites,const Eigen::MatrixXd & QS,const Eigen::VectorXd & scoeff,const Eigen::VectorXd & zcoeff0,const Eigen::VectorXd & rho_grid,
                   double ne,Eigen::VectorXd & tcoeff,Eigen::VectorXd & qcoeff0,Eigen::VectorXd & qcoeff,Eigen::VectorXd & hartree_grid);


struct DIIS
{
    std::deque<Eigen::MatrixXd> Dhist;

    double alpha;
    int pulay_number;
    size_t max_iteration;
    double scf_tolerance;
    double total_energy_tolerance;
    double density_matrix_norm_tolerance;

    void init_from_json(const Json::Value & root);
    // D get rewritten to the modified one
    void mixing(Eigen::MatrixXd & D);

    // D get rewritten to the modified one
    void mixing_simple(Eigen::MatrixXd & D);

    double error(const Eigen::MatrixXd & D);

};

void save_cube_from_qcoeff(const std::vector<Site> & sites,const Eigen::VectorXd & qcoeff,const Eigen::Vector3d &lo,const Eigen::Vector3d &hi,double delta,const std::string & fn);

struct Computation
{

    MultipleSites multisites;

    Eigen::VectorXd rho_grid,hartree_grid;

    Eigen::MatrixXd QS,K,H,H0,S,Uh,Uxc,Un,eC,D;
    Eigen::VectorXd occ,eE,qcoeff,dqcoeff,tcoeff,scoeff,qcoeff0,zcoeff0;

    double Etothist;  // previous
    double dEtot;
    double Etot;      // Total energy
    double Ekin;      // Kinetic energy
    double Ecoul;     // Electron-electron interaction energy (Hartree term)
    double Eenuc;     // Electron-nucleus interaction energy
    double Exc;       // Exchange-correlation energy
    double Eharris;   // Harris energy not sure the name of it!!!!!!!!!!!!!!
    double Eorb;      // total orbital energy;

    double dDensdiff;

    int ihomo,ilumo;

    size_t scf_iteration;

    DIIS diis;

    void compute_Hinit2();
    void compute_Hscf2();

    //size_t max_iteration;
    //double scf_tolerance;

    void init_from_json(const Json::Value & root)
    {
        diis.init_from_json(root);
        multisites.init_from_json(root);
    }


};
  //-80.77604615    53.05335499   -11.69075452     0.00001861     0.26528076







#endif // COMPUTATIONS_H
