#include "computations.h"
#include "solid_real_sphericalharmonics_maxl5.h"
#include "slater_poisson_maxn7.h"
#include <iomanip>
#include <fstream>
#include <iostream>
#include "timer.h"

void Computation::compute_Hinit2()
{

    Timer timethis;

    const std::vector<Site> & sites = multisites.sites;
    int bs = basis_size(sites);
    int qbs = qbasis_size(sites);
    int gs = grid_size(sites);
    double ne = multisites.total_electron;

    std::cout << "# Basic sites grid memory: " << gs*5*8/1024 << "kB = " << gs*5*8/1024/1024 << "MB" << std::endl;

    rho_grid.resize(gs);
    hartree_grid.resize(gs);
    std::cout << "# Auxiliary grid memory:   " << gs*2*8/1024 << "kB = " << gs*2*8/1024/1024 << "MB" <<std::endl;


    QS.resize(qbs,qbs);
    qcoeff.resize(qbs);
    tcoeff.resize(qbs);
    scoeff.resize(qbs);
    qcoeff0.resize(qbs);
    zcoeff0.resize(qbs);


    K.resize(bs,bs);
    H.resize(bs,bs);
    H0.resize(bs,bs);
    D.resize(bs,bs);
    S.resize(bs,bs);
    Uh.resize(bs,bs);
    Uxc.resize(bs,bs);
    Un.resize(bs,bs);
    eC.resize(bs,bs);
    eE.resize(bs);
    occ.resize(bs);

    std::cout << "# Matrix memory:           " << (qbs*qbs+5*qbs + 9*bs*bs+2*bs)*8/1024 << "kB = " << (qbs*qbs+5*qbs + 9*bs*bs+2*bs)*8/1024/1024  << "MB" << std::endl;




    // init part>>
    compute_K_Un_S(sites,K,Un,S);
    //std::cerr << S;
    H0 = K+Un;
    Timer timethisQS;
    compute_QS_scoeff(sites,QS,scoeff);
    //std::cerr << std::setprecision(5) << std::fixed << "compute_QS_scoeff time : " << timethisQS.elapsed() << "s" << std::endl;
    timethisQS.reset();
    zcoeff0 = QS.ldlt().solve(scoeff);
    std::cerr << std::setprecision(5) << std::fixed << "QS.ldlt().solve(scoeff) time : " << timethisQS.elapsed() << "s" << std::endl;

    // <<end inti part
        std::cerr << H0 << std::endl << S;

    solve_SE(H0,S,eE,eC);

        std::cerr << '!';

    compute_occcupancy_with_spindegeneracy(eE,ne,0.0,occ,ihomo,ilumo);
    compute_Eorb(eE,occ,Eorb);

    compute_D(occ,eC,D);

    std::cerr << D;
    /*compute_rho_grid_from_D(sites,D,rho_grid);

    //compute_rho_grid(sites,occ,eC,rho_grid);
    compute_tcoeff(sites,rho_grid,tcoeff);

    qcoeff0 = QS.ldlt().solve(tcoeff);
    qcoeff = qcoeff0 - (scoeff.dot(qcoeff0) - ne)/scoeff.dot(zcoeff0) * zcoeff0;*/

    //diis.reset(); //qcoeff.rows(),1);
    /*
    // denosty and chareg conservation check
    double integrated_density1 = 0.0;
    compute_volume_integral(sites,rho_grid,integrated_density1);
    double integrated_density2 = 0.0;
    double integrated_density3 = 0.0;
    auto rho2_grid =rho_grid;
    compute_rho_grid_from_qcoeff(sites,qcoeff,rho2_grid);
    compute_volume_integral(sites,rho2_grid,integrated_density2);
    compute_volume_integral(sites,(rho2_grid-rho_grid).array().square(),integrated_density3);
    //std::cerr << ne << ' '<< scoeff.dot(qcoeff) << ' ' << integrated_density1 << ' ' <<  integrated_density2  << ' ' << std::sqrt(integrated_density3) << std::endl;
    // end check
    */

    // we got the first rho
    Etot = Eorb;
    dDensdiff = diis.error(D);
    std::cerr << dDensdiff;

    scf_iteration = 0;

    std::cerr << std::setprecision(5) << std::fixed << "Total SCF Initialization time : " << timethis.elapsed() << "s" << std::endl;

}

void Computation::compute_Hscf2()
{
    const std::vector<Site> & sites = multisites.sites;
    double ne = multisites.total_electron;

    // now sfc starts
    diis.mixing_simple(D);
    std::cerr << dDensdiff;

    Timer timethis;
    compute_rho_grid_from_D(sites,D,rho_grid);
    solve_poisson(sites,QS,scoeff,zcoeff0,rho_grid,ne,tcoeff,qcoeff0,qcoeff,hartree_grid);
    timethis.reset();

    //compute_hartree_grid_from_qcoeff(sites,qcoeff,hartree_grid);
    //std::cerr << std::setprecision(5) << std::fixed << "compute_hartree_grid_from_qcoeff : " << timethis.elapsed() << "s" << std::endl;

    //compute_rho_grid_from_qcoeff(sites,qcoeff,rho_grid);
    timethis.reset();
    compute_Uh_Uxc_Ecoul_Exc_Eharris(sites,rho_grid,hartree_grid,Uh,Uxc,Ecoul,Exc,Eharris);
    //std::cerr << std::setprecision(5) << std::fixed << "compute_Uh_Uxc_Ecoul_Exc_Eharris : " << timethis.elapsed() << "s" << std::endl;

    H = H0+Uh+Uxc;

    timethis.reset();
    solve_SE(H,S,eE,eC);
    compute_occcupancy_with_spindegeneracy(eE,ne,0.0,occ,ihomo,ilumo);
    compute_Eorb(eE,occ,Eorb);
    //std::cerr << std::setprecision(5) << std::fixed << "solve_SE : " << timethis.elapsed() << "s" << std::endl;

    Etothist = Etot;
    Etot = Eorb - Ecoul - Eharris + Exc;
    dEtot = Etot - Etothist;

    timethis.reset();

    compute_D(occ,eC,D);
    dDensdiff = diis.error(D);

    //compute_rho_grid(sites,occ,eC,rho_grid);
    //std::cerr << std::setprecision(5) << std::fixed << "compute_rho_grid : " << timethis.elapsed() << "s" << std::endl;

    /*timethis.reset();
    compute_tcoeff(sites,rho_grid,tcoeff);
    //    std::cerr << std::setprecision(5) << std::fixed << "compute_tcoeff : " << timethis.elapsed() << "s" << std::endl;

    timethis.reset();
    qcoeff0 = QS.ldlt().solve(tcoeff);
    auto qcoeffold = qcoeff;
    qcoeff = qcoeff0 - (scoeff.dot(qcoeff0) - ne)/scoeff.dot(zcoeff0) * zcoeff0;
    dDensdiff = (qcoeffold-qcoeff).norm();
    //std::cerr << std::setprecision(5) << std::fixed << "QS.ldlt().solve(tcoeff) : " << timethis.elapsed() << "s" << std::endl;


    /*
    timethis.reset();
    // denosty and chareg conservation check
    double integrated_density1 = 0.0;
    compute_volume_integral(sites,rho_grid,integrated_density1);
    std::cerr << std::setprecision(5) << std::fixed << "compute_volume_integral : " << timethis.elapsed() << "s" << std::endl;

    double integrated_density2 = 0.0;
    double integrated_density3 = 0.0;
    auto rho2_grid =rho_grid;


    timethis.reset();
    compute_rho_grid_from_qcoeff(sites,qcoeff,rho2_grid);
    std::cerr << std::setprecision(5) << std::fixed << "compute_rho_grid_from_qcoeff : " << timethis.elapsed() << "s" << std::endl;

    timethis.reset();
    compute_volume_integral(sites,rho2_grid,integrated_density2);
    compute_volume_integral(sites,(rho2_grid-rho_grid).array().square(),integrated_density3);
    std::cerr << ne << ' '<< scoeff.dot(qcoeff) << ' ' << integrated_density1 << ' ' <<  integrated_density2  << ' ' << std::sqrt(integrated_density3) << std::endl;
    // end check
    std::cerr << std::setprecision(5) << std::fixed << "conservation check : " << timethis.elapsed() << "s" << std::endl;
    */

    scf_iteration++;
}


void save_cube_from_qcoeff(const std::vector<Site> & sites,const Eigen::VectorXd & qcoeff,const Eigen::Vector3d &lo,const Eigen::Vector3d &hi,double delta,const std::string & fn)
{

    solid_real_spherical_harmonics_maxl5 Ylma;
    double wx = (hi.x()-lo.x());
    double wy = (hi.y()-lo.y());
    double wz = (hi.z()-lo.z());

    int nx = wx/delta;
    int ny = wy/delta;
    int nz = wz/delta;

    std::ofstream os(fn.c_str());

    os << ' ' << std::endl;
    os << ' ' << std::endl;
    //os << 0 << ' ' << lo.x() << ' ' << lo.y() << ' ' << lo.z() << std::endl;
    os << 1 << ' ' << 0.00 << ' ' << 0.00 << ' ' << 0.00 << std::endl;
    os << nx << ' ' << wx << ' ' << 0.00 << ' ' << 0.00 << std::endl;
    os << ny << ' ' << 0.00 << ' ' << wy << ' ' << 0.00 << std::endl;
    os << nz << ' ' << 0.00 << ' ' << 0.00 << ' ' << wz << std::endl;
    os << 1 << ' ' << 1 << ' ' << std::endl;

    int ix, iy, iz;
    for (ix=0; ix<nx; ix++)
    {
        for (iy=0; iy<ny; iy++)
        {
            for (iz=0; iz<nz; iz++)
            {

                const Eigen::Vector3d r(lo.x() + delta*ix,lo.y()  + delta*iy,lo.z() + delta*iz);
                const double rnorm        = r.norm();
                //const double rnorm2       = rnorm*rnorm;
                //const double wp           = sites[p].weight[i]*sites[p].projector[i];

                double grid_value = 0.0;

                size_t ia = 0;

                for(size_t a = 0; a < sites.size(); a++)
                {
                    Eigen::Vector3d ra = r - sites[a].R;
                    double ranorm      = ra.norm();

                    int maxla = sites[a].qbasis_maxl;
                    Ylma.eval(ra.x(),ra.y(),ra.z(),ranorm,maxla);

                    for(size_t oa = 0; oa < sites[a].qbasis.size(); oa++)
                    {
                        int na       = sites[a].qbasis[oa].n;
                        int la       = sites[a].qbasis[oa].l;
                        int ma       = sites[a].qbasis[oa].m;
                        double zetaa = sites[a].qbasis[oa].zeta;

                        double factor1 = sites[a].qbasis[oa].norm;
                        factor1 *= std::pow(ranorm,na-1)*std::exp(-zetaa*ranorm);
                        factor1 *= Ylma.get(la,ma); //[la*(1 + la)+ma];

                        grid_value += qcoeff[ia]*factor1;

                        ia++;
                    }//oa
                }//a

                os << grid_value << ' ';
                if (iz % 6 == 5) os << std::endl;
            }
            os << std::endl;
        }
    }

    os.close();
}
