#include "computations.h"
#include "solid_real_sphericalharmonics_maxl5.h"
#include "slater_poisson_maxn7.h"
#include <iomanip>
//#include <omp.h>

//#include "fastonebigheader.h"
//#include "fmath/fmath.hpp"

const double PI = 3.1415926535897932384626433832795028841971693993751;

int factoriall(int n)
{
    return (n == 1 || n == 0) ? 1 : factoriall(n - 1) * n;
}

inline double func_exp(double x)
{
    return std::exp(x);
    //return fmath::expd(x);

    //fmath::PowGenerator pw(2.0);
}

size_t grid_size(const std::vector<Site> & sites)
{
    size_t counter = 0;
    for(size_t p = 0; p < sites.size(); p++) counter += sites[p].grid.size();
    return counter;
}

size_t qbasis_size(const std::vector<Site> & sites)
{
    size_t counter = 0;
    for(size_t p = 0; p < sites.size(); p++) counter += sites[p].qbasis.size();
    return counter;
}

size_t basis_size(const std::vector<Site> & sites)
{
    size_t counter = 0;
    for(size_t p = 0; p < sites.size(); p++) counter += sites[p].basis.size();
    return counter;
}

void compute_K_Un_S(const std::vector<Site> & sites,Eigen::MatrixXd & K,Eigen::MatrixXd & Un,Eigen::MatrixXd & S)
{

    K.setZero();
    Un.setZero();
    S.setZero();

    solid_real_spherical_harmonics_maxl5 Ylma,Ylmb;

    /*MatrixXi table(n,2);
    for(size_t p = 0; p < sites.size(); p++)
    {
        for(size_t i = p; i<sites[p].grid.size(); i++)
        {
            table(index++,0) = p;
            table(index++,1) = i;
        }
    }
    */

    size_t co = 0;
    for(size_t a = 0; a < sites.size(); a++) co+=sites[a].basis.size();


    /*stobasis.resize(K.rows());
          kernelUn_rnorm2_wp_stobasis.resize(K.rows());
          kernelK_rnorm2_wp_stobasis.resize(K.rows());
          kernelS_rnorm2_wp_stobasis.resize(K.rows());
    */
    ///#pragma omp parallel for private(Ylma,Ylmb) collapse(1)
    //#pragma omp parallel for private(Ylma,Ylmb,stobasis,kernelUn_rnorm2_wp_stobasis,kernelK_rnorm2_wp_stobasis,kernelS_rnorm2_wp_stobasis) collapse(1)
    for(size_t p = 0; p < sites.size(); p++)
    {
        for(size_t i = p; i<sites[p].grid.size(); i++)
        {


            Eigen::VectorXd stobasis(co),
                  kernelUn_rnorm2_wp_stobasis(co),
                  kernelK_rnorm2_wp_stobasis(co),
                  kernelS_rnorm2_wp_stobasis(co);
            //weight and projector
            const double wp           = sites[p].weight[i]*sites[p].projector[i];

            //radial r2 factor
            const Eigen::Vector3d & r = sites[p].grid[i];
            const double rnorm        = r.norm();
            const double rnorm2       = rnorm*rnorm;

            // overlap kernel empty
            double kernelS_rnorm2_wp = rnorm2*wp;

            // nuclear kernel sum(Zc/(r-Rc))*r^2
            double kernelUn_rnorm2_wp = 0.0;
            for(size_t c = 0; c < sites.size(); c++) if (c==p)
                {
                    kernelUn_rnorm2_wp += -sites[p].zindex*rnorm;
                }
                else
                {
                    double rcnorm = (r - sites[c].R + sites[p].R).norm();
                    if(rcnorm != 0.0) kernelUn_rnorm2_wp += -sites[c].zindex*rnorm2/rcnorm; // this can be because the weight from the projector should cancel infinity by def!
                }
            kernelUn_rnorm2_wp *= wp;
            // end

            // stobasis
            size_t ia = 0;
            for(size_t a = 0; a < sites.size(); a++)
            {
                Eigen::Vector3d ra = r - sites[a].R + sites[p].R;
                double ranorm      = ra.norm();

                int maxla = sites[a].basis_maxl;
                Ylma.eval(ra.x(),ra.y(),ra.z(),ranorm,maxla);

                for(size_t oa = 0; oa < sites[a].basis.size(); oa++)
                {
                    int na       = sites[a].basis[oa].n;
                    int la       = sites[a].basis[oa].l;
                    int ma       = sites[a].basis[oa].m;
                    double zetaa = sites[a].basis[oa].zeta;

                    double factor1 = sites[a].basis[oa].norm;
                    factor1 *= std::pow(ranorm,na-1)*func_exp(-zetaa*ranorm);
                    factor1 *= Ylma.get(la,ma);

                    stobasis(ia) = factor1;
                    ia++;
                }//oa
            }//a
            // end;

            // kernelUn_rnorm2_wp_stobasis
            size_t ib = 0;
            for(size_t b = 0; b < sites.size(); b++)
            {
                Eigen::Vector3d rb = r - sites[b].R + sites[p].R;
                double rbnorm      = rb.norm();

                double kinetic_factor1 = b == p ? rnorm : (rbnorm != 0.0?rnorm2/rbnorm:0.0);
                double kinetic_factor2 = b == p ? 1.0   : (rbnorm != 0.0?rnorm2/(rbnorm*rbnorm):0.0); // 0.0 can be because the weight from the projector should cancel infinity by def!

                int maxlb = sites[b].basis_maxl;
                Ylmb.eval(rb.x(),rb.y(),rb.z(),rbnorm,maxlb);

                for(size_t ob = 0; ob < sites[b].basis.size(); ob++)
                {
                    int nb       = sites[b].basis[ob].n;
                    int lb       = sites[b].basis[ob].l;
                    int mb       = sites[b].basis[ob].m;
                    double zetab = sites[b].basis[ob].zeta;

                    double kernelK_rnorm2_wp = wp * ( 0.5*(lb*(lb+1)-nb*(nb-1))*kinetic_factor2
                                                      + zetab*nb*kinetic_factor1
                                                      - 0.5*zetab*zetab*rnorm2 ) ;

                    double sto = sites[b].basis[ob].norm *
                                 std::pow(rbnorm,nb-1) *
                                 func_exp(-zetab*rbnorm) *
                                 Ylmb.get(lb,mb);

                    kernelUn_rnorm2_wp_stobasis(ib) = kernelUn_rnorm2_wp*sto;
                    kernelK_rnorm2_wp_stobasis(ib) = kernelK_rnorm2_wp*sto;
                    kernelS_rnorm2_wp_stobasis(ib) = kernelS_rnorm2_wp*sto;

                    ib++;
                }

            }
            // end

            for(int ib = 0; ib < stobasis.rows(); ib++)
                for(int ia = 0; ia <= ib; ia++)
                {
                    //#pragma omp atomic
                    Un(ia,ib) += stobasis(ia)*kernelUn_rnorm2_wp_stobasis(ib);
                    //#pragma omp atomic
                    S(ia,ib)  += stobasis(ia)*kernelS_rnorm2_wp_stobasis(ib);
                    //#pragma omp atomic
                    K(ia,ib)  += stobasis(ia)*kernelK_rnorm2_wp_stobasis(ib);
                }

        }//i
    }//p

    for(size_t row = 0; row<S.rows(); row++)
        for(size_t col = 0; col<row; col++) S(row,col) = S(col,row);
    for(size_t row = 0; row<Un.rows(); row++)
        for(size_t col = 0; col<row; col++) Un(row,col) = Un(col,row);
    for(size_t row = 0; row<K.rows(); row++)
        for(size_t col = 0; col<row; col++) K(row,col) = K(col,row);

}


void compute_K_Un_S_slow(const std::vector<Site> & sites,Eigen::MatrixXd & K,Eigen::MatrixXd & Un,Eigen::MatrixXd & S)
{

    K.setZero();
    Un.setZero();
    S.setZero();

    solid_real_spherical_harmonics_maxl5 Ylma,Ylmb;

    ////#pragma omp parallel for private(Ylma,Ylmb) collapse(1)
    for(size_t p = 0; p < sites.size(); p++)
    {
        for(size_t i = p; i<sites[p].grid.size(); i++)
        {
            const Eigen::Vector3d & r = sites[p].grid[i];
            const double rnorm        = r.norm();
            const double rnorm2       = rnorm*rnorm;
            const double wp           = sites[p].weight[i]*sites[p].projector[i];

            // all nuclear potential * |r|^2
            double nuclear_rnorm2 = -sites[p].zindex*rnorm;
            for(size_t c = 0; c < sites.size(); c++) if(c!=p)
                {
                    //std::cerr << "c" << c << std::endl;
                    double rcnorm = (r - sites[c].R + sites[p].R).norm();
                    if(rcnorm != 0.0) nuclear_rnorm2 += -sites[c].zindex*rnorm2/rcnorm;
                    // this can be because the weight from the projector should cancel infinity by def!
                }

            size_t ib = 0;
            for(size_t b = 0; b < sites.size(); b++)
            {
                Eigen::Vector3d rb = r - sites[b].R + sites[p].R;
                double rbnorm      = rb.norm();

                double kinetic_factor1 = b == p ? rnorm : (rbnorm != 0.0?rnorm2/rbnorm:0.0);
                double kinetic_factor2 = b == p ? 1.0   : (rbnorm != 0.0?rnorm2/(rbnorm*rbnorm):0.0); // 0.0 can be because the weight from the projector should cancel infinity by def!

                int maxlb = sites[b].basis_maxl;
                Ylmb.eval(rb.x(),rb.y(),rb.z(),rbnorm,maxlb);

                for(size_t ob = 0; ob < sites[b].basis.size(); ob++)
                {
                    int nb       = sites[b].basis[ob].n;
                    int lb       = sites[b].basis[ob].l;
                    int mb       = sites[b].basis[ob].m;
                    double zetab = sites[b].basis[ob].zeta;

                    double factor2 = wp*sites[b].basis[ob].norm;
                    factor2 *= std::pow(rbnorm,nb-1)*func_exp(-zetab*rbnorm);
                    factor2 *= Ylmb.get(lb,mb); //Ylmb[lb*(1 + lb)+mb];

                    double kinetic_rnorm2 = 0.5*(lb*(lb+1)-nb*(nb-1))*kinetic_factor2 + zetab*nb*kinetic_factor1 - 0.5*zetab*zetab*rnorm2;

                    double factor2Un = factor2*nuclear_rnorm2;
                    double factor2K  = factor2*kinetic_rnorm2;
                    double factor2S  = factor2*rnorm2;

                    size_t ia = 0;
                    for(size_t a = 0; a <= b; a++)
                    {
                        Eigen::Vector3d ra = r - sites[a].R + sites[p].R;
                        double ranorm      = ra.norm();

                        int maxla = sites[a].basis_maxl;
                        Ylma.eval(ra.x(),ra.y(),ra.z(),ranorm,maxla);

                        for(size_t oa = 0; oa < sites[a].basis.size(); oa++) if (ia<=ib)
                            {
                                int na       = sites[a].basis[oa].n;
                                int la       = sites[a].basis[oa].l;
                                int ma       = sites[a].basis[oa].m;
                                double zetaa = sites[a].basis[oa].zeta;

                                double factor1 = sites[a].basis[oa].norm;
                                factor1 *= std::pow(ranorm,na-1)*func_exp(-zetaa*ranorm);
                                factor1 *= Ylma.get(la,ma); //Ylma[la*(1 + la)+ma];

                                //#pragma omp critical
                                //{
                                #pragma omp atomic
                                Un.coeffRef(ia,ib) += factor1*factor2Un;
                                #pragma omp atomic
                                S.coeffRef(ia,ib)  += factor1*factor2S;
                                #pragma omp atomic
                                K.coeffRef(ia,ib)  += factor1*factor2K;
                                //}

                                ia++;
                            }//oa
                    }//a
                    ib++;
                }//ob
            }//b
        }//i
    }//p

    for(size_t row = 0; row<S.rows(); row++)
        for(size_t col = 0; col<row; col++) S(row,col) = S(col,row);
    for(size_t row = 0; row<Un.rows(); row++)
        for(size_t col = 0; col<row; col++) Un(row,col) = Un(col,row);
    for(size_t row = 0; row<K.rows(); row++)
        for(size_t col = 0; col<row; col++) K(row,col) = K(col,row);

}

void compute_QS_scoeff(const std::vector<Site> & sites,Eigen::MatrixXd & QS,Eigen::VectorXd & scoeff)
{

    scoeff.setZero();
    {
        size_t ia = 0;
        for(size_t a = 0; a < sites.size(); a++)
        {
            for(size_t oa = 0; oa < sites[a].qbasis.size(); oa++)
            {
                int na       = sites[a].qbasis[oa].n;
                int la       = sites[a].qbasis[oa].l;
                int ma       = sites[a].qbasis[oa].m;
                double zetaa = sites[a].qbasis[oa].zeta;

                if (la == 0) scoeff(ia) = sites[a].qbasis[oa].norm*2.0*std::sqrt(PI)*factoriall(na+1)/std::pow(zetaa,na+2);

                ia++;
            }//oa
        }//a
    }

    QS.setZero();

    size_t co = 0;
    for(size_t a = 0; a < sites.size(); a++) co+=sites[a].qbasis.size();

    solid_real_spherical_harmonics_maxl5 Ylma,Ylmb;

    ///#pragma omp parallel for private(Ylma,Ylmb) collapse(1)
    //#pragma omp parallel for private(Ylma,Ylmb,stobasis,kernelUn_rnorm2_wp_stobasis,kernelK_rnorm2_wp_stobasis,kernelS_rnorm2_wp_stobasis) collapse(1)
    for(size_t p = 0; p < sites.size(); p++)
    {
        for(size_t i = p; i<sites[p].grid.size(); i++)
        {

            Eigen::VectorXd stobasis(co),
                  kernelS_rnorm2_wp_stobasis(co);
            //weight and projector
            const double wp           = sites[p].weight[i]*sites[p].projector[i];

            //radial r2 factor
            const Eigen::Vector3d & r = sites[p].grid[i];
            const double rnorm        = r.norm();
            const double rnorm2       = rnorm*rnorm;

            // overlap kernel empty
            double kernelS_rnorm2_wp = rnorm2*wp;

            // stobasis
            size_t ia = 0;
            for(size_t a = 0; a < sites.size(); a++)
            {
                Eigen::Vector3d ra = r - sites[a].R + sites[p].R;
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
                    factor1 *= std::pow(ranorm,na-1)*func_exp(-zetaa*ranorm);
                    factor1 *= Ylma.get(la,ma);

                    stobasis(ia) = factor1;
                    ia++;
                }//oa
            }//a
            // end;

            // kernelUn_rnorm2_wp_stobasis
            size_t ib = 0;
            for(size_t b = 0; b < sites.size(); b++)
            {
                Eigen::Vector3d rb = r - sites[b].R + sites[p].R;
                double rbnorm      = rb.norm();

                int maxlb = sites[b].qbasis_maxl;
                Ylmb.eval(rb.x(),rb.y(),rb.z(),rbnorm,maxlb);

                for(size_t ob = 0; ob < sites[b].qbasis.size(); ob++)
                {
                    int nb       = sites[b].qbasis[ob].n;
                    int lb       = sites[b].qbasis[ob].l;
                    int mb       = sites[b].qbasis[ob].m;
                    double zetab = sites[b].qbasis[ob].zeta;

                    double sto = sites[b].qbasis[ob].norm *
                                 std::pow(rbnorm,nb-1) *
                                 func_exp(-zetab*rbnorm) *
                                 Ylmb.get(lb,mb);

                    kernelS_rnorm2_wp_stobasis(ib) = kernelS_rnorm2_wp*sto;

                    ib++;
                }

            }
            // end

            for(int ib = 0; ib < stobasis.rows(); ib++)
                for(int ia = 0; ia <= ib; ia++)
                {
                    #pragma omp atomic
                    QS(ia,ib)  += stobasis(ia)*kernelS_rnorm2_wp_stobasis(ib);
                }

        }//i
    }//p

    for(size_t row = 0; row<QS.rows(); row++)
        for(size_t col = 0; col<row; col++) QS(row,col) = QS(col,row);

}



void compute_QS_scoeff_slow(const std::vector<Site> & sites,Eigen::MatrixXd & QS,Eigen::VectorXd & scoeff)
{

    scoeff.setZero();
    {
        size_t ia = 0;
        for(size_t a = 0; a < sites.size(); a++)
        {
            for(size_t oa = 0; oa < sites[a].qbasis.size(); oa++)
            {
                int na       = sites[a].qbasis[oa].n;
                int la       = sites[a].qbasis[oa].l;
                int ma       = sites[a].qbasis[oa].m;
                double zetaa = sites[a].qbasis[oa].zeta;

                if (la == 0) scoeff(ia) = sites[a].qbasis[oa].norm*2.0*std::sqrt(PI)*factoriall(na+1)/std::pow(zetaa,na+2);

                ia++;
            }//oa
        }//a
    }

    QS.setZero();

    solid_real_spherical_harmonics_maxl5 Ylma,Ylmb;

    #pragma omp parallel for private(Ylma,Ylmb) collapse(1)
    for(size_t p = 0; p < sites.size(); p++)
    {
        /*int NCPU,tid,NPR,NTHR;
        // get the total number of CPUs/cores available for OpenMP
        NCPU = omp_get_num_procs();
        // get the current thread ID in the parallel region
        tid = omp_get_thread_num();
        // get the total number of threads available in this parallel region
        NPR = omp_get_num_threads();
        // get the total number of threads requested
        NTHR = omp_get_max_threads();
        // only execute this on the master thread!

        std::cerr << "NCPU tid NPR NTHR" <<NCPU << ' ' << tid << ' ' << NPR << ' ' << NTHR << std::endl;
        */

        for(size_t i = 0; i < sites[p].grid.size(); i++)
        {



            const Eigen::Vector3d & r = sites[p].grid[i];
            const double rnorm        = r.norm();
            const double rnorm2       = rnorm*rnorm;
            const double wp           = sites[p].weight[i]*sites[p].projector[i];

            ////size_t ib = 0;
            for(size_t b = 0; b < sites.size(); b++)
            {
                Eigen::Vector3d rb = r - sites[b].R + sites[p].R;
                double rbnorm      = rb.norm();

                int maxlb = sites[b].qbasis_maxl;
                Ylmb.eval(rb.x(),rb.y(),rb.z(),rbnorm,maxlb);

                for(size_t ob = 0; ob < sites[b].qbasis.size(); ob++)
                {

                    int ib = sites[b].qbai[ob];

                    int nb       = sites[b].qbasis[ob].n;
                    int lb       = sites[b].qbasis[ob].l;
                    int mb       = sites[b].qbasis[ob].m;
                    double zetab = sites[b].qbasis[ob].zeta;

                    double factor2 = wp*sites[b].qbasis[ob].norm;
                    factor2 *= std::pow(rbnorm,nb-1)*func_exp(-zetab*rbnorm);
                    factor2 *= Ylmb.get(lb,mb); //[lb*(1 + lb)+mb];

                    double factor2S  = factor2*rnorm2;

                    ////size_t ia = 0;
                    for(size_t a = 0; a <= b; a++)
                    {
                        Eigen::Vector3d ra = r - sites[a].R + sites[p].R;
                        double ranorm      = ra.norm();

                        int maxla = sites[a].qbasis_maxl;
                        Ylma.eval(ra.x(),ra.y(),ra.z(),ranorm,maxla);

                        for(size_t oa = 0; oa < sites[a].qbasis.size(); oa++)
                        {

                            int ia = sites[a].qbai[oa];

                            if (ia<=ib) ;
                            else continue;


                            int na       = sites[a].qbasis[oa].n;
                            int la       = sites[a].qbasis[oa].l;
                            int ma       = sites[a].qbasis[oa].m;
                            double zetaa = sites[a].qbasis[oa].zeta;

                            //TODOD controversial
                            //if ((a==b) && (la!=lb || ma != mb)) continue;

                            double factor1 = sites[a].qbasis[oa].norm;
                            factor1 *= std::pow(ranorm,na-1)*func_exp(-zetaa*ranorm);
                            factor1 *= Ylma.get(la,ma); //[la*(1 + la)+ma];

                            factor1 *= factor2S;
                            #pragma omp atomic
                            QS.coeffRef(ia,ib)  += factor1;

                            ////ia++;
                        }//oa
                    }//a
                    ////ib++;
                }//ob
            }//b
        }//i
    }//p

    for(size_t row = 0; row<QS.rows(); row++)
        for(size_t col = 0; col<row; col++) QS(row,col) = QS(col,row);
}

void solve_SE(const Eigen::MatrixXd & H,const Eigen::MatrixXd & S,Eigen::VectorXd & eE,Eigen::MatrixXd & eC)
{
    //std::cerr << "----!-----";
    //SS = Eigen::MatrixXd::Identity(S.rows(),S.cols());
    Eigen::MatrixXd SS = 0.5*(S+S.transpose());
    Eigen::MatrixXd HH = 0.5*(H+H.transpose());
    //Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> es1(SS*SS.transpose()); // it always decompose S,allocate memory etc, performance loss
    //std::cerr << es1.operatorInverseSqrt();

    Eigen::GeneralizedSelfAdjointEigenSolver<Eigen::MatrixXd> es(HH,SS); // it always decompose S,allocate memory etc, performance loss

    //Eigen::GeneralizedSelfAdjointEigenSolver<Eigen::MatrixXd> es(H,S); // it always decompose S,allocate memory etc, performance loss
    //Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> es(H); // it always decompose S,allocate memory etc, performance loss
    //std::cerr << "----!8-----";

    //es.compute(H,S);
    std::cerr << "----!-----";
    eE = es.eigenvalues();
    eC = es.eigenvectors(); // columns are eigen vectors norm is <w|S|w>=1
}

void compute_occcupancy_with_spindegeneracy(const Eigen::VectorXd & eE,double nelectrons, double smearing, Eigen::VectorXd & occupancy, int & ihomo, int & ilumo)
{
    occupancy.setZero();

    const double degeneracy_tolerance = 0.001;

    size_t degeneracy = 1;

    double nelectron_remains = nelectrons;

    size_t ei = 0;
    size_t oi = 0;
    do
    {
        degeneracy = 1;
        while (ei<eE.rows()-1)
        {
            if(std::abs(eE[ei]-eE[++ei])<degeneracy_tolerance) degeneracy++;
            else break;
        }

        if (nelectron_remains >= 2.0*degeneracy)
        {
            for(size_t k=0; k<degeneracy; k++) occupancy[oi++] = 2.0;
            nelectron_remains -= 2.0*degeneracy;
        }
        else
        {
            for(size_t k=0; k<degeneracy; k++) occupancy[oi++] = nelectron_remains/degeneracy;
            nelectron_remains = 0.0;
        }

    }
    while (nelectron_remains>0.0 && oi < occupancy.rows());

    ihomo = oi - 1;
    ilumo = oi;

}
void compute_D(const Eigen::VectorXd & occupancy,const Eigen::MatrixXd & eC,Eigen::MatrixXd & D)
{
    D.setZero();
    for(int ib = 0; ib < D.rows(); ib++)
        for(int ia = 0; ia <= ib; ia++)
        {
            for(int k = 0; k < occupancy.rows();k++) if (occupancy[k]>0.0)
                    D(ia,ib) += occupancy[k]*eC.col(k)[ia]*eC.col(k)[ib];
        }

    for(size_t row = 0; row<D.rows(); row++)
        for(size_t col = 0; col<row; col++) D(row,col) = D(col,row);
}

void compute_rho_grid_from_D(const std::vector<Site> & sites, const Eigen::MatrixXd & D, Eigen::VectorXd & rho_grid)
{
    solid_real_spherical_harmonics_maxl5 Ylma;

    size_t co = 0;
    for(size_t a = 0; a < sites.size(); a++) co+=sites[a].basis.size();

    ///#pragma omp parallel for private(Ylma) collapse(1)
    for(size_t p = 0; p < sites.size(); p++)
    {
        for(size_t i = 0; i < sites[p].grid.size(); i++)
        {
            Eigen::VectorXd stobasis(co);

            const int pi              = sites[p].gi[i];

            const Eigen::Vector3d & r = sites[p].grid[i];
            const double rnorm        = r.norm();

            rho_grid[pi] = 0.0;

            // stobasis
            size_t ia = 0;
            for(size_t a = 0; a < sites.size(); a++)
            {
                Eigen::Vector3d ra = r - sites[a].R + sites[p].R;
                double ranorm      = ra.norm();

                int maxla = sites[a].basis_maxl;
                Ylma.eval(ra.x(),ra.y(),ra.z(),ranorm,maxla);

                for(size_t oa = 0; oa < sites[a].basis.size(); oa++)
                {
                    int na       = sites[a].basis[oa].n;
                    int la       = sites[a].basis[oa].l;
                    int ma       = sites[a].basis[oa].m;
                    double zetaa = sites[a].basis[oa].zeta;

                    double factor1 = sites[a].basis[oa].norm;
                    factor1 *= std::pow(ranorm,na-1)*func_exp(-zetaa*ranorm);
                    factor1 *= Ylma.get(la,ma);

                    stobasis(ia) = factor1;
                    ia++;
                }//oa
            }//a
            // end


            rho_grid[pi]  += stobasis.dot(D*stobasis);

        }//i
    }//p

}

void compute_Eorb(const Eigen::VectorXd & eE,const Eigen::VectorXd & occupancy, double & Eorb)
{
    Eorb = (eE.array()*occupancy.array()).sum();
}

void compute_rho_grid(const std::vector<Site> & sites, const Eigen::VectorXd & occupancy, const Eigen::MatrixXd & eC, Eigen::VectorXd & rho_grid)
{
    solid_real_spherical_harmonics_maxl5 Ylma;

    size_t co = 0;
    for(size_t a = 0; a < sites.size(); a++) co+=sites[a].basis.size();

    ////size_t pi = 0;
    #pragma omp parallel for private(Ylma) collapse(1)
    for(size_t p = 0; p < sites.size(); p++)
    {
        for(size_t i = 0; i < sites[p].grid.size(); i++)
        {
            Eigen::VectorXd stobasis(co);

            const int pi              = sites[p].gi[i];

            const Eigen::Vector3d & r = sites[p].grid[i];
            const double rnorm        = r.norm();

            rho_grid[pi] = 0.0;

            // stobasis
            size_t ia = 0;
            for(size_t a = 0; a < sites.size(); a++)
            {
                Eigen::Vector3d ra = r - sites[a].R + sites[p].R;
                double ranorm      = ra.norm();

                int maxla = sites[a].basis_maxl;
                Ylma.eval(ra.x(),ra.y(),ra.z(),ranorm,maxla);

                for(size_t oa = 0; oa < sites[a].basis.size(); oa++)
                {
                    int na       = sites[a].basis[oa].n;
                    int la       = sites[a].basis[oa].l;
                    int ma       = sites[a].basis[oa].m;
                    double zetaa = sites[a].basis[oa].zeta;

                    double factor1 = sites[a].basis[oa].norm;
                    factor1 *= std::pow(ranorm,na-1)*func_exp(-zetaa*ranorm);
                    factor1 *= Ylma.get(la,ma);

                    stobasis(ia) = factor1;
                    ia++;
                }//oa
            }//a
            // end

            for(size_t k = 0; k<occupancy.rows(); k++) if (occupancy[k] > 0.0)
                {
                    double wf = 0.0;

                    for(size_t ia = 0; ia < stobasis.rows(); ia++)  wf += eC.col(k)[ia]*stobasis(ia);
                    //wf = eC.col(k).dot(stobasis);
                    rho_grid[pi]  += occupancy[k]*wf*wf;

                }



            /*for(size_t k = 0; k<occupancy.rows(); k++) if (occupancy[k] > 0.0)
                {
                    double wf = 0.0;
                    //size_t ia = 0;
                    for(size_t a = 0; a < sites.size(); a++)
                    {
                        Eigen::Vector3d ra = r - sites[a].R + sites[p].R;
                        double ranorm      = ra.norm();

                        int maxla = sites[a].basis_maxl;
                        Ylma.eval(ra.x(),ra.y(),ra.z(),ranorm,maxla);

                        for(size_t oa = 0; oa < sites[a].basis.size(); oa++)
                        {
                            int na       = sites[a].basis[oa].n;
                            int la       = sites[a].basis[oa].l;
                            int ma       = sites[a].basis[oa].m;
                            double zetaa = sites[a].basis[oa].zeta;

                            int ia = sites[a].bai[oa];

                            double factor1 = sites[a].basis[oa].norm;
                            factor1 *= std::pow(ranorm,na-1)*func_exp(-zetaa*ranorm);
                            factor1 *= Ylma.get(la,ma); //[la*(1 + la)+ma];

                            wf += eC.col(k)[ia]*factor1;

                            ////ia++;
                        }//oa
                    }//a

                    //#pragma omp atomic
                    //#pragma omp critical
                    //{
                    rho_grid[pi]  += occupancy[k]*wf*wf;
                    //}

                }//k(filling)*/




            ////pi++;
        }//i
    }//p

}

void compute_qcoeff_from_init(const std::vector<Site> & sites,Eigen::VectorXd & qcoeff)
{
    qcoeff.setZero();

    size_t ia = 0;
    for(size_t a = 0; a < sites.size(); a++)
    {
        for(size_t oa = 0; oa < sites[a].qbasis.size(); oa++)
        {
            qcoeff[ia] = sites[a].qbasis_coeff[oa];
            ia++;
        }

    }

}

void compute_rho_grid_from_qcoeff(const std::vector<Site> & sites,const Eigen::VectorXd & qcoeff,Eigen::VectorXd & rho_grid)
{
    solid_real_spherical_harmonics_maxl5 Ylma;

    #pragma omp parallel for private(Ylma) collapse(1)
    for(size_t p = 0; p < sites.size(); p++)
    {
        for(size_t i = 0; i < sites[p].grid.size(); i++)
        {
            const int pi              = sites[p].gi[i];
            const Eigen::Vector3d & r = sites[p].grid[i];
            const double rnorm        = r.norm();
            //const double rnorm2       = rnorm*rnorm;
            //const double wp           = sites[p].weight[i]*sites[p].projector[i];

            rho_grid[pi] = 0.0;

            size_t ia = 0;

            for(size_t a = 0; a < sites.size(); a++)
            {
                Eigen::Vector3d ra = r - sites[a].R + sites[p].R;
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
                    factor1 *= std::pow(ranorm,na-1)*func_exp(-zetaa*ranorm);
                    factor1 *= Ylma.get(la,ma); //[la*(1 + la)+ma];

                    //rho_grid[pi] += sites[a].qbasis_coeff[oa]*factor1;
                    //#pragma omp atomic
                    rho_grid[pi] += qcoeff[ia]*factor1;

                    ia++;
                }//oa
            }//a
            //std::cout << std::setprecision(3) << sites[p].rho_grid[i] << ' ';
            ////pi++;
        }//i
    }//p

}

void compute_volume_integral(const std::vector<Site> & sites,const Eigen::VectorXd & grid,double & total)
{
    solid_real_spherical_harmonics_maxl5 Ylma;

    total = 0.0;
    //size_t pi = 0;
///    #pragma omp parallel for private(Ylma) collapse(1) reduction(+:total)
    for(size_t p = 0; p < sites.size(); p++)
    {
        for(size_t i = 0; i < sites[p].grid.size(); i++)
        {

            const int pi              = sites[p].gi[i];

            const Eigen::Vector3d & r = sites[p].grid[i];
            const double rnorm        = r.norm();
            const double rnorm2       = rnorm*rnorm;
            const double wp           = sites[p].weight[i]*sites[p].projector[i];

            total += wp*rnorm2*grid[pi];

            ///pi++;
        }//i
    }//p

}

void compute_tcoeff(const std::vector<Site> & sites,const Eigen::VectorXd & rho_grid,Eigen::VectorXd & tcoeff)
{
    solid_real_spherical_harmonics_maxl5 Ylma;

    tcoeff.setZero();

    size_t co = 0;
    for(size_t a = 0; a < sites.size(); a++) co+=sites[a].qbasis.size();

    #pragma omp parallel for private(Ylma) collapse(1)
    for(size_t p = 0; p < sites.size(); p++)
    {
        for(size_t i = 0; i < sites[p].grid.size(); i++)
        {
            Eigen::VectorXd stobasis(co),
                  kernel_rho_grid_rnorm2_wp_stobasis(co);

            const int pi              = sites[p].gi[i];

            const Eigen::Vector3d & r = sites[p].grid[i];
            const double rnorm        = r.norm();
            const double rnorm2       = rnorm*rnorm;
            const double wp           = sites[p].weight[i]*sites[p].projector[i];

            double kernel_rho_grid_rnorm2_wp = rnorm2*wp*rho_grid[pi];

            // stobasis
            size_t ia = 0;
            for(size_t a = 0; a < sites.size(); a++)
            {
                Eigen::Vector3d ra = r - sites[a].R + sites[p].R;
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
                    factor1 *= std::pow(ranorm,na-1)*func_exp(-zetaa*ranorm);
                    factor1 *= Ylma.get(la,ma);

                    stobasis(ia) = factor1;
                    kernel_rho_grid_rnorm2_wp_stobasis(ia)=kernel_rho_grid_rnorm2_wp*factor1;
                    ia++;
                }//oa
            }//a
            // end;

            for(size_t ia = 0; ia < tcoeff.size(); ia++)
            {
                #pragma omp atomic
                tcoeff[ia] += kernel_rho_grid_rnorm2_wp_stobasis(ia);
            }

        }//i
    }//p


}

void compute_tcoeff_slow(const std::vector<Site> & sites,const Eigen::VectorXd & rho_grid,Eigen::VectorXd & tcoeff)
{
    solid_real_spherical_harmonics_maxl5 Ylma;

    tcoeff.setZero();
    //size_t pi = 0;

    #pragma omp parallel for private(Ylma) collapse(1)
    for(size_t p = 0; p < sites.size(); p++)
    {
        for(size_t i = 0; i < sites[p].grid.size(); i++)
        {

            const int pi              = sites[p].gi[i];

            const Eigen::Vector3d & r = sites[p].grid[i];
            const double rnorm        = r.norm();
            const double rnorm2       = rnorm*rnorm;
            const double wp           = sites[p].weight[i]*sites[p].projector[i];

            size_t ia = 0;
            for(size_t a = 0; a < sites.size(); a++)
            {
                Eigen::Vector3d ra = r - sites[a].R + sites[p].R;
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
                    factor1 *= std::pow(ranorm,na-1)*func_exp(-zetaa*ranorm);
                    factor1 *= Ylma.get(la,ma); //[la*(1 + la)+ma];

                    factor1 = wp*rnorm2*factor1*rho_grid[pi];

                    #pragma omp atomic
                    tcoeff[ia] += factor1;

                    ia++;
                }//oa
            }//a
            ////pi++;
        }//i
    }//p
}

void compute_hartree_grid_from_qcoeff( const std::vector<Site> & sites,const Eigen::VectorXd & qcoeff, Eigen::VectorXd & hartree_grid)
{
    solid_real_spherical_harmonics_maxl5 Ylma;

    #pragma omp parallel for private(Ylma) collapse(1)
    for(size_t p = 0; p < sites.size(); p++)
    {
        for(size_t i = 0; i < sites[p].grid.size(); i++)
        {
            const int pi              = sites[p].gi[i];

            const Eigen::Vector3d & r = sites[p].grid[i];
            const double rnorm        = r.norm();
            const double rnorm2       = rnorm*rnorm;
            const double wp           = sites[p].weight[i]*sites[p].projector[i];

            hartree_grid[pi] = 0.0;

            size_t ia = 0;

            for(size_t a = 0; a < sites.size(); a++)
            {
                Eigen::Vector3d ra = r - sites[a].R + sites[p].R;
                double ranorm      = ra.norm();

                int maxla = sites[a].qbasis_maxl;
                Ylma.eval(ra.x(),ra.y(),ra.z(),ranorm,maxla);

                for(size_t oa = 0; oa < sites[a].qbasis.size(); oa++)
                {
                    int na       = sites[a].qbasis[oa].n;
                    int la       = sites[a].qbasis[oa].l;
                    int ma       = sites[a].qbasis[oa].m;
                    double zetaa = sites[a].qbasis[oa].zeta;

                    double factor1 = sites[a].qbasis[oa].norm/std::pow(zetaa,na+1)*4.0*PI/(2.0*la+1.0);
                    factor1 *= slater_poisson_maxn7(na,la,ranorm*zetaa);
                    factor1 *= Ylma.get(la,ma); //[la*(1 + la)+ma]; //!!!!!!!

                    //#pragma omp atomic
                    hartree_grid[pi] += qcoeff[ia] * factor1 ;


                    ia++;
                }//oa
            }//a

            ////pi++;
        }//i
    }//p
}

inline void XC_Xa(double rh, double &epsilon_xc, double &mu_xc)
{
    const double PI = 3.1415926535897932384626433832795028841971693993751;

    /****************************************************
       Xalpha exchange-correlation potential by Slater
    ****************************************************/

    int i,j;
    double dum,alpha;

    alpha = 0.7000000000000;

    dum = 3.0/PI*rh;
    epsilon_xc = -3.0/2.0*alpha*std::pow(dum,0.33333333333333333);
    mu_xc = epsilon_xc;
}

/*
ITERATION: 19
  2.00000000   2.00000000   2.00000000   0.00000000
 -3.03288837  -1.46004831  -0.14705546  -0.09326319
Etot  = -22.94421063
Ecoul = 13.66422634
Exc   = -4.58908025
dEtot = -0.00027827
dDens = 0.00664730

*/
void compute_Uh_Uxc_Ecoul_Exc_Eharris(const std::vector<Site> & sites,const Eigen::VectorXd & rho_grid, const Eigen::VectorXd & hartree_grid,Eigen::MatrixXd & Uh,Eigen::MatrixXd & Uxc,double & Ecoul,double & Exc, double & Eharris)
{

    Uh.setZero();
    Uxc.setZero();

    size_t co = 0;
    for(size_t a = 0; a < sites.size(); a++) co+=sites[a].basis.size();


    Ecoul   = 0.0;
    Exc     = 0.0;
    Eharris = 0.0;
    // -Ecoul+Exc-Ehar

    solid_real_spherical_harmonics_maxl5 Ylma,Ylmb;

    ////#pragma omp parallel for private(Ylma,Ylmb) collapse(1) reduction(+:Ecoul,Exc,Eharris)
    for(size_t p = 0; p < sites.size(); p++)
    {
        for(size_t i = 0; i < sites[p].grid.size(); i++)
        {

            Eigen::VectorXd stobasis(co),
                  kernel_hartree_rnorm2_wp_stobasis(co),
                  kernel_xc_rnorm2_wp_stobasis(co);

            const int pi              = sites[p].gi[i];

            //weight and projector
            const double wp           = sites[p].weight[i]*sites[p].projector[i];

            //radial r2 factor
            const Eigen::Vector3d & r = sites[p].grid[i];
            const double rnorm        = r.norm();
            const double rnorm2       = rnorm*rnorm;

            // hartree kernel
            double kernel_hartree_rnorm2_wp = rnorm2*wp*hartree_grid[pi];


            //if (pi != piori) throw __LINE__;

            //double hartree_rnorm2 = hartree_grid[pi]*rnorm2;

            // #pragma omp atomic
            Ecoul += 0.5*rho_grid[pi]*hartree_grid[pi]*rnorm2*wp;

            double exc;
            double vxc;

            XC_Xa(rho_grid[pi],exc,vxc);
            //exc = 0.0;
            //vxc = 0.0;

            //#pragma omp atomic
            Exc += exc*rho_grid[pi]*rnorm2*wp;

            //#pragma omp atomic
            Eharris += vxc*rho_grid[pi]*rnorm2*wp;

            double xc_rnorm2 = vxc*rnorm2;

            // xc kernel
            double kernel_xc_rnorm2_wp = rnorm2*wp*vxc;

           // stobasis
            size_t ia = 0;
            for(size_t a = 0; a < sites.size(); a++)
            {
                Eigen::Vector3d ra = r - sites[a].R + sites[p].R;
                double ranorm      = ra.norm();

                int maxla = sites[a].basis_maxl;
                Ylma.eval(ra.x(),ra.y(),ra.z(),ranorm,maxla);

                for(size_t oa = 0; oa < sites[a].basis.size(); oa++)
                {
                    int na       = sites[a].basis[oa].n;
                    int la       = sites[a].basis[oa].l;
                    int ma       = sites[a].basis[oa].m;
                    double zetaa = sites[a].basis[oa].zeta;

                    double factor1 = sites[a].basis[oa].norm;
                    factor1 *= std::pow(ranorm,na-1)*func_exp(-zetaa*ranorm);
                    factor1 *= Ylma.get(la,ma);

                    stobasis(ia) = factor1;
                    kernel_hartree_rnorm2_wp_stobasis(ia) = kernel_hartree_rnorm2_wp*factor1;
                    kernel_xc_rnorm2_wp_stobasis(ia) = kernel_xc_rnorm2_wp*factor1;

                    ia++;
                }//oa
            }//a
            // end;


            for(int ib = 0; ib < stobasis.rows(); ib++)
                for(int ia = 0; ia <= ib; ia++)
                {
                    #pragma omp atomic
                    Uh(ia,ib) += stobasis(ia)*kernel_hartree_rnorm2_wp_stobasis(ib);
                    #pragma omp atomic
                    Uxc(ia,ib)  += stobasis(ia)*kernel_xc_rnorm2_wp_stobasis(ib);
                }

        }//i
    }//p

    for(size_t row = 0; row<Uh.rows(); row++)
        for(size_t col = 0; col<row; col++) Uh(row,col) = Uh(col,row);
    for(size_t row = 0; row<Uxc.rows(); row++)
        for(size_t col = 0; col<row; col++) Uxc(row,col) = Uxc(col,row);

}

void compute_Uh_Uxc_Ecoul_Exc_Eharris_slow(const std::vector<Site> & sites,const Eigen::VectorXd & rho_grid, const Eigen::VectorXd & hartree_grid,Eigen::MatrixXd & Uh,Eigen::MatrixXd & Uxc,double & Ecoul,double & Exc, double & Eharris)
{

    Uh.setZero();
    Uxc.setZero();

    Ecoul   = 0.0;
    Exc     = 0.0;
    Eharris = 0.0;
    // -Ecoul+Exc-Ehar

    ////size_t pi = 0;
    //#pragma omp parallel for
    //shared(Uh,Uxc,Ecoul,Exc,Eharris)
    //#pragma omp for

    solid_real_spherical_harmonics_maxl5 Ylma,Ylmb;

    ////#pragma omp parallel for private(Ylma,Ylmb) collapse(1) reduction(+:Ecoul,Exc,Eharris)
    for(size_t p = 0; p < sites.size(); p++)   // ??num of projectors around a and b
    {
        /*int NCPU,tid,NPR,NTHR;
        // get the total number of CPUs/cores available for OpenMP
        NCPU = omp_get_num_procs();
        // get the current thread ID in the parallel region
        tid = omp_get_thread_num();
        // get the total number of threads available in this parallel region
        NPR = omp_get_num_threads();
        // get the total number of threads requested
        NTHR = omp_get_max_threads();
        // only execute this on the master thread!

        std::cerr << "NCPU tid NPR NTHR" <<NCPU << ' ' << tid << ' ' << NPR << ' ' << NTHR << std::endl;*/

        for(size_t i = 0; i < sites[p].grid.size(); i++)
        {

            const Eigen::Vector3d & r = sites[p].grid[i];
            const double rnorm        = r.norm();
            const double rnorm2       = rnorm*rnorm;
            const double wp           = sites[p].weight[i]*sites[p].projector[i];

            const int pi              = sites[p].gi[i];

            //if (pi != piori) throw __LINE__;

            double hartree_rnorm2 = hartree_grid[pi]*rnorm2;

            // #pragma omp atomic
            Ecoul += 0.5*rho_grid[pi]*hartree_grid[pi]*rnorm2*wp;

            double exc;
            double vxc;

            XC_Xa(rho_grid[pi],exc,vxc);
            //exc = 0.0;
            //vxc = 0.0;

            //#pragma omp atomic
            Exc += exc*rho_grid[pi]*rnorm2*wp;

            //#pragma omp atomic
            Eharris += vxc*rho_grid[pi]*rnorm2*wp;

            double xc_rnorm2 = vxc*rnorm2;

            //size_t ib = 0;
            for(size_t b = 0; b < sites.size(); b++)
            {
                Eigen::Vector3d rb = r - sites[b].R + sites[p].R;
                double rbnorm      = rb.norm();

                int maxlb = sites[b].basis_maxl;
                Ylmb.eval(rb.x(),rb.y(),rb.z(),rbnorm,maxlb);

                for(size_t ob = 0; ob < sites[b].basis.size(); ob++)
                {
                    int nb       = sites[b].basis[ob].n;
                    int lb       = sites[b].basis[ob].l;
                    int mb       = sites[b].basis[ob].m;
                    double zetab = sites[b].basis[ob].zeta;

                    int ib       = sites[b].bai[ob];

                    double factor2 = sites[b].basis[ob].norm;
                    factor2 *= std::pow(rbnorm,nb-1)*func_exp(-zetab*rbnorm);
                    factor2 *= Ylmb.get(lb,mb); //[lb*(1 + lb)+mb];

                    double factor2Uh = wp*factor2*hartree_rnorm2;
                    double factor2Uxc = wp*factor2*xc_rnorm2;

                    //size_t ia = 0;
                    for(size_t a = 0; a <= b; a++)
                    {
                        Eigen::Vector3d ra = r - sites[a].R + sites[p].R;
                        double ranorm      = ra.norm();

                        int maxla = sites[a].basis_maxl;
                        Ylma.eval(ra.x(),ra.y(),ra.z(),ranorm,maxla);

                        for(size_t oa = 0; oa < sites[a].basis.size(); oa++)
                        {
                            int na       = sites[a].basis[oa].n;
                            int la       = sites[a].basis[oa].l;
                            int ma       = sites[a].basis[oa].m;
                            double zetaa = sites[a].basis[oa].zeta;

                            int ia       = sites[a].bai[oa];

                            if (ia<=ib) ;
                            else continue;

                            double factor1 = sites[a].basis[oa].norm;
                            factor1 *= std::pow(ranorm,na-1)*func_exp(-zetaa*ranorm);
                            factor1 *= Ylma.get(la,ma); // [la*(1 + la)+ma];

                            #pragma omp atomic
                            Uh(ia,ib) += factor1*factor2Uh;

                            #pragma omp atomic
                            Uxc(ia,ib) += factor1*factor2Uxc;

                            ////ia++;
                        }//oa
                    }//a
                    ////ib++;
                }//ob
            }//b
            ////pi++;
        }//i
    }//p

    for(size_t row = 0; row<Uh.rows(); row++)
        for(size_t col = 0; col<row; col++) Uh(row,col) = Uh(col,row);
    for(size_t row = 0; row<Uxc.rows(); row++)
        for(size_t col = 0; col<row; col++) Uxc(row,col) = Uxc(col,row);

}


void solve_poisson(const std::vector<Site> & sites,const Eigen::MatrixXd & QS,const Eigen::VectorXd & scoeff,const Eigen::VectorXd & zcoeff0,const Eigen::VectorXd & rho_grid,
                   double ne,Eigen::VectorXd & tcoeff,Eigen::VectorXd & qcoeff0,Eigen::VectorXd & qcoeff,Eigen::VectorXd & hartree_grid) {

 compute_tcoeff(sites,rho_grid,tcoeff);

 qcoeff0 = QS.ldlt().solve(tcoeff);
 qcoeff = qcoeff0 - (scoeff.dot(qcoeff0) - ne)/scoeff.dot(zcoeff0) * zcoeff0;

 compute_hartree_grid_from_qcoeff(sites,qcoeff,hartree_grid);


}

