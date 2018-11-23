#include "site_type.h"
#include "gausslegendre128_quadrature.h"
#include "lebedev59_quadrature.h"

void generate_radial_grid(std::vector<double> & r,std::vector<double> & w, double scale,int n)
{
    const double *xp;
    const double *wp;
    gausslegendre_points(n,&xp,&wp);

    w.resize(n);
    r.resize(n);

    scale = scale/std::log(2.0);
    for(int i = 0; i < n; i++) w[i] = wp[i]*scale*((1.0 + xp[i])/(1.0 - xp[i]) + std::log(2.0/(1.0 - xp[i])));
    for(int i = 0; i <  n; i++) r[i] = scale*(1.0+xp[i])*std::log(2.0/(1.0-xp[i]));
}

void generate_spherical_grid(std::vector<Eigen::Vector3d> & u, std::vector<double> & w,int l)
{
    const double *xp;
    const double *yp;
    const double *zp;
    const double *wp;
    int n;
    lebedev_points(l,&xp,&yp,&zp,&wp,&n);

    w.resize(n);
    u.resize(n);

    for(int i = 0; i <  n; i++) w[i] = wp[i];
    for(int i = 0; i <  n; i++) u[i] << xp[i],yp[i],zp[i];

}


void site_type::generate_grid(const std::vector<int> & lords,double scale)
{

    grid.clear();
    weight.clear();

    std::vector<double> radial_grid;
    std::vector<double> radial_weight;

    generate_radial_grid(radial_grid,radial_weight,scale,lords.size());

    for(size_t i = 0; i < lords.size(); i++)
    {
        std::vector<Eigen::Vector3d> spherical_grid;
        std::vector<double> spherical_weight;

        generate_spherical_grid(spherical_grid,spherical_weight,lords[i]);

        for(size_t k = 0; k < spherical_grid.size(); k++)
        {
            grid.push_back(radial_grid[i]*spherical_grid[k]);
            weight.push_back(radial_weight[i]*spherical_weight[k]);
        }
    }
}


/*
void compute_K_Un_S_1site(sites,Eigen::MatrixXd & K,Eigen::MatrixXd & Un,Eigen::MatrixXd & S)
{

    K.setZero();
    Un.setZero();
    S.setZero();

    solid_real_spherical_harmonics_maxl5 Ylma,Ylmb;

        #pragma omp parallel for
        for(size_t i = p; i<grid.size(); i++)
        {
            const Eigen::Vector3d & r = grid[i];
            const double rnorm        = r.norm();
            const double rnorm2       = rnorm*rnorm;
            const double wp           = weight[i];

            // all nuclear potential * |r|^2
            double nuclear_rnorm2 = -sites[p].zindex*rnorm;

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
                    factor2 *= std::pow(rbnorm,nb-1)*std::exp(-zetab*rbnorm);
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
                                factor1 *= std::pow(ranorm,na-1)*std::exp(-zetaa*ranorm);
                                factor1 *= Ylma.get(la,ma); //Ylma[la*(1 + la)+ma];

                                Un(ia,ib) += factor1*factor2Un;
                                S(ia,ib)  += factor1*factor2S;
                                K(ia,ib)  += factor1*factor2K;

                                ia++;
                            }//oa
                    }//a
                    ib++;
                }//ob
            }//b
        }//i

    for(size_t row = 0; row<S.rows(); row++)
        for(size_t col = 0; col<row; col++) S(row,col) = S(col,row);
    for(size_t row = 0; row<Un.rows(); row++)
        for(size_t col = 0; col<row; col++) Un(row,col) = Un(col,row);
    for(size_t row = 0; row<K.rows(); row++)
        for(size_t col = 0; col<row; col++) K(row,col) = K(col,row);

}

*/
