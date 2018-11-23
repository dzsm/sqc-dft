#include <iostream>
#include <fstream>
#include <iomanip>
#include "gausslegendre128_quadrature.h"
#include "lebedev59_quadrature.h"
#include "solid_real_sphericalharmonics_maxl5.h"
#include "slater_poisson_maxn7.h"
#include "sto_params.h"

#include "site_type.h"

#include "json/json.h"

#include <Eigen/Eigen>

#include "site.h"
#include "computations.h"
#include "timer.h"




//using namespace std;
bool cmdOptionExists(char** begin, char** end, const std::string& option)
{
    return std::find(begin, end, option) != end;
}

int main(int argc, char * argv[])
{
    try
    {

        bool test = false;        test = true;
        if(cmdOptionExists(argv, argv+argc, "-t")) test = true;

        std::string rho_cube_file = "rho.cube";
        bool cube = false;;
        if(cmdOptionExists(argv, argv+argc, "-c")) cube = true;


        //std::string input_file = "input.json";
        std::string input_file = "test/t1.json";

        std::cout << "# ################### #" << std::endl;
        std::cout << "# STO based DFT code! #" << std::endl;
        std::cout << "# ################### #" << std::endl;
        std::cout << std::endl;

        Timer general_time;
        std::cout << "# Loading input... " << std::endl;
        Json::Value root;

        std::string test_reference;
        if(test)
        {

            std::cout << "# LOADING TEST VERSION!" << std::endl;
            std::istringstream is( "{ \"labels\" : [ \"Li\" , \"Li\" ], "
                                   "\"positions\" : [ [0.0 , 0.0, 0.0 ] , [1.0, 0.0, 0.0 ] ], "
                                   "\"a1\" : [10.0,0.0,0.0], \"a2\" : [0.0,10.0,0.0], \"a3\" : [0.0,0.0,10.0], "
                                   "\"scf\" : { "
                                   "\"alpha\" : 0.5, "
                                   "\"pulay_number\" : 4, "
                                   "\"max_iteration\" : 50, "
                                   "\"total_energy_tolerance\" :0.0000000001, "
                                   "\"density_matrix_norm_tolerance\" : 0.0000000001 }, "
                                   " \"Li\" : { \"symbol\" : \"Li\" , \"zindex\" : 3.0 , \"nelectron\" : 3.0 , \"radius\" : 10.0 , "
                                               "\"mass\" : 3.0 , \"grid\" : { \"scale\" : 0.6, \"gauss\" : 32, \"lebedev\" : 29 }, "
                                               "\"name\" : \"Lithium\", \"title\" : \"Lithium (SZ)\", "
                                               "\"basis\" : [ [\"1S\" , 2.69], [\"2S\" , 0.80], [\"2P\" , 0.80] ], "
                                               "\"qbasis\" : [ [\"1S\" , 5.38 , 0.488016982152593E+00], [\"2S\" , 5.41 , -0.369857854120576E+00], [\"2S\" , 3.36 , -0.154262136775740E+00], [\"3S\" , 3.08 , -0.630098272767401E-01], [\"3S\" , 2.06 , 0.526058672124576E-01], [\"3S\" , 1.38 , 0.115711696734655E-01], [\"3S\" , 0.92 , -0.236602213283877E-03], [\"2P\" , 3.49], [\"2P\" , 1.96], [\"2P\" , 1.10], [\"3P\" , 0.92], [\"3D\" , 1.60], [\"3D\" , 1.21], [\"3D\" , 0.92], [\"4F\" , 5.00], [\"4F\" , 3.50], [\"5G\" , 3.50] ] } }");
            is >> root;
            test_reference = "# ###### TEST REFERENCE #### ITERATION:   19 ###### TEST REFERENCE ##### #\n"
                             "#       homo-2         homo-1           homo           lumo         lumo+1\n"
                             "    2.00000000     2.00000000     2.00000000     0.00000000     0.00000000\n"
                             "   -3.03288837    -1.46004831    -0.14705546    -0.09326319    -0.09326319\n"
                             "#         Etot          Ecoul            Exc          dEtot          dDens\n"
                             "  -22.94421063    13.66422634    -4.58908025    -0.00027827     0.00664729\n"
                             "# ---------------------------------------------------------------------- #\n"
                             "# Self-sonsitent iteration stopped in 6.53s = 0.11m = 0.00h\n"
                             "# ###### TEST REFERENCE #### ############### ###### TEST REFERENCE ##### #\n";
        }
        else
        {

            std::ifstream config_doc(input_file, std::ifstream::binary);
            config_doc >> root;

        }
        Computation computation;
        computation.init_from_json(root);

        std::cout << "# Number of atoms     = " << computation.multisites.sites.size() << std::endl;
        std::cout << "# Number of basis     = " << computation.multisites.basis_ref.size() << std::endl;
        std::cout << "# Number of qbasis    = " << computation.multisites.qbasis_ref.size() << std::endl;
        std::cout << "# Number of electrons = " << computation.multisites.total_electron << std::endl;

        for(int i = 0; i < computation.multisites.sites.size(); i++ )
        {
            std::cout << "# Site " << i+1 << " " << computation.multisites.sites[i].label << " " << std::endl;
        }
        std::cout << "# Input is loaded! ("<< general_time.elapsed() <<"s)"<< std::endl;

        general_time.reset();
        std::cout << "# Initializing the system... " << std::endl;
        computation.compute_Hinit2();
        std::cout << "# Number of grid points = " << computation.rho_grid.size() << std::endl;
        std::cout << "# System is initialized! ("<< general_time.elapsed() <<"s)"<< std::endl;

        std::cout << std::endl;
        std::cout << "# Starting self-consistent calculation... "<< std::endl;
        std::cout << computation.dDensdiff << std::endl;

        general_time.reset();
        for(size_t i = 0; i<computation.diis.max_iteration; i++)
        {
            computation.compute_Hscf2();
            std::cout << "# ########################## ITERATION: " << std::setw(4) << i << " ########################### #"<< std::endl;

            std::cout << '#';
            std::cout << std::setprecision(8) << std::fixed << std::setw(13)<< "homo-2" << ' ';
            std::cout << std::setprecision(8) << std::fixed << std::setw(14)<< "homo-1" << ' ';
            std::cout << std::setprecision(8) << std::fixed << std::setw(14)<< "homo" << ' ';
            std::cout << std::setprecision(8) << std::fixed << std::setw(14)<< "lumo"  << ' ';
            std::cout << std::setprecision(8) << std::fixed << std::setw(14)<< "lumo+1" << ' ';
            std::cout << std::endl;
            for(size_t k = std::max(0,int(computation.ihomo-2)); k < std::min(computation.ilumo+2,int(computation.occ.size())); k++)
            {
                std::cout << std::setprecision(8) << std::fixed << std::setw(14) << computation.occ[k] << ' ';
            }
            std::cout << std::endl;
            for(size_t k = std::max(0,int(computation.ihomo-2)); k < std::min(computation.ilumo+2,int(computation.occ.size())); k++)
            {
                std::cout << std::setprecision(8) << std::fixed << std::setw(14) << computation.eE[k] << ' ';
            }
            std::cout << std::endl;
            /*
            for(size_t k = 0; k < computation.occ.size(); k++) if (computation.occ[k==0?k:k-1]>0.0) {
                std::cout << std::setprecision(8) << std::setw(12) << std::fixed << computation.eE(k) << ' ';
            }
            std::cout << std::endl;*/

            std::cout << '#';
            std::cout << std::setprecision(8) << std::fixed << std::setw(13)<< "Etot" << ' ';
            std::cout << std::setprecision(8) << std::fixed << std::setw(14)<< "Ecoul" << ' ';
            std::cout << std::setprecision(8) << std::fixed << std::setw(14)<< "Exc" << ' ';
            std::cout << std::setprecision(8) << std::fixed << std::setw(14)<< "dEtot"  << ' ';
            std::cout << std::setprecision(8) << std::fixed << std::setw(14)<< "dDens" << ' ';
            std::cout << std::endl;

            std::cout << std::setprecision(8) << std::fixed << std::setw(14)<< computation.Etot << ' ';
            std::cout << std::setprecision(8) << std::fixed << std::setw(14)<< computation.Ecoul << ' ';
            std::cout << std::setprecision(8) << std::fixed << std::setw(14)<< computation.Exc << ' ';
            std::cout << std::setprecision(8) << std::fixed << std::setw(14)<< computation.dEtot << ' ';
            std::cout << std::setprecision(8) << std::fixed << std::setw(14)<< computation.dDensdiff << ' ';
            std::cout << std::endl;
            if(std::abs(computation.dEtot)<computation.diis.total_energy_tolerance) break;
            if(std::abs(computation.dDensdiff)<computation.diis.density_matrix_norm_tolerance) break;
            std::cout << "# ---------------------------------------------------------------------- #"<< std::endl;

        }

        std::cout << std::setprecision(2) << std::fixed << "# Self-sonsitent iteration stopped in " << general_time.elapsed() << "s = " << general_time.elapsed()/60.0<< "m = " << general_time.elapsed()/3600.0 << "h"<< std::endl;

        if (test)
        {
            std::cout << test_reference << std::endl;
        }
        if (cube)
        {
            general_time.reset();
            std::cout << "# Saving qcoeff to cube file... " << std::endl;
            save_cube_from_qcoeff(computation.multisites.sites,computation.qcoeff, {-2,-2,-2}, {3,2,2},0.1,rho_cube_file);
            std::cout << "# qcoeff is saved! ("<< general_time.elapsed() <<"s)"<< std::endl;
        }

    }
    catch (const char * e)
    {
        std::cerr << e << std::endl;
    }
    return 0;
}
