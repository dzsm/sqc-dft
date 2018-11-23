#include "site.h"
#include <iostream>
#include <unordered_map>

void decoding_sym(std::string sym, int & n, int & l)
{

    char c1 = sym[0];
    char c2 = sym[1];

    switch (c1)
    {
    case '1' :
        n=1;
        break;
    case '2' :
        n=2;
        break;
    case '3' :
        n=3;
        break;
    case '4' :
        n=4;
        break;
    case '5' :
        n=5;
        break;
    case '6' :
        n=6;
        break;
    case '7' :
        n=7;
        break;
    case '8' :
        n=8;
        break;
    case '9' :
        n=9;
        break;

    default :
        throw __FUNCTION__;
    }

    switch (c2)
    {
    case 'S' :
        l=0;
        break;
    case 'P' :
        l=1;
        break;
    case 'D' :
        l=2;
        break;
    case 'F' :
        l=3;
        break;
    case 'G' :
        l=4;
        break;
    case 'H' :
        l=5;
        break;

    default :
        throw __FUNCTION__;
    }


}

// r measured from site p
double projector_function_Becke(const std::vector<Eigen::Vector3d> & centersR,size_t p,const Eigen::Vector3d & r)
{
    double pp = 1.0;
    for(size_t b = 0; b < centersR.size(); b++) if (p!=b) {
            double rb = (r - centersR[b] + centersR[p]).norm();
            double mu = (r.norm()-rb)/(centersR[p] - centersR[b]).norm(); // possible diverging since a!=b)

            double f = mu;
            f =  0.5*f*(3.0-f*f);
            f =  0.5*f*(3.0-f*f);
            f =  0.5*f*(3.0-f*f);
            double s = 0.5*(1.0-f);

            pp *= s;
        }
    return pp;
}

void MultipleSites::init_from_json(const Json::Value & root) {

        std::vector<std::string> labels;
        const Json::Value & json_labels = root["labels"];
        for ( int i = 0; i < static_cast<int>(json_labels.size()); ++i )
            labels.push_back( json_labels[i].asString() );

        std::vector<Eigen::Vector3d> positions;
        const Json::Value & json_positions = root["positions"];
        for ( int i = 0; i < static_cast<int>(json_positions.size()); ++i )
            positions.push_back( { json_positions[i].get(0u,0.0).asDouble(),json_positions[i].get(1u,0.0).asDouble(),json_positions[i].get(2u,0.0).asDouble() } );


        std::vector<std::string> unique_labels;
        std::unique_copy( labels.begin(), labels.end(), std::back_inserter( unique_labels ) );
        //std::vector<site_type>    &    sitetypes;
        sitetypes.resize(unique_labels.size());
        for (int ul = 0; ul < static_cast<int>(unique_labels.size()); ++ul)
        {
            std::string label = unique_labels[ul];
            site_type & sitetype = sitetypes[ul];
            sitetype.label = label;

            int basis_maxl = 0;
            int basis_maxn = 1;

            const Json::Value & json_basis = root[label]["basis"];
            for ( int i = 0; i < static_cast<int>(json_basis.size()); ++i )
            {
                std::string nl = json_basis[i].get(0u,"1S").asString();
                double zeta = json_basis[i].get(1u,1.0).asDouble();

                int n,l;
                decoding_sym(nl,n,l);
                for(int m=-l; m<=l; m++)
                {
                    //cout << n << ' '<<l<<' '<<m<<' '<<zeta << endl;
                    sitetype.basis.push_back(sto_params(n,l,m,zeta));

                }
                basis_maxl = l>basis_maxl?l:basis_maxl;
                basis_maxn = n>basis_maxn?n:basis_maxn;
            }

            sitetype.basis_maxl = basis_maxl;
            sitetype.basis_maxn = basis_maxn;

            int qbasis_maxl = 0;
            int qbasis_maxn = 1;

            const Json::Value & json_qbasis = root[label]["qbasis"];
            for ( int i = 0; i < static_cast<int>(json_qbasis.size()); ++i )
            {
                std::string nl = json_qbasis[i].get(0u,"1S").asString();
                double zeta = json_qbasis[i].get(1u,1.0).asDouble();
                double coeff = json_qbasis[i].get(2u,0.0).asDouble();

                int n,l;
                decoding_sym(nl,n,l);
                for(int m=-l; m<=l; m++)
                {
                    //cout << n << ' '<<l<<' '<<m<<' '<<zeta << endl;
                    sitetype.qbasis.push_back(sto_params(n,l,m,zeta));
                    sitetype.qbasis_coeff.push_back(coeff);

                }
                qbasis_maxl = l>qbasis_maxl?l:qbasis_maxl;
                qbasis_maxn = n>qbasis_maxn?n:qbasis_maxn;
            }

            /*double qbssum = 0.0;
            for(int qbs = 0; qbs < sitetype.qbasis_coeff.size(); qbs++)  {
            qbssum += sitetype.qbasis_coeff[qbs]*sitetype.qbasis_coeff[qbs];
            std::cerr << sitetype.qbasis_coeff[qbs] << std::endl;
            }
            for(int qbs = 0; qbs < sitetype.qbasis_coeff.size(); qbs++) sitetype.qbasis_coeff[qbs]/=std::sqrt(qbssum);*/

            sitetype.qbasis_maxl = qbasis_maxl;
            sitetype.qbasis_maxn = qbasis_maxn;

            sitetype.nelectron = root[label].get("nelectron",0.0).asDouble();
            sitetype.mass = root[label].get("mass",0.0).asDouble();
            sitetype.radius = root[label].get("radius",0.0).asDouble();
            sitetype.symbol = root[label].get("symbol",label).asString();
            sitetype.zindex = root[label].get("zindex",0.0).asDouble();

            std::vector<int> grid_order;
            const Json::Value & json_grid = root[label]["grid"]["orders"];
            for ( int i = 0; i < static_cast<int>(json_grid.size()); ++i )
            grid_order.push_back( json_grid[i].asInt() );

            int gausspoints = root[label]["grid"].get("gauss",0).asInt();
            int lebedevpoints = root[label]["grid"].get("lebedev",3).asInt();
            for(size_t i = 0; i < gausspoints; i++ ) grid_order.push_back(lebedevpoints);

            sitetype.grid_order = grid_order;
            sitetype.grid_scale = root[label]["grid"].get("scale",0.6).asDouble();

            sitetype.generate_grid(sitetype.grid_order,sitetype.grid_scale);

        }


        // now we combine labels,positions and sitetpes to generate all the sites

        sites.resize(labels.size());
        //rho_grid.resize(labels.size());
        //hartree_grid.resize(labels.size());

        std::unordered_map<std::string,int> sitetype_map;
        for(size_t i = 0; i < sitetypes.size(); i++) {
            sitetype_map[sitetypes[i].label] = i;
        }

        total_electron = 0.0;

        int basis_counter = 0;
        int qbasis_counter = 0;
        int grid_counter = 0;
        for(size_t i = 0; i< sites.size(); i++) {
            *((site_type*)(&sites[i])) = sitetypes[sitetype_map[labels[i]]];
            sites[i].R = positions[i];
            sites[i].projector.resize(sites[i].grid.size());

            sites[i].gi.resize(sites[i].grid.size());
            for(size_t k = 0; k < sites[i].grid.size(); k++) {
                sites[i].gi[k] = grid_counter;
                grid_counter++;
            }

            //sites[i].rho_grid.resize(sites[i].grid.size());
            //sites[i].hartree_grid.resize(sites[i].grid.size());

            //rho_grid[i].resize(sites[i].grid.size());
            //hartree_grid[i].resize(sites[i].grid.size());
            sites[i].bai.resize(sites[i].basis.size());
            for(size_t k = 0; k < sites[i].basis.size(); k++) {

                sites[i].bai[k] = basis_counter;
                basis_counter++;

                basis_ref.push_back(&(sites[i].basis[k]));

            }
            sites[i].qbai.resize(sites[i].qbasis.size());
            for(size_t k = 0; k < sites[i].qbasis.size(); k++) {

                sites[i].qbai[k] = qbasis_counter;

                qbasis_counter++;
                qbasis_ref.push_back(&(sites[i].qbasis[k]));

            }


            total_electron += sites[i].nelectron;
        }

        for(size_t p = 0; p < sites.size(); p++) {
            for(size_t i = 0; i < sites[p].grid.size(); i++) {

                sites[p].projector[i] = projector_function_Becke(positions,p,sites[p].grid[i]);


                double proj_norm = 0.0;

                for(size_t c = 0; c < sites.size(); c++)
                    proj_norm += projector_function_Becke(positions,c,sites[p].grid[i]-sites[c].R);

                sites[p].projector[i] = sites[p].projector[i]/proj_norm;

                //std::cout << sites[p].projector[i] << "  ---  " << 1 << std::endl;

            }

        }

        for(size_t p = 0; p < sites.size(); p++) {

            std::cerr << sites[p].projector.size() << std::endl;
            std::cerr << sites[p].label << std::endl;
            std::cerr << sites[p].R << std::endl;
            std::cerr << sites[p].nelectron << std::endl;

            std::cerr << "------------" << std::endl;

        }


}

