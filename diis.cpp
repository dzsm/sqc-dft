#include "computations.h"

void DIIS::init_from_json(const Json::Value & root)
{
    Dhist.clear();
    this->alpha = root["scf"].get("alpha",1.0).asDouble();
    this->pulay_number = root["scf"].get("pulay_number",4).asInt();
    this->max_iteration = max_iteration = root["scf"].get("max_iteration",50).asInt();
    this->total_energy_tolerance = root["scf"].get("total_energy_tolerance",1.0e-8).asDouble();;
    this->density_matrix_norm_tolerance = root["scf"].get("density_matrix_norm_tolerance",1.0e-8).asDouble();;

    //dmntol


}

// D get rewritten to the modified one
void DIIS:: mixing(Eigen::MatrixXd & D)
{
    if (Dhist.size()>0)
    {

        std::deque<Eigen::MatrixXd> Rhist = Dhist;
        Dhist.push_back(D);

        //Eigen::MatrixXd Dtmp = Dhist.front();
        //Dhist.pop_front();

        for(int i = 0; i < Rhist.size(); i++) Rhist[i] -= Dhist[i+1];

        Eigen::MatrixXd M(Rhist.size(),Rhist.size());
        M.setZero();
        for(int i = 0; i < Rhist.size(); i++)
            for(int j = 0; j < Rhist.size(); j++) M(i,j) = (Rhist[i]*Rhist[j]).trace();

        Eigen::VectorXd u(Rhist.size());
        u.setOnes();


        Eigen::FullPivLU<Eigen::MatrixXd> lu(M);
        if (lu.rank()<M.rows())
        {
            Dhist.pop_back();
            Dhist.pop_front();
            mixing(D);
        };
        Eigen::VectorXd  g = lu.solve(u);

        //Eigen::VectorXd  g = M.ldlt().solve(u);
        //Eigen::MatrixXd iM = M.inverse(); /// avoid error from this
        //Eigen::VectorXd  c = iM*u/(u.dot(iM*u));
        Eigen::VectorXd  c = g/(u.dot(g));

        Eigen::MatrixXd d = D;
        d.setZero();
        for(int i = 0; i < Rhist.size(); i++) d += c(i)*(Dhist[i+1]);


        Eigen::MatrixXd newD = (1.0-alpha)*D+alpha*d;

        D = newD;
        Dhist.pop_back();
        Dhist.push_back(D);

        if(Dhist.size()>pulay_number) Dhist.pop_front();
    }
    else
    {
        Dhist.push_back(D);
        return;
    }

}
// D get rewritten to the modified one
void DIIS::mixing_simple(Eigen::MatrixXd & D)
{

    double alpha = 0.8;
    if (Dhist.size()>0)
    {
        D = (1.0-alpha)*D+alpha*(Dhist.back());
        Dhist.push_back(D);
        Dhist.pop_front();
    }
    else
    {
        Dhist.push_back(D);
        return;
    }

}

double DIIS::error(const Eigen::MatrixXd & D)
{
    if (Dhist.size()>0)
    {
       // Eigen::MatrixXd A = D;
      //  Eigen::MatrixXd B = Dhist.back();

       // Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> solver1(A),solver2(B); // it always decompose S,allocate memory etc, performance loss
      //  std::cout << solver1.eigenvalues().transpose() << std::endl << solver2.eigenvalues().transpose() << std::endl;
        // return (solver1.eigenvalues()-solver2.eigenvalues()).lpNorm<Eigen::Infinity>();;
        //eC = es.eigenvectors(); // columns are eigen vectors norm is <w|S|w>=1

        return (D-Dhist.back()).array().abs().maxCoeff(); //.norm()/D.rows()/D.cols();//lpNorm<1>();//.abs().max(); //.norm();
        // return (A-B).lpNorm<Eigen::Infinity>(); //.cwiseAbs().maxCoeff(); //.norm()/D.rows()/D.cols();//lpNorm<1>();//.abs().max(); //.norm();
    }
    else
    {
        return D.array().abs().maxCoeff();

    }

}


