#include <iostream>
#include <fstream>
#include <random>
#include <chrono>
#include <iomanip>  // for setprecision

#include "functions.hpp"
#include "model.hpp"
#include "montecarlo.hpp"
#include "structuralparameters.hpp"
#include "multilevelparameters.hpp"
#include "estimator.hpp"
#include "scheme.hpp"

#include <Eigen/Dense>

using namespace std;

int main() {

    /* ----------------------------------------------------------------------------------------------- */
    /*                  EIGEN EXAMPLE
    /* ----------------------------------------------------------------------------------------------- */

    typedef mt19937_64 generator;
    unsigned int seed = std::chrono::system_clock::now().time_since_epoch().count();
    generator gen = generator(seed);

    // We define the Black and Scholes model
    int sizeX = 1;
    int sizeW = 1;
    Eigen::VectorXd x0(sizeX); x0 << 100;
    Eigen::VectorXd weights(sizeX); weights << 1;
    double r=0.05;
    Eigen::MatrixXd sigma(sizeX,sizeW); sigma << 0.2;
    Eigen::MatrixXd rho(sizeW,sizeW); rho << 1.;
    double K = 100;
    double T = 5;

    //BlackAndScholes<Eigen::VectorXd,Eigen::MatrixXd, Eigen::MatrixXd,Eigen::VectorXd> BS(r, sigma, rho, x0, T);
    blackAndScholesfPtr<Eigen::VectorXd,Eigen::MatrixXd, Eigen::MatrixXd,Eigen::VectorXd> BS(new BlackAndScholes<Eigen::VectorXd,Eigen::MatrixXd, Eigen::MatrixXd,Eigen::VectorXd>(r, sigma, rho, x0, T));
    eulerPtr<Eigen::VectorXd,Eigen::MatrixXd, Eigen::VectorXd> eulerScheme(new Euler<Eigen::VectorXd,Eigen::MatrixXd, Eigen::VectorXd>(BS));

    auto call = [=](Eigen::VectorXd x) {
        double sum = weights.dot(x);
        return sum > K ? exp(-r*T)*(sum - K) : 0; };

    // We set the Monte Carlo object that we will use to inizialize the structural parameters
    // We recall that h = min (1,T). ... a verifier

    unsigned N = 1e6; // montecarlo simulations for the structural parameters computation

    double alpha = 1;
    double beta = 1;
    double c1 = 1;
    //    double H = min(1.,T);
    double H = T;

    // We define the structural parameters
    StructuralParameters structParam = StructuralParameters(alpha,beta,c1,H);
    structParam.computeParameters<Eigen::VectorXd,Eigen::MatrixXd, Eigen::VectorXd, Eigen::VectorXd>(gen, std::function<double(Eigen::VectorXd const &)>(call), eulerScheme, N);

    // We compute the multilevel parameters with a tolerance epsilon
    double epsilon = pow(2.0, -1);

    MultilevelParameters multilevelParamMC = MultilevelParameters(epsilon, structParam, MC);
    MultilevelParameters multilevelParamRR = MultilevelParameters(epsilon, structParam, RR);

    // We compute the estimators
    Estimator<Eigen::VectorXd,Eigen::MatrixXd, Eigen::VectorXd, Eigen::VectorXd> estimatorMC(gen, std::function<double(Eigen::VectorXd const &)>(call), eulerScheme, multilevelParamMC);
    Estimator<Eigen::VectorXd,Eigen::MatrixXd, Eigen::VectorXd, Eigen::VectorXd> estimatorRR(gen, std::function<double(Eigen::VectorXd const &)>(call), eulerScheme, multilevelParamRR);

    double estimatorValueMC = estimatorMC.compute();
    double estimatorValueRR = estimatorRR.compute();

    // We display and write the parameters in a file
    multilevelParamMC.displayParameters();
    multilevelParamRR.displayParameters();
    multilevelParamMC.writeParameters("parameters.txt");
    multilevelParamRR.writeParameters("parameters.txt");

    cout << "Estimator value MC " << estimatorValueMC << endl;
    cout << "Estimator value RR " << estimatorValueRR << endl;

    double callBS = call_black_scholes(100, 100, r, 0.2, T);

    cout << "Call price " << callBS << endl;

    return 0;


    /* ----------------------------------------------------------------------------------------------- */
    /*                  DOUBLE EXAMPLE
    /* ----------------------------------------------------------------------------------------------- */


//    // We define the Black and Scholes model
//    double x0 = 100;
//    double r=0.05;
//    double sigma = 0.2;
//    double K = 100;
//    double T = 5;

//    typedef mt19937_64 generator;
//    unsigned int seed = std::chrono::system_clock::now().time_since_epoch().count();
//    generator gen = generator(seed);

//    // We define the one dimensional Black and Scholes model
//    blackAndScholesfPtr<double, double, double, double> BS(new BlackAndScholes<double, double, double, double>(r, sigma, x0, T));
//    eulerPtr<double, double, double> eulerScheme(new Euler<double, double, double>(BS));

//    auto call = [=](double x) {return x > K ? exp(-r*T)*(x - K) : 0; };

//    unsigned N = 1e6; // montecarlo simulations for the structural parameters computation

//    double alpha = 1;
//    double beta = 1;
//    double c1 = 1;
//    //    double H = min(1.,T);
//    double H = T;

//    // We define the structural parameters
//    StructuralParameters structParam = StructuralParameters(alpha,beta,c1,H);
//    structParam.computeParameters<double, double, double>(gen, std::function<double(double const &)>(call), eulerScheme, N);

//    // We compute the multilevel parameters with a tolerance epsilon
//    double epsilon = pow(2.0, -1);

//    MultilevelParameters multilevelParamMC = MultilevelParameters(epsilon, structParam, MC);
//    MultilevelParameters multilevelParamRR = MultilevelParameters(epsilon, structParam, RR);

//    // We compute the estimators
//    Estimator<double, double, double> estimatorMC(gen, std::function<double(double const &)>(call), eulerScheme, multilevelParamMC);
//    Estimator<double, double, double> estimatorRR(gen, std::function<double(double const &)>(call), eulerScheme, multilevelParamRR);

//    double estimatorValueMC = estimatorMC.compute();
//    double estimatorValueRR = estimatorRR.compute();

//    // We display and write the parameters in a file
//    multilevelParamMC.displayParameters();
//    multilevelParamRR.displayParameters();
//    multilevelParamMC.writeParameters("parameters.txt");
//    multilevelParamRR.writeParameters("parameters.txt");

//    cout << "Estimator value MC " << estimatorValueMC << endl;
//    cout << "Estimator value RR " << estimatorValueRR << endl;

//    double callBS = call_black_scholes(100, 100, r, 0.2, T);

//    cout << "Call price " << callBS << endl;

//    vector<double> result = vector<double>();
//    result.push_back(estimatorValueMC);
//    result.push_back(estimatorValueRR);
//    result.push_back(callBS);

//    return 0;
}
