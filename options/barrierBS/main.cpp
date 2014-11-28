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

    // We define the Black and Scholes model
    double x0 = 100;
    double r = 0.;
    double sigma = 0.15;
    double K = 100;
    double T = 1;
    double L = 130;

    typedef mt19937_64 generator;
    unsigned int seed = std::chrono::system_clock::now().time_since_epoch().count();
    generator gen = generator(seed);

    // We define the one dimensional Black and Scholes model
    blackAndScholesfPtr<double, double> BS(new BlackAndScholes<double, double>(r, sigma, x0, T));
    eulerPtr<double, double, double> eulerScheme(new Euler<double, double, double>(BS));

    std::function<double(double, double)> myMax = [](double x,double y) {return max(x,y);};

    phiSchemePtr<double,double,double> phiScheme(new PhiScheme<double,double,double>(BS, eulerScheme, myMax) ) ;

    function<double(std::pair<double, double> const &)> up_out_call =
            [=](std::pair<double, double> const & x) -> double { return (x.second < L && x.first > K) ? exp(-r*T)*(x.first - K): 0; };


    unsigned N = 1e6; // montecarlo simulations for the structural parameters computation

    double alpha = 0.5;
    double beta = 0.5;
    double c1 = 1;
    //    double H = min(1.,T);
    double H = T;

    // We define the structural parameters
    StructuralParameters structParam = StructuralParameters(alpha,beta,c1,H);
    structParam.computeParameters<double, double, double,std::pair<double, double>>(gen, std::function<double(std::pair<double, double> const &)>(up_out_call), phiScheme, N);

    // We compute the multilevel parameters with a tolerance epsilon
    double epsilon = pow(2.0, -1);

    MultilevelParameters multilevelParamMC = MultilevelParameters(epsilon, structParam, MC);
    MultilevelParameters multilevelParamRR = MultilevelParameters(epsilon, structParam, RR);

    // We compute the estimators
    Estimator<double, double, double, std::pair<double, double>> estimatorMC(gen, std::function<double(std::pair<double, double> const &)>(up_out_call), phiScheme, multilevelParamMC);
    Estimator<double, double, double, std::pair<double, double>> estimatorRR(gen, std::function<double(std::pair<double, double> const &)>(up_out_call), phiScheme, multilevelParamRR);

    double estimatorValueMC = estimatorMC.compute();
    double estimatorValueRR = estimatorRR.compute();

    // We display and write the parameters in a file
    multilevelParamMC.displayParameters();
    multilevelParamRR.displayParameters();
    multilevelParamMC.writeParameters("parameters.txt");
    multilevelParamRR.writeParameters("parameters.txt");

    cout << "Estimator value MC " << estimatorValueMC << endl;
    cout << "Estimator value RR " << estimatorValueRR << endl;

    double true_value = 3.869693;

    cout << "True value " << 3.869693 << endl;

    vector<double> result = vector<double>();
    result.push_back(estimatorValueMC);
    result.push_back(estimatorValueRR);
    result.push_back(true_value);

    return 0;
}
