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

    typedef double VectorType, MatrixType, StateType, VolType, RandomType, TransitionType ;

    // We define the Black and Scholes model
    double x0 = 100;
    double r=0.06;
    double sigma = 0.4;
    double K = 80;
    double T = 1;

    typedef mt19937_64 generator;
    unsigned int seed = std::chrono::system_clock::now().time_since_epoch().count();
    generator gen = generator(seed);

    // We define the one dimensional Black and Scholes model
    blackAndScholesfPtr<VectorType, MatrixType> BS(new BlackAndScholes<VectorType, MatrixType>(r, sigma, x0, T));
    eulerPtr<StateType, VolType, RandomType> eulerScheme(new Euler<StateType, VolType, RandomType>(BS));

    auto call = [=](double x) {return x > K ? exp(-r*T)*(x - K) : 0; };

    unsigned N = 1e6; // montecarlo simulations for the structural parameters computation

    double alpha = 1;
    double beta = 1;
    double c_tilde = 1;
    double H = T;

    // We define the structural parameters
    StructuralParameters structParam = StructuralParameters(alpha,beta,c_tilde,H);
    structParam.computeParameters<StateType, VolType, RandomType, TransitionType>(gen, std::function<double(TransitionType const &)>(call), eulerScheme, N);

    string filenameMLMC = "MLMC_callBS.txt";
    string filenameML2R = "ML2R_callBS.txt";
    structParam.displayParameters();
    structParam.writeParameters(filenameMLMC);
    structParam.writeParameters(filenameML2R);

    // We compute the multilevel parameters with a tolerance epsilon
    for (int i=1; i<9; ++i){
        double epsilon = pow(2.0, -i);

        double callBS = call_black_scholes(x0, K, r, sigma, T);

        MultilevelParameters multilevelParam2R = MultilevelParameters(epsilon, structParam, RR);
        MultilevelParameters multilevelParamMC = MultilevelParameters(epsilon, structParam, MC);

        multilevelParam2R.writeParameters(filenameML2R);
        multilevelParamMC.writeParameters(filenameMLMC);

        // We compute the estimators
        Estimator<StateType, VolType, RandomType, TransitionType> estimator2R(gen, std::function<double(TransitionType const &)>(call), eulerScheme, multilevelParam2R);
        Estimator<StateType, VolType, RandomType, TransitionType> estimatorMC(gen, std::function<double(TransitionType const &)>(call), eulerScheme, multilevelParamMC);

        estimator2R.L2Error(256, callBS);
        estimatorMC.L2Error(256, callBS);

        estimator2R.write(filenameML2R);
        estimatorMC.write(filenameMLMC);

        // We display and write the parameters in a file
        multilevelParam2R.displayParameters();
        estimator2R.display();
        cout << "Call price " << callBS << endl << endl;
        multilevelParamMC.displayParameters();
        estimatorMC.display();
        cout << "Call price " << callBS << endl << endl;
    }

//    vector<double> result = vector<double>();
//    result.push_back(estimatorValueMC);
//    result.push_back(estimatorValueRR);
//    result.push_back(callBS);

    return 0;
}
