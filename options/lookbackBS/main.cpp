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

    typedef double VectorType, MatrixType, StateType, VolType, RandomType;
    typedef std::pair<double, double> TransitionType;

    // We define the Black and Scholes model
    double x0 = 100;
    double r=0.15;
    double sigma = 0.1;
    double T = 1;
    double lambda = 1.1;

    typedef mt19937_64 generator;
    unsigned int seed = std::chrono::system_clock::now().time_since_epoch().count();
    generator gen = generator(seed);

    // We define the one dimensional Black and Scholes model
    blackAndScholesfPtr<VectorType, MatrixType> BS(new BlackAndScholes<VectorType, MatrixType>(r, sigma, x0, T));
    eulerPtr<StateType, VolType, RandomType> eulerScheme(new Euler<StateType, VolType, RandomType>(BS));

    std::function<double(double, double)> myMin = [](double x,double y) {return min(x,y);};

    phiSchemePtr<StateType, VolType, RandomType> phiScheme(new PhiScheme<StateType, VolType, RandomType>(BS, eulerScheme, myMin) ) ;

    function<double(std::pair<double, double> const &)> lookback_call =
            [=](std::pair<double, double> const & x) -> double { return x.first > lambda * x.second ? exp(-r*T)*(x.first - lambda*x.second) : 0; };

    double true_value = x0 * call_black_scholes(1, lambda, r, sigma, T)
            + lambda * sigma*sigma / (2.*r) * x0 *
             put_black_scholes(pow(lambda, 2.*r/(sigma*sigma)), 1., r, 2.*r/sigma, T);

    unsigned N = 1e6; // montecarlo simulations for the structural parameters computation

    double alpha = 0.5;
    double beta = 1;
    double c1 = 1;
    double H = T;

    // We define the structural parameters  
    StructuralParameters structParam = StructuralParameters(alpha,beta,c1,H);
    structParam.computeParameters<StateType, VolType, RandomType,TransitionType>(gen, std::function<double(TransitionType const &)>(lookback_call), phiScheme, N);

    string filenameMLMC = "MLMC_callBS.txt";
    string filenameML2R = "ML2R_callBS.txt";
    structParam.displayParameters();
    structParam.writeParameters(filenameMLMC);
    structParam.writeParameters(filenameML2R);

    // We compute the multilevel parameters with a tolerance epsilon
    for (int i=1; i<9; ++i){
        double epsilon = pow(2.0, -i);

        MultilevelParameters multilevelParam2R = MultilevelParameters(epsilon, structParam, RR);
        MultilevelParameters multilevelParamMC = MultilevelParameters(epsilon, structParam, MC);

        multilevelParam2R.writeParameters(filenameML2R);
        multilevelParamMC.writeParameters(filenameMLMC);

        // We compute the estimators
        Estimator<StateType, VolType, RandomType,TransitionType> estimator2R(gen, std::function<double(TransitionType const &)>(lookback_call), phiScheme, multilevelParam2R);
        Estimator<StateType, VolType, RandomType,TransitionType> estimatorMC(gen, std::function<double(TransitionType const &)>(lookback_call), phiScheme, multilevelParamMC);

        estimator2R.L2Error(256, true_value);
        estimatorMC.L2Error(256, true_value);

        estimator2R.write(filenameML2R);
        estimatorMC.write(filenameMLMC);

        // We display and write the parameters in a file
        multilevelParam2R.displayParameters();
        estimator2R.display();
        cout << "Look back price " << true_value << endl << endl;
        multilevelParamMC.displayParameters();
        estimatorMC.display();
        cout << "Look back price " << true_value << endl << endl;
    }

//    vector<double> result = vector<double>();
//    result.push_back(estimatorValueMC);
//    result.push_back(estimatorValueRR);
//    result.push_back(callBS);

    return 0;
}
