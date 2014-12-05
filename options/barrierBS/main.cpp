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
    double r = 0.;
    double sigma = 0.15;
    double K = 100;
    double T = 1;
    double L = 120;

    typedef mt19937_64 generator;
    unsigned int seed = std::chrono::system_clock::now().time_since_epoch().count();
    generator gen = generator(seed);

    // We define the one dimensional Black and Scholes model
    blackAndScholesfPtr<VectorType, MatrixType> BS(new BlackAndScholes<VectorType, MatrixType>(r, sigma, x0, T));
    eulerPtr<StateType, VolType, RandomType> eulerScheme(new Euler<StateType, VolType, RandomType>(BS));

    std::function<double(double, double)> myMax = [](double x,double y) {return max(x,y);};

    phiSchemePtr<StateType, VolType, RandomType> phiScheme(new PhiScheme<StateType, VolType, RandomType>(BS, eulerScheme, myMax) ) ;

    function<double(std::pair<double, double> const &)> up_out_call =
            [=](std::pair<double, double> const & x) -> double { return (x.second < L && x.first > K) ? exp(-r*T)*(x.first - K): 0; };

    unsigned N = 1e6; // montecarlo simulations for the structural parameters computation

    double alpha = 0.5;
    double beta = 0.5;
    double c1 = 1;
    double H = T;

    // We define the structural parameters
    StructuralParameters structParam = StructuralParameters(alpha,beta,c1,H);
    structParam.computeParameters<StateType, VolType, RandomType,TransitionType>(gen, std::function<double(TransitionType const &)>(up_out_call), phiScheme, N);

    string filenameMLMC = "MLMC_barrierBS.txt";
    string filenameML2R = "ML2R_barrierBS.txt";
    structParam.displayParameters();
    structParam.writeParameters(filenameMLMC);
    structParam.writeParameters(filenameML2R);

    // We compute the multilevel parameters with a tolerance epsilon
    for (int i=1; i<9; ++i){
        double epsilon = pow(2.0, -i);

        double true_value = 1.855225;

        MultilevelParameters multilevelParam2R = MultilevelParameters(epsilon, structParam, RR);
        MultilevelParameters multilevelParamMC = MultilevelParameters(epsilon, structParam, MC);

        multilevelParam2R.writeParameters(filenameML2R);
        multilevelParamMC.writeParameters(filenameMLMC);

        // We compute the estimators
        Estimator<StateType, VolType, RandomType,TransitionType> estimator2R(gen, std::function<double(TransitionType const &)>(up_out_call), phiScheme, multilevelParam2R);
        Estimator<StateType, VolType, RandomType,TransitionType> estimatorMC(gen, std::function<double(TransitionType const &)>(up_out_call), phiScheme, multilevelParamMC);

        estimator2R.L2Error(256, true_value);
        estimatorMC.L2Error(256, true_value);

        estimator2R.write(filenameML2R);
        estimatorMC.write(filenameMLMC);

        // We display and write the parameters in a file
        multilevelParamMC.displayParameters();
        estimatorMC.display();
        cout << "True value  " << true_value << endl << endl;
        multilevelParam2R.displayParameters();
        estimator2R.display();
        cout << "True value " << true_value << endl << endl;
    }

//    vector<double> result = vector<double>();
//    result.push_back(estimatorValueMC);
//    result.push_back(estimatorValueRR);
//    result.push_back(true_value);

    return 0;
}
