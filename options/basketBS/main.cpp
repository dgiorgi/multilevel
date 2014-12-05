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

    typedef Eigen::VectorXd VectorType, StateType, RandomType, TransitionType ;
    typedef Eigen::MatrixXd MatrixType, VolType;

    // We define the Black and Scholes model
    VectorType x0(10);
    VectorType weights(10);
    double r = 0.05;
    VectorType sigma(10);
    double K = 100;
    double T = 1;
    MatrixType rho(10,10);

    for (int i=0; i<10; i++){
        x0[i] = 100;
        weights[i] = 1./10.;
        sigma[i] = 0.3;
        rho(i,i) = 1;
                for (int j=0; j<i; j++){
            rho(i,j) = 0.5;
            rho(j,i) = 0.5;
        }
    }

    typedef mt19937_64 generator;
    unsigned int seed = std::chrono::system_clock::now().time_since_epoch().count();
    generator gen = generator(seed);

    // We define the one dimensional Black and Scholes model
    blackAndScholesfPtr<VectorType, MatrixType> BS(new BlackAndScholes<VectorType, MatrixType>(r, sigma, rho, x0, T));
    eulerPtr<StateType, VolType, RandomType> eulerScheme(new Euler<StateType, VolType, RandomType>(BS));

    auto call = [=](StateType x) {
        double sum = weights.dot(x);
        return sum > K ? exp(-r*T)*(sum - K) : 0; };

    unsigned N = 1e6; // montecarlo simulations for the structural parameters computation

    double alpha = 0.5;
    double beta = 0.5;
    double c1 = 1;
    double H = T;

    // We define the structural parameters
    StructuralParameters structParam = StructuralParameters(alpha,beta,c1,H);
    structParam.computeParameters<StateType, VolType, RandomType, TransitionType>(gen, std::function<double(TransitionType)>(call), eulerScheme, N);

    // We compute the multilevel parameters with a tolerance epsilon
    double epsilon = pow(2.0, -1);

    MultilevelParameters multilevelParamMC = MultilevelParameters(epsilon, structParam, MC);
    MultilevelParameters multilevelParamRR = MultilevelParameters(epsilon, structParam, RR);

    // We compute the estimators
    Estimator<StateType, VolType, RandomType, TransitionType> estimatorMC(gen, std::function<double(TransitionType)>(call), eulerScheme, multilevelParamMC);
    Estimator<StateType, VolType, RandomType, TransitionType> estimatorRR(gen, std::function<double(TransitionType)>(call), eulerScheme, multilevelParamRR);

    double estimatorValueMC = estimatorMC.compute();
    double estimatorValueRR = estimatorRR.compute();

    // We display and write the parameters in a file
    multilevelParamMC.displayParameters();
    multilevelParamRR.displayParameters();
    multilevelParamMC.writeParameters("parameters.txt");
    multilevelParamRR.writeParameters("parameters.txt");

    cout << "Estimator value MC " << estimatorValueMC << endl;
    cout << "Estimator value RR " << estimatorValueRR << endl;

    double callBS = 8.880440;   // a remplacer avec un prix calculé par premia

    cout << "Call price " << callBS << endl;

    vector<double> result = vector<double>();
    result.push_back(estimatorValueMC);
    result.push_back(estimatorValueRR);
    result.push_back(callBS);

    return 0;
}
