#include <iostream>
#include <fstream>
#include <random>
#include <chrono>
#include <iomanip>  // for setprecision

#include <functions.hpp>
#include <model.hpp>
#include <montecarlo.hpp>
#include <structuralparameters.hpp>
#include <multilevelparameters.hpp>
#include <estimator.hpp>

using namespace std;

int main() {
    typedef mt19937_64 generator;
    unsigned int seed = std::chrono::system_clock::now().time_since_epoch().count();
    generator gen = generator(seed);

    // We define the Black and Scholes model
    double x0=100, r=0.05, sigma=0.2, K=100, T=5;
    blackAndScholesfPtr BS = blackAndScholesfPtr(new BlackAndScholes(r,sigma, x0, T));
    
    auto call = [=](double x) { return x > K ? exp(-r*T)*(x - K) : 0; };
    
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
    structParam.computeParameters(gen, std::function<double(double const &)>(call), BS, N);

    // We compute the multilevel parameters with a tolerance epsilon
    double epsilon = pow(2.0, -1);

    MultilevelParameters multilevelParamMC = MultilevelParameters(epsilon, structParam, MC);
    MultilevelParameters multilevelParamRR = MultilevelParameters(epsilon, structParam, RR);

//    // We search for the M which give the best computational cost
//    double bestCostMC = MY_INF;
//    unsigned bestMMC = 2;
//    double bestCostRR = MY_INF;
//    unsigned bestMRR = 2;

//    for (unsigned M=2; M<10; ++M){
//        MultilevelParameters multilevelParamMC = MultilevelParameters(epsilon, structParam, MC);
//        multilevelParamMC.computeOptimalParameters();
//        double tempCostMC = multilevelParamMC.computeCost();
//        bestMMC = tempCostMC < bestCostMC ? M : bestMMC;
//        bestCostMC = min(tempCostMC, bestCostMC);

//        MultilevelParameters multilevelParamRR = MultilevelParameters(epsilon, structParam, RR);
//        multilevelParamRR.computeOptimalParameters();
//        double tempCostRR = multilevelParamRR.computeCost();
//        bestMRR = tempCostRR < bestCostRR ? M : bestMRR;
//        bestCostRR = min(tempCostRR, bestCostRR);
//    }

//    // Compute the multilevel parameters for the best M
//    MultilevelParameters multilevelParamMC = MultilevelParameters(epsilon, structParam, MC, bestMMC);
//    multilevelParamMC.computeOptimalParameters();

//    MultilevelParameters multilevelParamRR = MultilevelParameters(epsilon, structParam, RR, bestMRR);
//    multilevelParamRR.computeOptimalParameters();

    // We compute the estimators
    Estimator estimatorMC(gen, std::function<double(double const &)>(call), BS, multilevelParamMC);
    Estimator estimatorRR(gen, std::function<double(double const &)>(call), BS, multilevelParamRR);

    double estimatorValueMC = estimatorMC.compute();
    double estimatorValueRR = estimatorRR.compute();

    // We display and write the parameters in a file
    multilevelParamMC.displayParameters();
    multilevelParamRR.displayParameters();
    multilevelParamMC.writeParameters("parameters.txt");
    multilevelParamRR.writeParameters("parameters.txt");

    cout << "Estimator value MC " << estimatorValueMC << endl;
    cout << "Estimator value RR " << estimatorValueRR << endl;

    double callBS = call_black_scholes(x0, K, r, sigma, T);

    cout << "Call price " << callBS << endl;
}
