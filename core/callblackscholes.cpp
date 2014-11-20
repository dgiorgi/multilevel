#include "callblackscholes.hpp"

using namespace std;

vector<double> callBlackScholes(double x0, double r, double sigma, double K, double T, unsigned int seed) {
    typedef mt19937_64 generator;
    if (seed == 0)
        seed = std::chrono::system_clock::now().time_since_epoch().count();

    generator gen = generator(seed);

    // We define the Black and Scholes model
    //double x0=100, r=0.05, sigma=0.2, K=100, T=5;
    blackAndScholesfPtr BS = blackAndScholesfPtr(new BlackAndScholes(r,sigma, x0, T));
    eulerPtr eulerScheme(new Euler(BS));
    
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
    structParam.computeParameters(gen, std::function<double(double const &)>(call), eulerScheme, N);

    // We compute the multilevel parameters with a tolerance epsilon
    double epsilon = pow(2.0, -1);

    MultilevelParameters multilevelParamMC = MultilevelParameters(epsilon, structParam, MC);
    MultilevelParameters multilevelParamRR = MultilevelParameters(epsilon, structParam, RR);

    // We compute the estimators
    Estimator estimatorMC(gen, std::function<double(double const &)>(call), eulerScheme, multilevelParamMC);
    Estimator estimatorRR(gen, std::function<double(double const &)>(call), eulerScheme, multilevelParamRR);

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

    vector<double> result = vector<double>();
    result.push_back(estimatorValueMC); 
    result.push_back(estimatorValueRR);
    result.push_back(callBS);

    return result;
}

