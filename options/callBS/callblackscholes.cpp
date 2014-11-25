#include "callblackscholes.hpp"

using namespace std;

vector<double> callBlackScholes(double x0, double r, double sigma, double K, double T, unsigned int seed) {
    typedef mt19937_64 generator;
    if (seed == 0)
        seed = std::chrono::system_clock::now().time_since_epoch().count();

    generator gen = generator(seed);

    // We define the one dimensional Black and Scholes model
    int sizeX = 1;
    int sizeW = 1;
    Eigen::VectorXd eigen_x0(sizeX); eigen_x0 << x0;
    Eigen::VectorXd eigen_weights(sizeX); eigen_weights << 1;
    Eigen::MatrixXd eigen_sigma(sizeX,sizeW); eigen_sigma << sigma;
    Eigen::MatrixXd eigen_rho(sizeW,sizeW); eigen_rho << 1.;

    blackAndScholesfPtr<Eigen::VectorXd,Eigen::MatrixXd, Eigen::MatrixXd,Eigen::VectorXd> BS(new BlackAndScholes<Eigen::VectorXd,Eigen::MatrixXd, Eigen::MatrixXd,Eigen::VectorXd>(r, eigen_sigma, eigen_rho, eigen_x0, T));
    eulerPtr<Eigen::VectorXd,Eigen::MatrixXd, Eigen::VectorXd> eulerScheme(new Euler<Eigen::VectorXd,Eigen::MatrixXd, Eigen::VectorXd>(BS));
    
    auto call = [=](Eigen::VectorXd x) {
        double sum = eigen_weights.dot(x);
        return sum > K ? exp(-r*T)*(sum - K) : 0; };
    
    // We recall that h = min (1,T). ... a verifier

    unsigned N = 1e6; // montecarlo simulations for the structural parameters computation

    double alpha = 1;
    double beta = 1;
    double c1 = 1;
    //    double H = min(1.,T);
    double H = T;

    // We define the structural parameters
    StructuralParameters structParam = StructuralParameters(alpha,beta,c1,H);
    structParam.computeParameters<Eigen::VectorXd,Eigen::MatrixXd, Eigen::VectorXd>(gen, std::function<double(Eigen::VectorXd const &)>(call), eulerScheme, N);

    // We compute the multilevel parameters with a tolerance epsilon
    double epsilon = pow(2.0, -1);

    MultilevelParameters multilevelParamMC = MultilevelParameters(epsilon, structParam, MC);
    MultilevelParameters multilevelParamRR = MultilevelParameters(epsilon, structParam, RR);

    // We compute the estimators
    Estimator<Eigen::VectorXd,Eigen::MatrixXd, Eigen::VectorXd> estimatorMC(gen, std::function<double(Eigen::VectorXd const &)>(call), eulerScheme, multilevelParamMC);
    Estimator<Eigen::VectorXd,Eigen::MatrixXd, Eigen::VectorXd> estimatorRR(gen, std::function<double(Eigen::VectorXd const &)>(call), eulerScheme, multilevelParamRR);

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

