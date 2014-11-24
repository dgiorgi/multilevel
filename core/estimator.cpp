#include "estimator.hpp"

using namespace std;

/**
 * @param gen Generator for the random variable
 * @param f Function of the random variable that we simulate \f$f(X_T)\f$
 * @param scheme Scheme for the simulation of \f$X_t\f$
 * @param multilevelParams Multilevel parameters associated to this estimator.
 */
Estimator::Estimator(mt19937_64& gen,
                     std::function<double(Eigen::VectorXd)> f,
                     const schemePtr scheme,
                     const MultilevelParameters multilevelParams):
    m_gen(gen), m_f(f), m_scheme(scheme), m_multilevelParams(multilevelParams)
{
}

/**
 * @return Computed value of the estimator \f$\bar{Y}_{h,\underline{n}}^{N,q}\f$.
 */
double Estimator::compute()
{
    estimator_type estimatorType= m_multilevelParams.getEstimatorType();

    unsigned int R = m_multilevelParams.getOrder();
    double N = m_multilevelParams.getSimulationsNumber();
    vector<double> q = m_multilevelParams.getStratification();
    Refiners n = m_multilevelParams.getRefiners();
    unsigned int hInverse = m_multilevelParams.getBiaisInverse();

    double sum = 0.0;

    // We first compute the first term by single Monte Carlo simulation
    // This term is the same for both methods
    unsigned int N_0 = ceil(N*q[0]);
    MonteCarlo monteCarlo = MonteCarlo(m_gen, m_f, m_scheme, hInverse);
    sum += monteCarlo(N_0);

    for (unsigned int j=1; j<R; ++j){
        unsigned int N_j = ceil(N*q[j]);

        function<double(Eigen::VectorXd)> newF;

        // We switch between the two cases, MC and RR
        // and in the RR case we modify the function by making the product with the weight W[j]
        switch(estimatorType){
        case RR: {
            double W_j = m_multilevelParams.getWeights()[j];
            newF = [=](Eigen::VectorXd x){return W_j*m_f(x);};
            break;
        }
        case MC:{
            newF = m_f;
            break;
        }
        }
        // We make the Monte Carlo
        DoubleMonteCarlo monteCarlo = DoubleMonteCarlo(m_gen, newF, m_scheme, hInverse*n[j-1], hInverse*n[j]);
        sum += monteCarlo(N_j);
    }

    return sum;
}
