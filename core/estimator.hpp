#ifndef ESTIMATOR_HPP
#define ESTIMATOR_HPP

#include "multilevelparameters.hpp"

using namespace std;

/**
 * @brief The Estimator class
 *
 * Class to build the Multilevel estimators.
 *
 * \f[\bar{Y}_{h,\underline{n}}^{N,q} =
 * \sum_{j=1}^R \frac{1}{N_j} \sum_{k=1}^{N_j}\langle \mathbf{T}^j, Y_{h,\underline{n}}^{(j),k}\rangle \f]
 * where
 * \f$\mathbf{T}^1 = e_1\f$ and for \f$j=2, \ldots,R\f$ \f$\mathbf{T}^j = -e_{j-1} + e_j\f$ for Monte Carlo and
 * \f$\mathbf{T}^j = -\mathbf{W}_je_{j-1} + \mathbf{W}_je_j\f$ for Richardson-Romberg.
 *
 * We remember that \f$Y_{h,n_0}^{(j),k}=0\f$, which leads to
 * \f[\frac{1}{N_1} \sum_{k=1}^{N_1} Y_h^{(1),k} +
 * \sum_{j=2}^R \frac{1}{N_j} \sum_{k=1}^{N_j} (Y_{\frac{h}{n_j}}^{(j),k} - Y_{\frac{h}{n_{j-1}}}^{(j),k})\f]
 * for Monte Carlo and
 * \f[\frac{1}{N_1} \sum_{k=1}^{N_1} Y_h^{(1),k} +
 * \sum_{j=2}^R \frac{1}{N_j} \sum_{k=1}^{N_j} \mathbf{W}_j(Y_{\frac{h}{n_j}}^{(j),k} - Y_{\frac{h}{n_{j-1}}}^{(j),k})\f]
 * for Richardson-Romberg.
 */
template <typename StateType, typename VolType, typename RandomType, typename TransitionType>
class Estimator
{
public:
    /** Constructor. */
    Estimator(mt19937_64& gen,
              std::function<double(TransitionType)> f,
              const schemePtr<StateType, VolType, RandomType, TransitionType> scheme,
              const MultilevelParameters multilevelParams);
    /** Method to compute the estimator */
    double compute();

protected:
    /** Generator for the random variable. */
    mt19937_64& m_gen;
    /** Function of the random variable that we simulate: \f$ f(X_T)\f$ */
    std::function<double(TransitionType)> m_f;
    /** Scheme for simulation.*/
    schemePtr<StateType, VolType, RandomType, TransitionType> m_scheme;
    /** Multilevel parameters which define the estimator. */
    MultilevelParameters m_multilevelParams;
};


/**
 * @param gen Generator for the random variable
 * @param f Function of the random variable that we simulate \f$f(X_T)\f$
 * @param scheme Scheme for the simulation of \f$X_t\f$
 * @param multilevelParams Multilevel parameters associated to this estimator.
 */
template <typename StateType, typename VolType, typename RandomType, typename TransitionType>
Estimator<StateType, VolType, RandomType,TransitionType>::Estimator(mt19937_64& gen,
                                                                std::function<double(TransitionType)> f,
                                                                const schemePtr<StateType, VolType, RandomType, TransitionType> scheme,
                                                                const MultilevelParameters multilevelParams):
    m_gen(gen), m_f(f), m_scheme(scheme), m_multilevelParams(multilevelParams)
{
}

/**
 * @return Computed value of the estimator \f$\bar{Y}_{h,\underline{n}}^{N,q}\f$.
 */
template <typename StateType, typename VolType, typename RandomType, typename TransitionType>
double Estimator<StateType, VolType, RandomType, TransitionType>::compute()
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
    MonteCarlo<StateType, VolType, RandomType, TransitionType> monteCarlo = MonteCarlo<StateType, VolType, RandomType, TransitionType>(m_gen, m_f, m_scheme, hInverse);
    sum += monteCarlo(N_0);

    for (unsigned int j=1; j<R; ++j){
        unsigned int N_j = ceil(N*q[j]);

        function<double(TransitionType)> newF;

        // We switch between the two cases, MC and RR
        // and in the RR case we modify the function by making the product with the weight W[j]
        switch(estimatorType){
        case RR: {
            double W_j = m_multilevelParams.getWeights()[j];
            newF = [=](TransitionType x){return W_j*m_f(x);};
            break;
        }
        case MC:{
            newF = m_f;
            break;
        }
        }
        // We make the Monte Carlo
        DoubleMonteCarlo<StateType, VolType, RandomType, TransitionType> monteCarlo = DoubleMonteCarlo<StateType, VolType, RandomType, TransitionType>(m_gen, newF, m_scheme, hInverse*n[j-1], hInverse*n[j]);
        sum += monteCarlo(N_j);
    }

    return sum;
}

#endif // ESTIMATOR_HPP
