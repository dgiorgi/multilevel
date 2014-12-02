#ifndef ESTIMATOR_HPP
#define ESTIMATOR_HPP

#include <chrono>

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
    /** Method to compute the L2 error of the estimator */
    double L2Error(int N, double trueValue);
    /** Method to display the L2 error, biais and variance */
    void display();
    /** Method to write the estimator results in a file */
    void write(const string fileName);

protected:
    /** Generator for the random variable. */
    mt19937_64& m_gen;
    /** Function of the random variable that we simulate: \f$ f(X_T)\f$ */
    std::function<double(TransitionType)> m_f;
    /** Scheme for simulation.*/
    schemePtr<StateType, VolType, RandomType, TransitionType> m_scheme;
    /** Multilevel parameters which define the estimator. */
    MultilevelParameters m_multilevelParams;
    /** Sum of estimator values. */
    double m_sum;
    /** Empirical variance of each estimator sum. */
    double m_varSum;
    /** Empirical biais. */
    double m_biais;
    /** Empirical variance. */
    double m_var;
    /** L2 error. */
    double m_L2Error;
    /** Total number of estimators. */
    int m_totalL;
    /** Total time passed in estimator computation. */
    double m_totalTime;
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
    m_gen(gen), m_f(f), m_scheme(scheme), m_multilevelParams(multilevelParams), m_sum(0.), m_varSum(0.), m_totalL(0), m_totalTime(0.)
{
}

/**
 * @return Computed value of the estimator \f$\bar{Y}_{h,\underline{n}}^{N,q}\f$.
 */
template <typename StateType, typename VolType, typename RandomType, typename TransitionType>
double Estimator<StateType, VolType, RandomType, TransitionType>::compute()
{
    std::chrono::steady_clock::time_point start = std::chrono::steady_clock::now();
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
    MonteCarlo<StateType, VolType, RandomType, TransitionType> monteCarlo(m_gen, m_f, m_scheme, hInverse);
    sum += monteCarlo(N_0);
    double var = monteCarlo.var()/(double)N_0;

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
        DoubleMonteCarlo<StateType, VolType, RandomType, TransitionType> doubleMonteCarlo(m_gen, newF, m_scheme, hInverse*n[j-1], hInverse*n[j]);
        sum += doubleMonteCarlo(N_j);
        var += doubleMonteCarlo.var()/(double)N_j;
    }

    m_sum += sum;
    m_varSum += var;

    std::chrono::steady_clock::time_point  end = std::chrono::steady_clock::now();

    m_totalTime +=  std::chrono::duration_cast<std::chrono::duration<double>>(end-start).count() ;
    return sum;
}


/**
 * @return Compute empirical biais, variance and L2 error of the estimator \f$\bar{Y}_{h,\underline{n}}^{N,q}\f$.
 */
template <typename StateType, typename VolType, typename RandomType, typename TransitionType>
double Estimator<StateType, VolType, RandomType, TransitionType>::L2Error(int L, double trueValue)
{
    for (int i=0; i<L; ++i)
        this->compute();

    double N = m_multilevelParams.getSimulationsNumber();

    m_totalL += L;

    m_biais = m_sum/(double)m_totalL - trueValue;
    m_var = m_varSum/(double)(m_totalL);
    m_L2Error = sqrt(m_biais*m_biais + m_var);

    return m_L2Error;
}

/**
 * @return Display empirical biais, variance and L2 error of the estimator \f$\bar{Y}_{h,\underline{n}}^{N,q}\f$.
 */
template <typename StateType, typename VolType, typename RandomType, typename TransitionType>
void Estimator<StateType, VolType, RandomType, TransitionType>::display()
{
    estimator_type estimatorType= m_multilevelParams.getEstimatorType();
    if (estimatorType == 0)
        cout << "MLMC ESTIMATOR : " << endl;
    else
        cout << "ML2R ESTIMATOR : " << endl;
    cout << "L2 error : " << m_L2Error << endl;
    cout << "Empirical biais : " << m_biais << endl;
    cout << "Empirical variance : " << m_var << endl;
    cout << "Mean value of estimator : " << m_sum/m_totalL << endl;
    cout << "Mean time for one estimator : " << m_totalTime/m_totalL << endl << endl;
}

/**
 * @brief Method to write the estimator result in a file.
 */
template <typename StateType, typename VolType, typename RandomType, typename TransitionType>
void Estimator<StateType, VolType, RandomType, TransitionType>::write(const string fileName)
{
    // Open the file where we want to write the results
    ofstream file_out(fileName.c_str(), ios::out | ios::binary | ios::app);

    file_out.precision(6);

    if (!file_out)  { // if the opening succeeded
        cerr << endl << " Write parameters in file error! " << endl;
        cerr << " Cannot open solution file (" << fileName << ") to write the multilevel parameters." << endl << endl;
        exit(1);
    }

    file_out << endl;
    file_out << "ESTIMATOR : " << endl << endl;

    file_out << "L2 error : " << m_L2Error << endl;
    file_out << "Empirical biais : " << m_biais << endl;
    file_out << "Empirical variance : " << m_var << endl << endl;
    file_out << "Mean value of estimator : " << m_sum/m_totalL << endl;
    file_out << "Mean time for one estimator : " << m_totalTime/m_totalL << endl;
    file_out << endl;
    file_out << "-------------------------------" << endl << endl;
}

#endif // ESTIMATOR_HPP
