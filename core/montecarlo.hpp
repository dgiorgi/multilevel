#ifndef MONTECARLO_HPP
#define MONTECARLO_HPP

#include <functional>
#include <math.h>
#include <omp.h>

#include "model.hpp"
#include "scheme.hpp"

/**
 * @brief The MonteCarlo class
 *
 * Class used to make Monte Carlo simulation \f$\frac{1}{N} \sum_{i=1}^N f(X_T^{(i)})\f$.
 */
template<typename StateType, typename VolType, typename RandomType, typename TransitionType>
class MonteCarlo
{
public:
    /** Constructor. */
    MonteCarlo(mt19937_64& gen,
               std::function<double(TransitionType)> f,
               const schemePtr<StateType, VolType, RandomType, TransitionType> scheme,
               const unsigned int modelSize);

    /** Method to compute the empirical mean \f$\frac 1 N \sum_{i=1}^N f(X_T^{(i)})\f$. */
    double mean() const { return m_sum/(double)m_totalN; }
    /** Method to compute the empirical mean of squares \f$\frac 1 N \sum_{i=1}^N f(X_T^{(i)})^2\f$. */
    double meanOfSquares() const { return m_sumSquares/(double)m_totalN; }
    /** Method to compute the empirical variance \f$\frac{\sum_{i=1}^N f(X_T^{(i)})^2-(\sum_{i=1}^N f(X_T^{(i)})^2/N}{N-1}\f$. */
    double var() const { return (m_sumSquares - m_sum*m_sum/(double)m_totalN)/((double)m_totalN-1.0); }
    /** Method to compute the standard deviation. */
    double stDev() const { return sqrt(var()); }

    /** Method to compute the Monte Carlo simulation.
     *
     * It makes the sum of \f$N\f$ new simulations (of the random variable and of its squares),
     * it adds this quantity to the previously computed sum (if the method was called at least once),
     * it stores the sum and the sum of the squares and returns the mean. */
    double operator()(const unsigned int N);

    /** @name Getters
     * @{ */
    /** Generator getter. */
    mt19937_64& getGenerator() const {return m_gen;}
    /** Function \f$f\f$ getter. */
    std::function<double(TransitionType)> getF() const {return m_f;}
    /** Scheme getter. */
    schemePtr<StateType, VolType, RandomType, TransitionType> getScheme() const {return m_scheme;}
    /** Model discretization size getter. */
    unsigned int getModelSize() const {return m_modelSize;}
    /** Total number of simulations getter. */
    unsigned int getTotalSize() const {return m_totalN;}
    /** @} */

protected:
    /** Generator for the random variable. */
    mt19937_64& m_gen;
    /** Function of the random variable that we simulate: \f$ f(X_T)\f$ */
    std::function<double(TransitionType)> m_f;
    /** Scheme for the simulation.*/
    schemePtr<StateType, VolType, RandomType, TransitionType> m_scheme;
    /** Discretization size of the random variable. */
    unsigned int m_modelSize;
    /** Total number of simulations. It's incremented every time we call of the simulation method. */
    unsigned int m_totalN;

    /** Total sum of the simulations. For one call this is equal to \f$\sum_{i=1}^N f(X_T^{(i)})\f$.*/
    double m_sum;
    /** Total sum of the squares of the simulations. For one call this is equal to \f$\sum_{i=1}^N {f(X_T^{(i)})}^2\f$.*/
    double m_sumSquares;
};

/**
 * @brief The DoubleMonteCarlo class
 *
 * Class used to compute \f$\mathbf{E}(Y_{\frac{h}{n1}} - Y_{\frac{h}{n2}})\f$.
 */
template<typename StateType, typename VolType, typename RandomType, typename TransitionType>
class DoubleMonteCarlo
{
public:
    /** Constructor. */
    DoubleMonteCarlo(mt19937_64& gen,
                     std::function<double(TransitionType)> f,
                     const schemePtr<StateType, VolType, RandomType, TransitionType> scheme,
                     const unsigned int modelDisc1,
                     const unsigned int modelDisc2);

    /** Method to compute the empirical mean \f$\frac 1 N \sum_{i=1}^N f(X_T^{(i)})\f$. */
    double mean() const { return m_sum/(double)m_totalN; }
    /** Method to compute the empirical mean of squares \f$\frac 1 N \sum_{i=1}^N f(X_T^{(i)})^2\f$. */
    double meanOfSquares() const { return m_sumSquares/(double)m_totalN; }
    /** Method to compute the empirical variance \f$\frac{\sum_{i=1}^N f(X_T^{(i)})^2-(\sum_{i=1}^N f(X_T^{(i)})^2/N}{N-1}\f$. */
    double var() const { return (m_sumSquares - m_sum*m_sum/(double)m_totalN)/((double)m_totalN-1); }
    /** Method to compute the standard deviation. */
    double stDev() const { return sqrt(var()); }

    /** Method to compute the Double Monte Carlo simulation.
     *
     * It makes the sum of \f$N\f$ new simulations
     * (of the difference of the random variables \f$Y_{\frac{h}{n1}} - Y_{\frac{h}{n2}}\f$ and of its squares),
     * it adds this quantity to the previously computed sum (if the method was called at least once),
     * it stores the sum and the sum of the squares and returns the mean. */
    double operator()(const unsigned int N);

    /** @name Getters
     * @{ */
    /** Generator getter. */
    mt19937_64& getGenerator() const {return m_gen;}
    /** Function \f$f\f$ getter. */
    std::function<double(TransitionType)> getF() const {return m_f;}
    /** Scheme getter. */
    schemePtr<StateType, VolType, RandomType, TransitionType> getScheme() const {return m_scheme;}
    /** First discretization size getter. */
    unsigned int getModelSize1() const {return m_modelSize1;}
    /** Second discretization size getter. */
    unsigned int getModelSize2() const {return m_modelSize2;}
    /** Total number of simulations getter. */
    unsigned int getTotalSize() const {return m_totalN;}
    /** @} */

protected:
    /** Generator for the random variable. */
    mt19937_64& m_gen;
    /** Function of the random variable that we simulate: \f$ f(X_T)\f$ */
    std::function<double(TransitionType)> m_f;
    /** Scheme for the simulation.*/
    schemePtr<StateType, VolType, RandomType, TransitionType> m_scheme;
    /** First discretization size of the random variable. */
    unsigned int m_modelSize1;
    /** Second discretization size of the random variable. */
    unsigned int m_modelSize2;
    /** Total number of simulations. It's incremented every time we call of the simulation method. */
    unsigned int m_totalN;

    /** Total sum of the simulations. For one call this is equal to \f$\sum_{i=1}^N f(X_T^{(i)})\f$.*/
    double m_sum;
    /** Total sum of the squares of the simulations. For one call this is equal to \f$\sum_{i=1}^N {f(X_T^{(i)})}^2\f$.*/
    double m_sumSquares;
};

/**
 * @param gen Generator for the random variable
 * @param f Function of the random variable that we simulate \f$f(X_T)\f$
 * @param scheme Scheme for the simulation of \f$X_t\f$
 * @param modelSize Discretization size of the SDE
 */
template<typename StateType, typename VolType, typename RandomType, typename TransitionType>
MonteCarlo<StateType, VolType, RandomType, TransitionType>::MonteCarlo(mt19937_64& gen,
                                                       std::function<double(TransitionType)> f,
                                                       const schemePtr<StateType, VolType, RandomType, TransitionType> scheme,
                                                       const unsigned int modelSize):
    m_gen(gen), m_f(f), m_scheme(scheme), m_modelSize(modelSize), m_totalN(0), m_sum(0.), m_sumSquares(0.)
{
}

/**
 * @param N Number of simulations done.
 * @return Empirical mean of the random variable, computed with the N simulations and the previous calls.
 */
template<typename StateType, typename VolType, typename RandomType, typename TransitionType>
double MonteCarlo<StateType, VolType, RandomType, TransitionType>::operator()(const unsigned int N){

    // We make N calls to the simulator,
    // we make the sum and the sum of squares and add them to the previously computed ones.
    for (unsigned int i=0; i< N; ++i) {
        TransitionType x = m_scheme->singleSimulation(m_gen, m_modelSize);
        m_sum += m_f(x);
        m_sumSquares += m_f(x)*m_f(x);
    }
    // We actualize the total number of simulations
    m_totalN += N;

    return mean();
}


/**
 * @param gen Generator for the random variable
 * @param f Function of the random variable that we simulate \f$f(X_T)\f$
 * @param scheme Scheme for the simulation of \f$X_t\f$
 * @param modelSize1 First discretization size of the SDE
 * @param modelSize2 Second discretization size of the SDE
 */
template<typename StateType, typename VolType, typename RandomType, typename TransitionType>
DoubleMonteCarlo<StateType, VolType, RandomType, TransitionType>::DoubleMonteCarlo(mt19937_64& gen,
                                                                   std::function<double(TransitionType)> f,
                                                                   const schemePtr<StateType, VolType, RandomType, TransitionType> scheme,
                                                                   const unsigned int modelSize1,
                                                                   const unsigned int modelSize2):
    m_gen(gen), m_f(f), m_scheme(scheme), m_modelSize1(modelSize1), m_modelSize2(modelSize2), m_totalN(0), m_sum(0.), m_sumSquares(0.)
{
}

/**
 * @param N Number of simulations done.
 * @return Empirical mean of the difference between the two random variables,
 * computed with the N simulations and the previous calls.
 */
template<typename StateType, typename VolType, typename RandomType, typename TransitionType>
double DoubleMonteCarlo<StateType, VolType, RandomType, TransitionType>::operator()(const unsigned int N){

    // We make N calls to the simulator,
    // we make the sum and the sum of squares and add them to the previously computed ones.

//#pragma omp parallel for       //making compile errors
    for (unsigned int i=0; i< N; ++i) {
        pair<TransitionType, TransitionType> x = m_scheme->doubleSimulation(m_gen, m_modelSize1, m_modelSize2 );
        double diff = m_f(x.second) - m_f(x.first) ;
        m_sum += diff;
        m_sumSquares += diff*diff;
    }
    // We actualize the total number of simulations
    m_totalN += N;

    return mean();
}


#endif // MONTECARLO_HPP
