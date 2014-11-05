#include "montecarlo.hpp"

/**
 * @param gen Generator for the random variable
 * @param f Function of the random variable that we simulate \f$f(X_T)\f$
 * @param model Model for \f$X_t\f$
 * @param modelSize Discretization size of the SDE
 */
MonteCarlo::MonteCarlo(mt19937_64& gen,
                       std::function<double(double)> f,
                       const modelfPtr model,
                       const unsigned int modelSize):
    m_gen(gen), m_f(f), m_model(model), m_modelSize(modelSize), m_totalN(0), m_sum(0.), m_sumSquares(0.)
{
}

/**
 * @param N Number of simulations done.
 * @return Empirical mean of the random variable, computed with the N simulations and the previous calls.
 */
double MonteCarlo::operator()(const unsigned int N){

    // We make N calls to the simulator,
    // we make the sum and the sum of squares and add them to the previously computed ones.
    for (unsigned int i=0; i< N; ++i) {
        double x = m_model->singleSimulation(m_gen, m_modelSize);
        m_sum += m_f(x);
        m_sumSquares += m_f(x)*m_f(x);
    }
    // We actualize the total number of simulations
    m_totalN += N;

    return mean();
};


/**
 * @param gen Generator for the random variable
 * @param f Function of the random variable that we simulate \f$f(X_T)\f$
 * @param model Model for \f$X_t\f$
 * @param modelSize1 First discretization size of the SDE
 * @param modelSize2 Second discretization size of the SDE
 */
DoubleMonteCarlo::DoubleMonteCarlo(mt19937_64& gen,
                                   std::function<double(double)> f,
                                   const modelfPtr model,
                                   const unsigned int modelSize1,
                                   const unsigned int modelSize2):
    m_gen(gen), m_f(f), m_model(model), m_modelSize1(modelSize1), m_modelSize2(modelSize2), m_totalN(0), m_sum(0.), m_sumSquares(0.)
{
}

/**
 * @param N Number of simulations done.
 * @return Empirical mean of the difference between the two random variables,
 * computed with the N simulations and the previous calls.
 */
double DoubleMonteCarlo::operator()(const unsigned int N){

    // We make N calls to the simulator,
    // we make the sum and the sum of squares and add them to the previously computed ones.
    for (unsigned int i=0; i< N; ++i) {
        pair<double, double> x = m_model->doubleSimulation(m_gen, m_modelSize1, m_modelSize2 );
        double diff = m_f(x.second) - m_f(x.first) ;
        m_sum += diff;
        m_sumSquares += diff*diff;
    }
    // We actualize the total number of simulations
    m_totalN += N;

    return mean();
};
