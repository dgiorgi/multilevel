#ifndef MONTECARLO_HPP
#define MONTECARLO_HPP

#include <functional>
#include <model.hpp>
#include <math.h>

/**
 * @brief The MonteCarlo class
 *
 * Class used to make Monte Carlo simulation \f$\frac{1}{N} \sum_{i=1}^N f(X_T^{(i)})\f$.
 */
class MonteCarlo
{
public:
    /** Constructor. */
    MonteCarlo(mt19937_64& gen,
               std::function<double(double)> f,
               const modelfPtr model,
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
    std::function<double(double)> getF() const {return m_f;}
    /** Model getter. */
    modelfPtr getModel() const {return m_model;}
    /** Model discretization size getter. */
    unsigned int getModelSize() const {return m_modelSize;}
    /** Total number of simulations getter. */
    unsigned int getTotalSize() const {return m_totalN;}
    /** @} */

protected:
    /** Generator for the random variable. */
    mt19937_64& m_gen;
    /** Function of the random variable that we simulate: \f$ f(X_T)\f$ */
    std::function<double(double)> m_f;
    /** Model for the random variable. For the moment a SDE.*/
    modelfPtr m_model;
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
class DoubleMonteCarlo
{
public:
    /** Constructor. */
    DoubleMonteCarlo(mt19937_64& gen,
                     std::function<double(double)> f,
                     const modelfPtr model,
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
    std::function<double(double)> getF() const {return m_f;}
    /** Model getter. */
    modelfPtr getModel() const {return m_model;}
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
    std::function<double(double)> m_f;
    /** Model for the random variable. For the moment a SDE.*/
    modelfPtr m_model;
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

#endif // MONTECARLO_HPP
