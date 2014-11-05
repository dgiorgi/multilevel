#ifndef MODEL_HPP
#define MODEL_HPP

#include <iostream>
#include <string>
#include <random>
#include <memory>

using namespace std;

/**
 * @brief The Model class
 *
 * Generic class to define a model described by a SDE.
 *
 * \f[dX_t = b(t, X_t)dt + \sigma(t,X_t)dW_t\f]
 * \f[X_0=x_0\f]
 */
class Model
{
public:
    /** Constructor. */
    Model(const double x0, const double tF);

    /** Virtual function to simulate a process. */
    virtual double singleSimulation(mt19937_64& gen, const unsigned int n) = 0;

    /** Virtual function to simulate two processes. */
    virtual pair<double, double> doubleSimulation(mt19937_64& gen, const unsigned int n1, const unsigned int n2) = 0;

    /** @name Getters
     * @{ */
    /** Starting point getter. */
    double getStartingPoint() const {return m_x0;}
    /** Maturity getter. */
    double getMaturity() const {return m_tF;}
    /** @} */

protected:
    /** Starting point \f$x_0\f$. */
    double m_x0;
    /** Maturity \f$T\f$. */
    double m_tF; // This is set to 1.
    /** Generic normal distribution. */
    std::normal_distribution<double> m_Gaussian;
};

/**
 * @brief The BlackAndScholes class
 *
 * Black and Scholes model.
 *
 * \f[dX_t = bX_tdt + \sigma X_tdW_t\f]
 * \f[X_0=x_0\f]
 */
class BlackAndScholes : public Model
{
public:
    /** Constructor */
    BlackAndScholes(const double b,
                    const double s,
                    const double x0,
                    const double tF);

    /** Drift function of the Black and Scholes model: \f$b(t,x)=b*x\f$. */
    double drift(const double t, const double x){ return m_b*x; }
    /** Volatility function of the Black and Scholes model: \f$\sigma(t,x)=\sigma*x\f$. */
    double sigma(const double t, const double x){ return m_s*x; }

    double singleSimulation(mt19937_64& gen, const unsigned int n);
    pair<double, double> doubleSimulation(mt19937_64& gen, const unsigned int n1, const unsigned int n2);

protected :
    /** Drift constant */
    double m_b;
    /** Volatility constant */
    double m_s;
};

// Smart pointers to Model and BlackAndScholes objects
typedef shared_ptr<Model> modelfPtr;
typedef shared_ptr<BlackAndScholes> blackAndScholesfPtr;

#endif // MODEL_HPP
