#ifndef STRUCTURALPARAMETERS_HPP
#define STRUCTURALPARAMETERS_HPP

#include <math.h>
#include <iostream>
#include <fstream>

#include "montecarlo.hpp"

/**
 * @brief The StructuralParameters class
 *
 * We consider a family \f$(Y_h)_{h\in \mathcal{H}}\f$ of real-valued random variables associated
 * to a random variable \f$Y_0\f$, indexed by \f$\mathcal{H}\subset (0,\mathbf{h}]\f$.
 * The family satisfies the strong and the weak rates of approximation of \f$Y_0\f$ by \f$Y_h\f$
 * when \f$h\to 0\f$.
 *
 * <b>Bias error expansion (weak error rate)</b>
 *
 * \f$\exists \alpha > 0, R\geq 1, \; \mathbf{E}(Y_h)= \mathbf{E}(Y_0) + \sum_{k=1}^R c_k h^{\alpha k}
 * + h^{\alpha R} \eta_R(h), \; lim_{h\to 0}\eta_R(h) = 0\f$
 *
 * <b>Strong approximation error assumption</b>
 *
 * \f$\exists \beta > 0, \; \|Y_h-Y_0\|_2^2 = \mathbf{E}(|Y_h-Y_0|^2) \leq V_1 h^\beta\f$
 *
 * We set \f[\theta = \sqrt{\frac{V_1}{var(Y_0)}}\f]
 *
 * This class describes the structural parameters related to these approximations.
 *
 */
class StructuralParameters
{
public:
    StructuralParameters(const double alpha,
                         const double beta,
                         const double c1,
                         const double h);

    /** @name Compute parameters functions*/
    template <typename StateType, typename VolType, typename RandomType>
    void computeParameters(mt19937_64& gen,
                           std::function<double(StateType)> f,
                           const schemePtr<StateType, VolType, RandomType> scheme,
                           const unsigned int N);
    template <typename StateType, typename VolType, typename RandomType>
    void computeV1(mt19937_64& gen,
                   std::function<double(StateType)> f,
                   const schemePtr<StateType, VolType, RandomType> scheme,
                   const unsigned int N);
    template <typename StateType, typename VolType, typename RandomType>
    void computeVarY0(mt19937_64& gen,
                      std::function<double(StateType)> f,
                      const schemePtr<StateType, VolType, RandomType> scheme,
                      const unsigned int N);
    void computeTheta();
    void displayParameters();
    void writeParameters(const string fileName);

    /** @name Getters
     * @{ */
    /** \f$\alpha\f$ getter. */
    double getAlpha() const {return m_alpha;}
    /** \f$\beta\f$ getter. */
    double getBeta() const {return m_beta;}
    /** \f$V_1\f$ getter. */
    double getV1() const {return m_V1;}
    /** \f$c_1\f$ getter. */
    double getC1() const {return m_c1;}
    /** \f$var(Y_0)\f$ getter. */
    double getVarY0() const {return m_varY0;}
    /** \f$h\f$ getter. */
    double getHBold() const {return m_hBold;}
    /** \f$\vartheta\f$ getter. */
    double getTheta() const {return m_theta;}
    /** @} */

private:
    /** @name Given to the constructor
      * @{ */
    /** \f$\alpha\f$ */
    double m_alpha;
    /** \f$\beta\f$ */
    double m_beta;
    /** \f$c_1\f$ */
    double m_c1;
    /** \f$\mathbf{h}\f$ */
    double m_hBold;
    /** @} */

    /** @name Computed
      * @{ */
    /** \f$V_1\f$ */
    double m_V1;
    /** \f$var(Y_0)\f$ */
    double m_varY0;
    /** \f$\vartheta\f$ */
    double m_theta;
    /** @} */
};


/**
 * @brief Generic method to compute all the structural parameters
 * \f$V_1, var(Y_0)\f$ and \f$\theta\f$.
 *
 * @param gen Generator for the random variable.
 * @param f Function \f$f\f$ of the random variable.
 * @param scheme Scheme for the simulation.
 * @param N Number of simulations.
 */
template <typename StateType, typename VolType, typename RandomType>
void StructuralParameters::computeParameters(mt19937_64& gen,
                                             std::function<double(StateType)> f,
                                             const schemePtr<StateType, VolType, RandomType> scheme,
                                             const unsigned int N)
{
    // First we compute the var(Y0) and V1
    computeVarY0(gen, f, scheme, N);
    computeV1(gen, f, scheme, N);

    // And then theta
    computeTheta();

    // We display the parameters
    displayParameters();
    writeParameters("parameters.txt");
}


/**
 * @brief Method to compute the parameter \f$V_1\f$. We consider the following estimator:
 * \f[V_1(h) = (1 + M_{max}^{\frac{\beta}{2}})^{-2}h^{-\beta}\|Y_h-Y_{h/M_{max}}\|_2^2\f]
 * where we set \f$M_{max} = 10\f$ and \f$h = \mathbf{h}.\f$.
 *
 * @param gen Generator for the random variable.
 * @param f Function \f$f\f$ of the random variable.
 * @param scheme Scheme for the simulation.
 * @param N Number of simulations.
 */
template <typename StateType, typename VolType, typename RandomType>
void StructuralParameters::computeV1(mt19937_64& gen,
                                     std::function<double(StateType)> f,
                                     const schemePtr<StateType, VolType, RandomType> scheme,
                                     const unsigned int N)
{
    unsigned int M = 10;
    double T = scheme->getModel()->getMaturity();


    // Throw error if the rest of T/m_hBold is not 0
    double epsilon = 1e-5;
    if (fabs(T/m_hBold-(int)(T/m_hBold)) > epsilon){
        cerr << "The discretization step is not an integer divisor of the maturity. Please choose a good discretization step." << endl;
        exit(1);
    }

    int n = T/m_hBold;
    double beta = m_beta;

    // We instanciate a double monte carlo
    DoubleMonteCarlo<StateType, VolType, RandomType> Y(gen, f, scheme, n, M*n);
    // We make N simulations
    Y(N);

    double L2 = Y.meanOfSquares();

    m_V1 = pow(1.0+pow((double)M, -beta*0.5), -2) * pow(m_hBold,-beta) * L2;
}

/**
 * @brief Method to compute the variance \f$var(Y_0) \sim var(Y_{\mathbf{h}})\f$.
 *
 * We take here the easiest choice h=min(1,T) i.e. 1 discretization.
 *
 * @param gen Generator for the random variable.
 * @param f Function \f$f\f$ of the random variable.
 * @param scheme Scheme for the simulation.
 * @param N Number of simulations.
 */
template <typename StateType, typename VolType, typename RandomType>
void StructuralParameters::computeVarY0(mt19937_64& gen,
                                        std::function<double(StateType)> f,
                                        const schemePtr<StateType, VolType, RandomType> scheme,
                                        const unsigned int N)
{
    double T = scheme->getModel()->getMaturity();

    // Throw error if the rest of T/m_hBold is not 0
    double epsilon = 1e-5;
    if (fabs(T/m_hBold-(int)(T/m_hBold)) > epsilon){
        cerr << "The discretization step is not an integer divisor of the maturity. Please choose a good discretization step." << endl;
        exit(1);
    }

    int n = T/m_hBold;

    // We instanciate a single Monte Carlo
    MonteCarlo<StateType, VolType, RandomType> Y(gen, f, scheme, n);
    // We make N simulations
    Y(N);
    m_varY0 = Y.var();
}

#endif // STRUCTURALPARAMETERS_HPP
