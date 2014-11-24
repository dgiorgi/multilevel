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
    void computeParameters(mt19937_64& gen,
                           std::function<double(Eigen::VectorXd)> f,
                           const schemePtr scheme,
                           const unsigned int N);
    void computeV1(mt19937_64& gen,
                   std::function<double(Eigen::VectorXd)> f,
                   const schemePtr scheme,
                   const unsigned int N);
    void computeVarY0(mt19937_64& gen,
                      std::function<double(Eigen::VectorXd)> f,
                      const schemePtr scheme,
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

#endif // STRUCTURALPARAMETERS_HPP
