#ifndef ESTIMATOR_HPP
#define ESTIMATOR_HPP

#include <multilevelparameters.hpp>

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
class Estimator
{
public:
    /** Constructor. */
    Estimator(mt19937_64& gen,
              std::function<double(double)> f,
              const modelfPtr model,
              const MultilevelParameters multilevelParams);
    /** Method to compute the estimator */
    double compute();

protected:
    /** Generator for the random variable. */
    mt19937_64& m_gen;
    /** Function of the random variable that we simulate: \f$ f(X_T)\f$ */
    std::function<double(double)> m_f;
    /** Model for the random variable. For the moment a SDE.*/
    modelfPtr m_model;
    /** Multilevel parameters which define the estimator. */
    MultilevelParameters m_multilevelParams;
};

#endif // ESTIMATOR_HPP
