#include "model.hpp"

#include <chrono>


/**
 * @param alpha \f$\alpha\f$
 * @param beta \f$\beta\f$
 * @param sigma \f$\sigma\f$ Volatility
 * @param rho \f$\rho\f$ Correlation
 * @param x0 \f$x_0\f$ Starting point
 * @param tF \f$T_F\f$ Maturity
 */
SABR::SABR(const double alpha,
           const double beta,
           const Eigen::Matrix2d sigma,
           const double rho,
           const Eigen::Vector2d x0,
           const double tF):
    Model<Eigen::Vector2d, Eigen::Matrix2d, Eigen::Vector2d>(x0, tF), m_alpha(alpha), m_beta(beta), m_s(sigma), m_rho(rho), m_gaussian(0.0,1.0)
{
    // Make the Cholesky decomposition of the correlation matrix
    m_L << 1, 0, m_rho, sqrt(1-rho*rho);
}

/**
 * @param t Current time
 * @param x Current state
 * @return Volatility
 */
Eigen::Matrix2d SABR::sigma(const double t, const Eigen::Vector2d x)
{
    Eigen::Matrix2d vol;
    vol << x[0]*pow(x[1], m_beta), 0, 0, m_alpha*x[1];
    return vol;
}

/**
 * @param gen Random generator
 * @return Random realization of \f$ (dW,dZ)\f$ with correlation \f$\rho\f$
 */
Eigen::Vector2d SABR::random(mt19937_64 &gen)
{
    Eigen::Vector2d gauss(m_gaussian(gen), m_gaussian(gen));
    return m_L*gauss;
}
