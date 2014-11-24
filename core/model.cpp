#include "model.hpp"

#include <chrono>

/**
 * @param x0 Starting point.
 * @param tF Maturity.
 */
Model::Model(const Eigen::VectorXd x0, const double tF) :
    m_x0(x0), m_tF(tF)
{
}


/**
 * @param b Drift constant.
 * @param s Volatility matrix.
 * @param Correlation matrix.
 * @param x0 Starting point.
 * @param tF Maturity.
 */
BlackAndScholes::BlackAndScholes(const double b,
                                 const Eigen::MatrixXd s,
                                 const Eigen::MatrixXd rho,
                                 const Eigen::VectorXd x0,
                                 const double tF):
    Model(x0, tF), m_b(b), m_s(s), m_rho(rho), m_gaussian(0.0,1.0)
{
    // We make the Cholesky decomposition of the correlation matrix
    Eigen::LLT<Eigen::MatrixXd> cholesky = m_rho.llt();
    m_L = cholesky.matrixL();
}

Eigen::VectorXd BlackAndScholes::random(mt19937_64 &gen)
{
    int size = m_x0.rows();
    Eigen::VectorXd randomRealization(size);
    for (int i=0; i<size; ++i)
        randomRealization << m_gaussian(gen);
    return randomRealization;
}
