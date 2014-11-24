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
    // Make the Cholesky decomposition of the correlation matrix
    Eigen::LLT<Eigen::MatrixXd> cholesky = m_rho.llt();
    m_L = cholesky.matrixL();
}

/**
 * @brief Generate the random realization of \f$(dW_1, \ldots, dW_p)\f$ with correlation matrix \f$\rho_{i,j}dt = d\langle W_i, w_j \rangle \f$
 *
 * @param gen Random generator
 * @return randomRealization \f$(dW_1, \ldots, dW_p)\f$
 */
Eigen::VectorXd BlackAndScholes::random(mt19937_64 &gen)
{
    // Generate n realizations of indipendent gaussian laws
    int size = m_L.rows();
    Eigen::VectorXd randomRealization(size);
    for (int i=0; i<size; ++i)
        randomRealization << m_gaussian(gen);

    // Make the product with the matrix L coming from the Cholesky decomposition of the correlation matrix
    randomRealization = m_L*randomRealization;

    return randomRealization;
}
