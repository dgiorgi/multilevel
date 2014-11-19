#include "model.hpp"

#include <chrono>

/**
 * @param x0 Starting point.
 * @param tF Maturity.
 */
Model::Model(const double x0, const double tF) :
    m_x0(x0), m_tF(tF)
{
}


/**
 * @param b Drift constant.
 * @param s Volatility constant.
 * @param x0 Starting point.
 * @param tF Final time
 */
BlackAndScholes::BlackAndScholes(const double b,
                                 const double s,
                                 const double x0,
                                 const double tF):
    Model(x0, tF), m_b(b), m_s(s), m_gaussian(0.0,1.0)
{
}
