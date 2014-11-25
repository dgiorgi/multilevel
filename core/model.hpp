#ifndef MODEL_HPP
#define MODEL_HPP

#include <iostream>
#include <string>
#include <random>
#include <memory>

#include <Eigen/Dense>

using namespace std;

/**
 * @brief The Model class
 *
 * Generic class to define a model described by a SDE.
 *
 * \f[dX_t = b(t, X_t)dt + \sigma(t,X_t)dW_t\f]
 * \f[X_0=x_0\f]
 */
template<typename StateType, typename VolType, typename RandomType> class Model
{
public:
    /** Constructor. */
    Model(const StateType x0, const double tF);

    /** Pure virtual drift function of the model: \f$b(t,x)\f$. */
    virtual StateType drift(const double t, const StateType x)  = 0;
    /** Pure virtual volatility function of the model: \f$\sigma(t,x)\f$. */
    virtual VolType sigma(const double t, const StateType x) = 0;
    /** Pure virtual random variable realization */
    virtual RandomType random(mt19937_64& gen) = 0;

    /** @name Getters
     * @{ */
    /** Starting point getter. */
    StateType getStartingPoint() const {return m_x0;}
    /** Maturity getter. */
    double getMaturity() const {return m_tF;}
    /** @} */

protected:
    /** Starting point \f$x_0\f$. */
    StateType m_x0;
    /** Maturity \f$T\f$. */
    double m_tF;
};

/**
 * @param x0 Starting point.
 * @param tF Maturity.
 */
template<typename StateType, typename VolType, typename RandomType>
Model<StateType, VolType, RandomType>::Model(const StateType x0, const double tF) :
    m_x0(x0), m_tF(tF)
{
}

/**
 * @brief The BlackAndScholes class
 *
 * Black and Scholes model.
 *
 * \f[dX_t = bX_tdt + \sigma X_tdW_t\f]
 * \f[X_0=x_0\f]
 */
template<typename StateType, typename VolType, typename CorrType, typename RandomType>
class BlackAndScholes : public Model<StateType, VolType, RandomType>
{
public:
    /** Constructor */
    BlackAndScholes(const double b,
                    const VolType s,
                    const CorrType rho,
                    const StateType x0,
                    const double tF);

    /** Drift function of the Black and Scholes model: \f$b(t,x)=b*x\f$. */
    StateType drift(const double t, const StateType x){ return m_b*x; }
    /** Volatility function of the Black and Scholes model: \f$\sigma(t,x)=\sigma*x\f$. */
    VolType sigma(const double t, const StateType x){ return m_s*x; }
    /** Random variable realization */
    RandomType random(mt19937_64& gen);

    /** @name Getters
     * @{ */
    /** Starting point getter. */
    std::normal_distribution<double> getRandom() const {return m_gaussian;}
    /** @} */

protected :
    /** Drift constant */
    double m_b;
    /** Volatility matrix */
    VolType m_s;
    /** Correlation matrix */
    CorrType m_rho;
    /** Cholesky decomposition of correlation matrix */
    CorrType m_L;
    /** Generic normal distribution. */
    std::normal_distribution<double> m_gaussian;
};


/**
 * @brief Specialization of the BlackAndScholes class, for the double case.
 */
template<>
class BlackAndScholes<double, double, double, double> : public Model<double, double, double>
{
public:
    /** Constructor */
    BlackAndScholes(const double b,
                    const double s,
                    const double x0,
                    const double tF): Model<double, double, double>(x0, tF), m_b(b), m_s(s), m_gaussian(0.0,1.0){}

    /** Drift function of the Black and Scholes model: \f$b(t,x)=b*x\f$. */
    double drift(const double t, const double x){ return m_b*x; }
    /** Volatility function of the Black and Scholes model: \f$\sigma(t,x)=\sigma*x\f$. */
    double sigma(const double t, const double x){ return m_s*x; }
    /** Random variable realization */
    double random(mt19937_64& gen){ return m_gaussian(gen); }

    /** @name Getters
     * @{ */
    /** Starting point getter. */
    std::normal_distribution<double> getRandom() const {return m_gaussian;}
    /** @} */

protected :
    /** Drift constant */
    double m_b;
    /** Volatility matrix */
    double m_s;
    /** Generic normal distribution. */
    std::normal_distribution<double> m_gaussian;
};



/**
 * @param b Drift constant.
 * @param s Volatility matrix.
 * @param Correlation matrix.
 * @param x0 Starting point.
 * @param tF Maturity.
 */
template<typename StateType, typename VolType, typename CorrType, typename RandomType>
BlackAndScholes<StateType, VolType, CorrType, RandomType>::BlackAndScholes(const double b,
                                                                           const VolType s,
                                                                           const CorrType rho,
                                                                           const StateType x0,
                                                                           const double tF):
    Model<StateType, VolType, RandomType>(x0, tF), m_b(b), m_s(s), m_rho(rho), m_gaussian(0.0,1.0)
{
    // Make the Cholesky decomposition of the correlation matrix
//    string T = typeid(CorrType).name();
//    std::size_t found = T.find("Matrix");

//    if (found!=std::string::npos){
        Eigen::LLT<CorrType> cholesky = m_rho.llt();
        m_L = cholesky.matrixL();
//    }
}

/**
 * @brief Generate the random realization of \f$(dW_1, \ldots, dW_p)\f$ with correlation matrix \f$\rho_{i,j}dt = d\langle W_i, w_j \rangle \f$
 *
 * @param gen Random generator
 * @return randomRealization \f$(dW_1, \ldots, dW_p)\f$
 */
template<typename StateType, typename VolType, typename CorrType, typename RandomType>
RandomType BlackAndScholes<StateType, VolType, CorrType, RandomType>::random(mt19937_64 &gen)
{
    // Generate n realizations of indipendent gaussian laws
    int size = m_L.rows();
    RandomType randomRealization(size);
    for (int i=0; i<size; ++i)
        randomRealization << m_gaussian(gen);

    // Make the product with the matrix L coming from the Cholesky decomposition of the correlation matrix
    randomRealization = m_L*randomRealization;

    return randomRealization;
}


/**
 * @brief The SABR class
 *
 * SABR model.
 *
 * \f[dF_t = \sigma_t F_t^{\beta}dW_t\f]
 * \f[d\sigma_t = \alpha \sigma_t dZ_t\f]
 *
 * with \f$d\langle W_t, Z_t \rangle = \rho dt \f$
 * \f[F_0=f_0\f]
 * \f[Z_0=z_0\f]
 */
class SABR : public Model<Eigen::Vector2d, Eigen::Matrix2d, Eigen::Vector2d>
{
public:
    /** Constructor */
    SABR(const double alpha,
         const double beta,
         const Eigen::Matrix2d sigma,
         const double rho,
         const Eigen::Vector2d x0,
         const double tF);

    /** Drift function of the SABR model: \f$b(t,x)=0\f$. */
    Eigen::Vector2d drift(const double t, const Eigen::Vector2d x){ return Eigen::Vector2d(0,0); }
    /** Volatility function of the SABR model: \f$\sigma(t,x)=\sigma*x\f$. */
    Eigen::Matrix2d sigma(const double t, const Eigen::Vector2d x);
    /** Random variable realization */
    Eigen::Vector2d random(mt19937_64& gen);

protected :
    /** Alpha */
    double m_alpha;
    /** Beta */
    double m_beta;
    /** Volatility constant */
    Eigen::Matrix2d m_s;
    /** Correlation */
    double m_rho;
    /** Cholesky decomposition of correlation matrix */
    Eigen::Matrix2d m_L;
    /** Generic normal distribution. */
    std::normal_distribution<double> m_gaussian;
};

// Smart pointer to Model
template<typename StateType, typename VolType, typename RandomType>
using modelPtr = shared_ptr<Model<StateType, VolType, RandomType>> ;

// Smart pointer to Black And Scholes Model
template<typename StateType, typename VolType, typename CorrType, typename RandomType>
using blackAndScholesfPtr = shared_ptr<BlackAndScholes<StateType, VolType, CorrType, RandomType>> ;

// Smart pointer to SABR Model
using sabrPtr = shared_ptr<SABR> ;

#endif // MODEL_HPP
