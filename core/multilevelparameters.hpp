#ifndef MULTILEVELPARAMETERS_HPP
#define MULTILEVELPARAMETERS_HPP

#include <vector>
#include <valarray>
#include <math.h>
#include <iostream>
#include <fstream>

#include "structuralparameters.hpp"

const double MY_INF = 2e19;

enum estimator_type { MC, RR};

/**
 * @brief The Refiners class
 *
 * We consider here the case \f$n_1=1\f$ and \f$n_i=M^{(i-1)}\f$ for \f$i=2, \ldots,R\f$.
 */
class Refiners {
public:
    /** Constructor.*/
    Refiners(unsigned int R = 2, unsigned int M = 2) : m_M(M), m_data(R) {
        m_data[0] = 1;
        for (unsigned int i = 1; i < R; ++i)
            m_data[i] = m_data[i-1] * m_M;
    }
    /** Operator []. */
    unsigned int operator[](unsigned int i) const { return m_data[i]; }
    /** Root of the refiner getter. */
    unsigned int getM() const { return m_M; }
    /** Clean method */
    void clear();
private:
    /** Root of the refiner. */
    unsigned int m_M;
    /** List of refiners in ascending order. */
    std::vector<unsigned int> m_data;
};

/**
 * @brief The MultilevelParameters class
 *
 * This class summarizes the asymptotic optimal value of the parameters \f$q, R, h\f$ and \f$n\f$
 * for the Multilevel Richardson-Romberg estimator (MLRR) and the Multilevel Monte Carlo estimator
 * (MLMC), for a fixed \f$\varepsilon > 0\f$.
 */
class MultilevelParameters
{
public:
    /** Constructor */
    MultilevelParameters(const double epsilon,
                         const StructuralParameters& structParam,
                         const estimator_type type);

    /** @name Methods to compute the parameters
     * @{
     */
    void initialize();
    void computeOptimalParametersForM(const unsigned M);
    void computeOptimalParameters();
    void computeOrder();
    void computeBiais();
    void computeWeights();
    void computeStratification();
    void computeN();
    void displayParameters();
    void writeParameters(const string fileName);
    /** @} */

    /** Method to compute the cost of the estimator **/
    double computeCost();

    /** @name Getters
     * @{
     */
    /** Root \f$M\f$ getter. */
    unsigned int getRoot()const {return m_M;}
    /** Order \f$R\f$ getter. */
    unsigned int getOrder() const {return m_R;}
    /** Biais \f$h\f$ getter. */
    double getBiais() const {return m_h;}
    /** Inverse of biais \f$h\f$ getter. */
    unsigned int getBiaisInverse() const {return m_hInverse;}
    /** Stratification \f$q=(q_1, \ldots, q_R)\f$ getter. */
    std::vector<double> getStratification() const {return m_q;}
    /** Total number of simulations \f$N\f$ getter. */
    double getSimulationsNumber() const {return m_N;}
    /** Allocation matrix weights \f$W_i\f$ getter. */
    std::vector<double> getWeights() const {return m_W;}
    /** Precision \f$\varepsilon\f$ getter. */
    double getPrecision() const {return m_epsilon;}
    /** Refiners \f$n=(n_1,\ldots,n_R)\f$ getter. */
    Refiners getRefiners() const {return m_n;}
    /** Structural parameters getter. */
    StructuralParameters getStructParameters() const {return m_structParam;}
    /** Estimator type getter. */
    estimator_type getEstimatorType() const {return m_type;}
    /** @} */

protected:
    /** Root: \f$M\f$. */
    unsigned int m_M;
    /** Order of the estimator: \f$R\f$. */
    unsigned int m_R;
    /** Biais parameter: \f$h\f$. */
    double m_h;
    /** Inverse of biais parameter: \f$h^{-1}\f$. */
    unsigned int m_hInverse;
    /** Stratification strategy: \f$q=(q_1, \ldots, q_R)\f$ with \f$q_i>0\f$ and \f$\sum_jq_j=1\f$. */
    std::vector<double> m_q;
    /** Number of simulations: \f$N\f$. */
    double m_N;
    /** Cost of the estimator */
    double m_cost;
    /** Weights of the allocation matrix (computed only in the RR case).*/
    std::vector<double> m_W;
    /** Precision of the estimator: \f$\varepsilon\f$.*/
    double m_epsilon;
    /** Refiners: \f$n=(n_1,\ldots,n_R)\f$ such that \f$n_1=1\f$ and \f$n_i=M^{(i-1)}\f$ for \f$i=2, \ldots,R\f$. */
    Refiners m_n;
    /** Structural parameters in the formulas for weak and strong error. See StructuralParameters class.*/
    StructuralParameters m_structParam;
    /** Type of multilevel estimator, can be Richardson-Romberg (RR) or Monte Carlo (MC). */
    estimator_type m_type; // can be RR or MC
    /** Flag for computed parameters */
    bool m_parametersComputationDone;
};

#endif // MULTILEVELPARAMETERS_HPP
