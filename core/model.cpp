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
    Model(x0, tF), m_b(b), m_s(s), m_Gaussian(0.0,1.0)
{
}

///**
// * @brief Method to simulate a realization of \f$X_T\f$.
// *
// * @param gen Generator.
// * @param n Discretization over time.
// * @return Simulated variable \f$X_T\f$.
// */
//double BlackAndScholes::singleSimulation(mt19937_64& gen, const unsigned int n)
//{
//    double X_i = m_x0;
//    double h = m_tF/(double)n;
//    for (unsigned int i=0; i<n; i++){
//        X_i += m_b*X_i*h + m_s*X_i*sqrt(h)*m_Gaussian(gen);
//    }

//    return X_i;
//}

///**
// * @brief Method to simulate two realizations of \f$X_T\f$.
// *
// * The two discretizations used are \f$\frac{T-t_0}{n_1}\f$ and \f$\frac{T-t_0}{n_2}\f$.
// *
// * @param gen Generator.
// * @param n1 Discretization over time for first simulation.
// * @param n2 Discretization over time for second simulation.
// * @return Pair of simulations of \f$X_T\f$.
// */
//pair<double, double> BlackAndScholes::doubleSimulation(mt19937_64& gen, const unsigned int n1, const unsigned int n2)
//{
//    // We throw an error if n2 is not a multiple of n1
//    if (n2%n1 != 0){
//        cout << "Making a double simulation with two incompatible discretizations." << endl;
//        exit(1);
//    }

//    // We throw an error if n2 is not a bigger than n1
//    if (n2<n1){
//        cout << "Making a double simulation with two incompatible discretizations, it must be n1 < n2." << endl;
//        exit(1);
//    }

//    unsigned div = n2/n1;

//    double X_n1_i = m_x0;
//    double X_n2_i = m_x0;

//    double h1 = m_tF/(double)n1;
//    double h2 = m_tF/(double)n2;

//    // Variables to store the brownian increments and their sum
//    double gauss;
//    double gaussSum;

//    // We simulate two variables using an euler scheme with two different discretizations
//    for (unsigned int i=0; i<n1; ++i){
//        gaussSum = 0;
//        // We first make n2/n1 increments in the euler scheme with the highest discretization
//        // Keeping memory of the brownian increments that we used
//        for (unsigned int j=0; j<div; ++j){
//            gauss = m_Gaussian(gen)*sqrt(h2);
//            X_n2_i += m_b*X_n2_i*h2 + m_s*X_n2_i*gauss;
//            gaussSum += gauss;
//        }
//        // Then we make one increment in the euler scheme with the lowest discretization
//        // using the sum of the brownian increments used before
//        X_n1_i += m_b*X_n1_i*h1 + m_s*X_n1_i*gaussSum;
//    }

//    return pair<double, double>(X_n1_i,X_n2_i);
//}
