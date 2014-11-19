#include "scheme.hpp"

/**
 * @param model Pointer to the model.
 */
Scheme::Scheme(const modelPtr model): m_model(model){}




/**
 * @param model Pointer to the model.
 */
Euler::Euler(const modelPtr model): Scheme(model){}

/**
 * @brief Method to simulate a realization of \f$X_T\f$.
 *
 * @param gen Generator.
 * @param n Discretization over time.
 * @return Simulated variable \f$X_T\f$.
 */
double Euler::singleSimulation(mt19937_64& gen, const unsigned int n)
{
    double currentTime = 0.;

    double X_i = m_model->getStartingPoint();
    double h = m_model->getMaturity()/(double)n;
    for (unsigned int i=0; i<n; i++){
        X_i += m_model->drift(currentTime,X_i)*h + m_model->sigma(currentTime, X_i)*sqrt(h)*m_model->random(gen);
        currentTime += h;
    }

    return X_i;
}

/**
 * @brief Method to simulate two realizations of \f$X_T\f$.
 *
 * The two discretizations used are \f$\frac{T-t_0}{n_1}\f$ and \f$\frac{T-t_0}{n_2}\f$.
 *
 * @param gen Generator.
 * @param n1 Discretization over time for first simulation.
 * @param n2 Discretization over time for second simulation.
 * @return Pair of simulations of \f$X_T\f$.
 */
pair<double, double> Euler::doubleSimulation(mt19937_64& gen, const unsigned int n1, const unsigned int n2)
{
    // We throw an error if n2 is not a multiple of n1
    if (n2%n1 != 0){
        cout << "Making a double simulation with two incompatible discretizations." << endl;
        exit(1);
    }

    // We throw an error if n2 is not a bigger than n1
    if (n2<n1){
        cout << "Making a double simulation with two incompatible discretizations, it must be n1 < n2." << endl;
        exit(1);
    }

    double currentTime = 0.;

    unsigned div = n2/n1;

    double X_n1_i = m_model->getStartingPoint();
    double X_n2_i = m_model->getStartingPoint();

    double h1 = m_model->getMaturity()/(double)n1;
    double h2 = m_model->getMaturity()/(double)n2;

    // Variables to store the brownian increments and their sum
    double randomRealization;
    double randomRealizationSum;

    // We simulate two variables using an euler scheme with two different discretizations
    for (unsigned int i=0; i<n1; ++i){
        randomRealizationSum = 0;

        // We first make n2/n1 increments in the euler scheme with the highest discretization
        // Keeping memory of the brownian increments that we used
        for (unsigned int j=0; j<div; ++j){
            randomRealization = m_model->random(gen)*sqrt(h2);
            X_n2_i += m_model->drift(currentTime + j*h2, X_n2_i)*h2 + m_model->sigma(currentTime + j*h2, X_n2_i)*randomRealization;
            randomRealizationSum += randomRealization;
        }
        // Then we make one increment in the euler scheme with the lowest discretization
        // using the sum of the brownian increments used before
        X_n1_i += m_model->drift(currentTime + i*h1, X_n1_i)*h1 + m_model->sigma(currentTime + i*h1, X_n1_i)*randomRealizationSum;
    }

    return pair<double, double>(X_n1_i,X_n2_i);
}
