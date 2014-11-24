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
 * @brief Method to make a transition of step h, at time t, starting from x and with a random realization random.
 *
 * @param t Starting time.
 * @param x Starting point.
 * @param h Time step.
 * @param r Random realization.
 */
Eigen::VectorXd Euler::transition(double t, Eigen::VectorXd x, double h, Eigen::VectorXd r)
{
    return x + m_model->drift(t,x)*h + m_model->sigma(t, x)*r;
}

/**
 * @brief Method to make a transition for two different discretizations of the same process, over step h, one with discretization 1 and one with discretization n,
 * starting from x at time t.
 *
 * @param gen Random generator.
 * @param t Starting time.
 * @param x Starting point.
 * @param h Time step.
 * @param n Discretization for the finest process.
 */
pair<Eigen::VectorXd, Eigen::VectorXd> Euler::pairTransition(mt19937_64 &gen, double t, pair<Eigen::VectorXd, Eigen::VectorXd> pairProcess, double h, unsigned int n)
{
    // We get the finest step
    double h_n = h/n;

    // We get the starting points
    Eigen::VectorXd X_0 = pairProcess.first;
    Eigen::VectorXd X_n = pairProcess.second;

    Eigen::VectorXd randomRealization;
    Eigen::VectorXd randomRealizationSum;

    // We make n transitions on the finest process, keeping memory of the random realizations
    for (unsigned int i=0; i<n; ++i){
        randomRealization = m_model->random(gen)*sqrt(h_n);
        X_n = transition(t + i*h_n, X_n, h_n, randomRealization);
        randomRealizationSum += randomRealization;
    }

    // We make one transition on the biggest process
    X_0 = transition(t, X_0, h, randomRealizationSum);

    return pair<Eigen::VectorXd, Eigen::VectorXd>(X_0, X_n);
}


/**
 * @brief Method to simulate a realization of \f$X_T\f$.
 *
 * @param gen Generator.
 * @param n Discretization over time.
 * @return Simulated variable \f$X_T\f$.
 */
Eigen::VectorXd Euler::singleSimulation(mt19937_64& gen, const unsigned int n)
{
    double currentTime = 0.;

    Eigen::VectorXd X_i = m_model->getStartingPoint();
    double h = m_model->getMaturity()/(double)n;
    for (unsigned int i=0; i<n; i++){
        Eigen::VectorXd random = sqrt(h)*m_model->random(gen);
        X_i = transition(currentTime, X_i, h, random);
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
pair<Eigen::VectorXd, Eigen::VectorXd> Euler::doubleSimulation(mt19937_64& gen, const unsigned int n1, const unsigned int n2)
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

    unsigned int div = n2/n1;

    Eigen::VectorXd X_n1= m_model->getStartingPoint();
    Eigen::VectorXd X_n2 = m_model->getStartingPoint();

    double h1 = m_model->getMaturity()/(double)n1;

    pair<Eigen::VectorXd, Eigen::VectorXd> pairProcess(X_n1, X_n2);

    // We make a loop over the largest discretization and we make a pair transition over each step
    for (unsigned int i=0; i<n1; ++i){
        pairProcess = pairTransition(gen, currentTime, pairProcess, h1, div);
        currentTime += h1;
    }

    return pairProcess;
}

///**
// * @param model Pointer to the model
// * @param scheme Pointer to the scheme
// */
//PhiScheme::PhiScheme(const modelPtr model, const schemePtr scheme, std::function<double(double)> f): Scheme(model), m_scheme(scheme), m_f(f){}

///**
// * @brief Method to make a transition of step h, at time t, starting from x and with a random realization random.
// *
// * @param t Starting time.
// * @param x Starting point.
// * @param h Time step.
// * @param r Random realization.
// */
//double PhiScheme::transition(double t, double x, double h, double r)
//{
////    return m_model->drift(t,x)*h + m_model->sigma(t, x)*sqrt(h)*m_model->random(gen);
//    return m_model->drift(t,x)*h + m_model->sigma(t, x)*sqrt(h)*r;
//}

///**
// * @brief Method to make a transition for two different discretizations of the same process, over step h, one with discretization 1 and one with discretization n,
// * starting from x at time t.
// *
// * @param gen Random generator.
// * @param t Starting time.
// * @param x Starting point.
// * @param h Time step.
// * @param n Discretization for the finest process.
// */
//pair<double, double> PhiScheme::pairTransition(mt19937_64 &gen, double t, pair<double, double> pairProcess, double h, unsigned int n)
//{
//    // We get the finest step
//    double h_n = h/n;

//    // We get the starting points
//    double X_0 = pairProcess.first;
//    double X_n = pairProcess.second;

//    double randomRealization = 0;
//    double randomRealizationSum = 0;

//    // We make n transitions on the finest process, keeping memory of the random realizations
//    for (unsigned int i=0; i<n; ++i){
//        randomRealization = m_model->random(gen)*sqrt(h_n);
//        X_n = transition(t + i*h_n, X_n, h, randomRealization);
//        randomRealizationSum += randomRealization;
//    }

//    // We make one transition on the biggest process
//    X_0 = transition(t, X_0, h, randomRealizationSum);

//    return pair<double, double>(X_0, X_n);
//}


///**
// * @brief PhiScheme::singleSimulation
// * @param gen
// * @param n
// * @return
// */
//double PhiScheme::singleSimulation(mt19937_64 &gen, const unsigned int n)
//{
////    double currentTime = 0.;

////    double X_i = m_model->getStartingPoint();
////    double h = m_model->getMaturity()/(double)n;

////    for (unsigned int i=0; i<n; i++){
////        X_i += m_model->drift(currentTime,X_i)*h + m_model->sigma(currentTime, X_i)*sqrt(h)*m_model->random(gen);
////        currentTime += h;
////    }

////    return X_i;
//    return 0;
//}

//pair<double, double> PhiScheme::doubleSimulation(mt19937_64 &gen, const unsigned int n1, const unsigned int n2)
//{
//    return pair<double,double>(0,0);
//}
