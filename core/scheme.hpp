#ifndef SCHEME_HPP
#define SCHEME_HPP

#include "model.hpp"

/**
 * @brief The Scheme class
 *
 * Generic class to define a scheme on the model that will be used to simulate the process.
 *
 */
template<typename StateType, typename VolType, typename RandomType>
class Scheme
{
public:
    /** Constructor. */
    Scheme(const modelPtr<StateType, VolType, RandomType> model);

    /** Pure virtual function to make a transion. */
    virtual StateType transition(double t, StateType x, double h, RandomType random) = 0;

    /** Pure virtual function to make a coupled transion. */
    virtual pair<StateType,StateType> pairTransition(mt19937_64& gen, double t, pair<StateType, StateType> x, double h, unsigned int n) = 0;

    /** Pure virtual function to simulate a process. */
    virtual StateType singleSimulation(mt19937_64& gen, const unsigned int n) = 0;

    /** Pure virtual function to simulate two processes. */
    virtual pair<StateType, StateType> doubleSimulation(mt19937_64& gen, const unsigned int n1, const unsigned int n2) = 0;

    modelPtr<StateType, VolType, RandomType> getModel() const {return m_model;}

protected:
    modelPtr<StateType, VolType, RandomType> m_model;
};


/**
 * @brief The Euler class
 *
 * Generic class to define an Euler scheme on the model that will be used to simulate the process.
 *
 */
template<typename StateType, typename VolType, typename RandomType>
class Euler : public Scheme<StateType, VolType, RandomType>
{
public:
    /** Constructor. */
    Euler(const modelPtr<StateType, VolType, RandomType> model);

    StateType transition(double t, StateType x, double h, RandomType random);
    pair<StateType,StateType> pairTransition(mt19937_64& gen, double t, pair<StateType, StateType> x, double h, unsigned int n);

    StateType singleSimulation(mt19937_64& gen, const unsigned int n);
    pair<StateType, StateType> doubleSimulation(mt19937_64& gen, const unsigned int n1, const unsigned int n2);
};



/**
 * @param model Pointer to the model.
 */
template<typename StateType, typename VolType, typename RandomType>
Scheme<StateType, VolType, RandomType>::Scheme(const modelPtr<StateType, VolType, RandomType> model): m_model(model){}


/**
 * @param model Pointer to the model.
 */
template<typename StateType, typename VolType, typename RandomType>
Euler<StateType, VolType, RandomType>::Euler(const modelPtr<StateType, VolType, RandomType> model): Scheme<StateType, VolType, RandomType>(model){}


/**
 * @brief Method to make a transition of step h, at time t, starting from x and with a random realization random.
 *
 * @param t Starting time.
 * @param x Starting point.
 * @param h Time step.
 * @param r Random realization.
 */
template<typename StateType, typename VolType, typename RandomType>
StateType Euler<StateType, VolType, RandomType>::transition(double t, StateType x, double h, RandomType r)
{
    return x + this->m_model->drift(t,x)*h + this->m_model->sigma(t, x)*r;  // on ne passe pas sqrt(h) ici à cause de la double transition
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
template<typename StateType, typename VolType, typename RandomType>
pair<StateType, StateType> Euler<StateType, VolType, RandomType>::pairTransition(mt19937_64 &gen, double t, pair<StateType, StateType> pairProcess, double h, unsigned int n)
{
    // We get the finest step
    double h_n = h/n;

    // We get the starting points
    StateType X_0 = pairProcess.first;
    StateType X_n = pairProcess.second;

    RandomType randomRealization;
    RandomType randomRealizationSum;

    // We make n transitions on the finest process, keeping memory of the random realizations
    for (unsigned int i=0; i<n; ++i){
        randomRealization = this->m_model->random(gen)*sqrt(h_n);
        X_n = transition(t + i*h_n, X_n, h_n, randomRealization);
        randomRealizationSum += randomRealization;
    }

    // We make one transition on the biggest process
    X_0 = transition(t, X_0, h, randomRealizationSum);

    return pair<StateType, StateType>(X_0, X_n);
}


/**
 * @brief Method to simulate a realization of \f$X_T\f$.
 *
 * @param gen Generator.
 * @param n Discretization over time.
 * @return Simulated variable \f$X_T\f$.
 */
template<typename StateType, typename VolType, typename RandomType>
StateType Euler<StateType, VolType, RandomType>::singleSimulation(mt19937_64& gen, const unsigned int n)
{
    double currentTime = 0.;

    StateType X_i = this->m_model->getStartingPoint();
    double h = this->m_model->getMaturity()/(double)n;
    for (unsigned int i=0; i<n; i++){
        RandomType random = sqrt(h)*this->m_model->random(gen);
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
template<typename StateType, typename VolType, typename RandomType>
pair<StateType, StateType> Euler<StateType, VolType, RandomType>::doubleSimulation(mt19937_64& gen, const unsigned int n1, const unsigned int n2)
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

    StateType X_n1= this->m_model->getStartingPoint();
    StateType X_n2 = this->m_model->getStartingPoint();

    double h1 = this->m_model->getMaturity()/(double)n1;

    pair<StateType, StateType> pairProcess(X_n1, X_n2);

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




//class PhiScheme : public Scheme
//{
//public:
//    /** Constructor. */
//    PhiScheme(const modelPtr model, const schemePtr scheme, std::function<double(double)> f);

//    double transition(double t, double x, double h, double random);
//    pair<double,double> pairTransition(mt19937_64& gen, double t, pair<double, double> x, double h, unsigned int n);

//    double singleSimulation(mt19937_64& gen, const unsigned int n);
//    pair<double, double> doubleSimulation(mt19937_64& gen, const unsigned int n1, const unsigned int n2);

//protected:
//    schemePtr m_scheme;
//    std::function<double(double)> m_f;
//};

//// Smart pointers to PhiScheme objects
//typedef shared_ptr<PhiScheme> phiSchemePtr;

//template<typename state> void saySomething(state thisstate)
//{
//    cout << thisstate << endl;
//}

// Smart pointer to Scheme objects
template<typename StateType, typename VolType, typename RandomType>
using schemePtr = shared_ptr<Scheme<StateType, VolType, RandomType>> ;

// Smart pointer to Euler objects
template<typename StateType, typename VolType, typename RandomType>
using eulerPtr = shared_ptr<Euler<StateType, VolType, RandomType>>;

#endif // SCHEME_HPP
