#ifndef SCHEME_HPP
#define SCHEME_HPP

#include "model.hpp"

/**
 * @brief The Scheme class
 *
 * Generic class to define a scheme on the model that will be used to simulate the process.
 *
 */
template<typename StateType, typename VolType, typename RandomType, typename TransitionType>
class Scheme
{
public:
    /** Constructor. */
    Scheme(const modelPtr<StateType, VolType, RandomType> model);
    /** Destructor. */
    virtual ~Scheme(){}

    /** Pure virtual function to make a transion. */
    virtual TransitionType transition(double t, TransitionType x, double h, RandomType random) = 0;

    /** Pure virtual function to make a coupled transion. */
    virtual pair<TransitionType,TransitionType> pairTransition(mt19937_64& gen, double t, pair<TransitionType, TransitionType> x, double h, unsigned int n) = 0;

    /** Pure virtual function to simulate a process. */
    virtual TransitionType singleSimulation(mt19937_64& gen, const unsigned int n) = 0;

    /** Pure virtual function to simulate two processes. */
    virtual pair<TransitionType, TransitionType> doubleSimulation(mt19937_64& gen, const unsigned int n1, const unsigned int n2) = 0;

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
class Euler : public Scheme<StateType, VolType, RandomType, StateType>
{
public:
    /** Constructor. */
    Euler(const modelPtr<StateType, VolType, RandomType> model);
    /** Destructor. */
    ~Euler(){}

    StateType transition(double t, StateType x, double h, RandomType random);
    pair<StateType,StateType> pairTransition(mt19937_64& gen, double t, pair<StateType, StateType> x, double h, unsigned int n);

    StateType singleSimulation(mt19937_64& gen, const unsigned int n);
    pair<StateType, StateType> doubleSimulation(mt19937_64& gen, const unsigned int n1, const unsigned int n2);
};



/**
 * @param model Pointer to the model.
 */
template<typename StateType, typename VolType, typename RandomType, typename TransitionType>
Scheme<StateType, VolType, RandomType, TransitionType>::Scheme(const modelPtr<StateType, VolType, RandomType> model): m_model(model){}


/**
 * @param model Pointer to the model.
 */
template<typename StateType, typename VolType, typename RandomType>
Euler<StateType, VolType, RandomType>::Euler(const modelPtr<StateType, VolType, RandomType> model): Scheme<StateType, VolType, RandomType, StateType>(model){}


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
    StateType drift = this->m_model->drift(t,x);
    VolType sigma = this->m_model->sigma(t, x);
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

// Smart pointer to Scheme objects
template<typename StateType, typename VolType, typename RandomType, typename TransitionType>
using schemePtr = shared_ptr<Scheme<StateType, VolType, RandomType, TransitionType>> ;

// Smart pointer to Euler objects
template<typename StateType, typename VolType, typename RandomType>
using eulerPtr = shared_ptr<Euler<StateType, VolType, RandomType>>;



/*------------------------------------------------------------------------------*/
/*                              PHI SCHEME                                      */
/*------------------------------------------------------------------------------*/

template <typename StateType, typename VolType, typename RandomType>
class PhiScheme : public Scheme<StateType, VolType, RandomType, pair<StateType,StateType>>
{
public:
    /** Constructor. */
    PhiScheme(const modelPtr<StateType, VolType, RandomType> model,
              const schemePtr<StateType, VolType, RandomType, StateType> scheme,
              std::function<StateType(StateType, StateType)> phi);
    /** Destructor. */
    ~PhiScheme(){}

    pair<StateType,StateType>  transition(double t, pair<StateType,StateType> x, double h, RandomType random);
    pair<pair<StateType,StateType>,pair<StateType,StateType>> pairTransition(mt19937_64& gen, double t, pair<pair<StateType,StateType>,pair<StateType,StateType>>  x, double h, unsigned int n);

    pair<StateType,StateType> singleSimulation(mt19937_64& gen, const unsigned int n);
    pair<pair<StateType,StateType>,pair<StateType,StateType>> doubleSimulation(mt19937_64& gen, const unsigned int n1, const unsigned int n2);

protected:
    schemePtr<StateType, VolType, RandomType, StateType> m_scheme;
    std::function<StateType(StateType, StateType)> m_phi;
};

/**
 * @param model Pointer to the model
 * @param scheme Pointer to the scheme
 */
template <typename StateType, typename VolType, typename RandomType>
PhiScheme<StateType, VolType, RandomType>::PhiScheme(const modelPtr<StateType, VolType, RandomType> model,
                                                     const schemePtr<StateType, VolType, RandomType, StateType> scheme,
                                                     std::function<StateType(StateType, StateType)> phi):
    Scheme<StateType, VolType, RandomType, pair<StateType,StateType>>(model), m_scheme(scheme), m_phi(phi){}


/**
 * @brief Method to make a transition of step h, at time t, starting from x and with a random realization random, applying at each step the functional \f$\phi\f$.
 *
 * @param t Starting time.
 * @param x Starting point.
 * @param h Time step.
 * @param random Random realization.
 * @return
 */
template <typename StateType, typename VolType, typename RandomType>
pair<StateType,StateType>  PhiScheme<StateType, VolType, RandomType>::transition(double t, pair<StateType,StateType> x, double h, RandomType random)
{
    StateType S_t = x.first;
    StateType phi_S_t = x.second;

    S_t = m_scheme->transition(t, S_t, h, random);
    phi_S_t = m_phi(S_t, phi_S_t);

    return pair<StateType,StateType>(S_t, phi_S_t);
}

/**
 * @brief Method to make a transition for two different discretizations of the same process, over step h, one with discretization 1 and one with discretization n,
 * starting from x at time t, applying at each step the functional \f$\phi\f$.
 *
 * @param gen Random generator.
 * @param t Starting time.
 * @param x Starting point.
 * @param h Time step.
 * @param n Discretization for the finest process.
 */
template<typename StateType, typename VolType, typename RandomType>
pair<pair<StateType,StateType>, pair<StateType,StateType>> PhiScheme<StateType, VolType, RandomType>::pairTransition(mt19937_64 &gen, double t,
                                                                                                                     pair<pair<StateType,StateType>, pair<StateType,StateType>> pairProcess,
                                                                                                                     double h,
                                                                                                                     unsigned int n)
{
    // We get the finest step
    double h_n = h/n;

    // We get the starting points
    pair<StateType,StateType> X_phiX_0 = pairProcess.first;
    pair<StateType,StateType> X_phiX_n = pairProcess.second;

    RandomType randomRealization;
    RandomType randomRealizationSum;

    // We make n transitions on the finest process, keeping memory of the random realizations
    for (unsigned int i=0; i<n; ++i){
        randomRealization = this->m_model->random(gen)*sqrt(h_n);
        X_phiX_n = transition(t + i*h_n, X_phiX_n, h_n, randomRealization);
        randomRealizationSum += randomRealization;
    }

    // We make one transition on the biggest process
    X_phiX_0 = transition(t, X_phiX_0, h, randomRealizationSum);

    return pair<pair<StateType,StateType>, pair<StateType,StateType>>(X_phiX_0, X_phiX_n);
}

/**
 * @brief Method to simulate a realization of \f$X_T\f$.
 *
 * @param gen Generator.
 * @param n Discretization over time.
 * @return Simulated variable \f$X_T\f$.
 */
template<typename StateType, typename VolType, typename RandomType>
pair<StateType, StateType> PhiScheme<StateType, VolType, RandomType>::singleSimulation(mt19937_64& gen, const unsigned int n)
{
    double currentTime = 0.;

    StateType X_i = this->m_model->getStartingPoint();
    StateType phi_X_i = m_phi(X_i, X_i);

    pair<StateType,StateType> pair_i(X_i,phi_X_i);

    double h = this->m_model->getMaturity()/(double)n;
    for (unsigned int i=0; i<n; i++){
        RandomType random = sqrt(h)*this->m_model->random(gen);
        pair_i = transition(currentTime, pair_i, h, random);
        currentTime += h;
    }

    return pair_i;
}

/**
 * @brief Method to simulate two realizations of \f$(X_T, \phi(X_0, \ldots, X_T))\f$.
 *
 * The two discretizations used are \f$\frac{T-t_0}{n_1}\f$ and \f$\frac{T-t_0}{n_2}\f$.
 *
 * @param gen Generator.
 * @param n1 Discretization over time for first simulation.
 * @param n2 Discretization over time for second simulation.
 * @return Pair of simulations of \f$(X_T, \phi(X_0, \ldots, X_T))\f$.
 */
template<typename StateType, typename VolType, typename RandomType>
pair<pair<StateType,StateType>, pair<StateType,StateType>> PhiScheme<StateType, VolType, RandomType>::doubleSimulation(mt19937_64& gen, const unsigned int n1, const unsigned int n2)
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

    StateType phi_X_n1= m_phi(X_n1,X_n1);
    StateType phi_X_n2 = m_phi(X_n2,X_n2);

    pair<StateType,StateType> pair_n1(X_n1, phi_X_n1);
    pair<StateType,StateType> pair_n2(X_n2, phi_X_n2);

    double h1 = this->m_model->getMaturity()/(double)n1;

    pair<pair<StateType,StateType>, pair<StateType,StateType>> pairProcess(pair_n1, pair_n2);

    // We make a loop over the largest discretization and we make a pair transition over each step
    for (unsigned int i=0; i<n1; ++i){
        pairProcess = pairTransition(gen, currentTime, pairProcess, h1, div);
        currentTime += h1;
    }

    return pairProcess;
}

// Smart pointer to Scheme objects
template<typename StateType, typename VolType, typename RandomType>
using phiSchemePtr = shared_ptr<PhiScheme<StateType, VolType, RandomType>> ;

#endif // SCHEME_HPP
