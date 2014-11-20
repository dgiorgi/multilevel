#ifndef SCHEME_HPP
#define SCHEME_HPP

#include "model.hpp"

class Scheme
{
public:
    /** Constructor. */
    Scheme(const modelPtr model);

    /** Virtual function to simulate a process. */
    virtual double singleSimulation(mt19937_64& gen, const unsigned int n) = 0;

    /** Virtual function to simulate two processes. */
    virtual pair<double, double> doubleSimulation(mt19937_64& gen, const unsigned int n1, const unsigned int n2) = 0;

    modelPtr getModel() const {return m_model;}

protected:
    modelPtr m_model;
};

// Smart pointers to Scheme objects
typedef shared_ptr<Scheme> schemePtr;

class Euler : public Scheme
{
public:
    /** Constructor. */
    Euler(const modelPtr model);

    double singleSimulation(mt19937_64& gen, const unsigned int n);
    pair<double, double> doubleSimulation(mt19937_64& gen, const unsigned int n1, const unsigned int n2);
};

// Smart pointers to Euler objects
typedef shared_ptr<Euler> eulerPtr;

class PhiScheme : public Scheme
{
public:
    /** Constructor. */
    PhiScheme(const modelPtr model, const schemePtr scheme, std::function<double(double)> f);
    double singleSimulation(mt19937_64& gen, const unsigned int n);
    pair<double, double> doubleSimulation(mt19937_64& gen, const unsigned int n1, const unsigned int n2);

protected:
    schemePtr m_scheme;
    std::function<double(double)> m_f;
};



//template<typename state> void saySomething(state thisstate)
//{
//    cout << thisstate << endl;
//}



#endif // SCHEME_HPP
