#ifndef SCHEME_HPP
#define SCHEME_HPP

#include <model.hpp>

class Scheme
{
public:
    /** Constructor. */
    Scheme(const Model& model);

    /** Virtual function to simulate a process. */
    virtual double singleSimulation(mt19937_64& gen, const unsigned int n) = 0;

    /** Virtual function to simulate two processes. */
    virtual pair<double, double> doubleSimulation(mt19937_64& gen, const unsigned int n1, const unsigned int n2) = 0;

protected:
    Model m_model;
};

class Euler : public Scheme
{
public:
    /** Constructor. */
    Euler(const Model& model);

    double singleSimulation(mt19937_64& gen, const unsigned int n);
    pair<double, double> doubleSimulation(mt19937_64& gen, const unsigned int n1, const unsigned int n2);

protected:

};


#endif // SCHEME_HPP
