#ifndef SCHEME_HPP
#define SCHEME_HPP

#include "model.hpp"

class Scheme
{
public:
    /** Constructor. */
    Scheme(const modelPtr model);

    /** Pure virtual function to make a transion. */
    virtual Eigen::VectorXd transition(double t, Eigen::VectorXd x, double h, Eigen::VectorXd random) = 0;

    /** Pure virtual function to make a coupled transion. */
    virtual pair<Eigen::VectorXd,Eigen::VectorXd> pairTransition(mt19937_64& gen, double t, pair<Eigen::VectorXd, Eigen::VectorXd> x, double h, unsigned int n) = 0;

    /** Pure virtual function to simulate a process. */
    virtual Eigen::VectorXd singleSimulation(mt19937_64& gen, const unsigned int n) = 0;

    /** Pure virtual function to simulate two processes. */
    virtual pair<Eigen::VectorXd, Eigen::VectorXd> doubleSimulation(mt19937_64& gen, const unsigned int n1, const unsigned int n2) = 0;

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

    Eigen::VectorXd transition(double t, Eigen::VectorXd x, double h, Eigen::VectorXd random);
    pair<Eigen::VectorXd,Eigen::VectorXd> pairTransition(mt19937_64& gen, double t, pair<Eigen::VectorXd, Eigen::VectorXd> x, double h, unsigned int n);

    Eigen::VectorXd singleSimulation(mt19937_64& gen, const unsigned int n);
    pair<Eigen::VectorXd, Eigen::VectorXd> doubleSimulation(mt19937_64& gen, const unsigned int n1, const unsigned int n2);
};

// Smart pointers to Euler objects
typedef shared_ptr<Euler> eulerPtr;

//class PhiScheme
//{
//public:
//    /** Constructor. */
//    PhiScheme(const modelPtr model, const schemePtr scheme, std::function<Eigen::VectorXd(Eigen::VectorXd)> phi);

//    pair<Eigen::VectorXd,Eigen::VectorXd>  transition(double t, pair<Eigen::VectorXd,Eigen::VectorXd> x, double h, Eigen::VectorXd random);
//    pair<pair<Eigen::VectorXd,Eigen::VectorXd>,pair<Eigen::VectorXd,Eigen::VectorXd>> pairTransition(mt19937_64& gen, double t, pair<pair<Eigen::VectorXd,Eigen::VectorXd>,pair<Eigen::VectorXd,Eigen::VectorXd>>  x, double h, unsigned int n);

//    pair<Eigen::VectorXd, Eigen::VectorXd> singleSimulation(mt19937_64& gen, const unsigned int n);
//    pair<pair<Eigen::VectorXd,Eigen::VectorXd>,pair<Eigen::VectorXd,Eigen::VectorXd>> doubleSimulation(mt19937_64& gen, const unsigned int n1, const unsigned int n2);

//protected:
//    modelPtr m_model;
//    schemePtr m_scheme;
//    std::function<double(double)> m_phi;
//};

// Smart pointers to PhiScheme objects
typedef shared_ptr<PhiScheme> phiSchemePtr;

//template<typename state> void saySomething(state thisstate)
//{
//    cout << thisstate << endl;
//}



#endif // SCHEME_HPP
