#include "structuralparameters.hpp"

using namespace std;

/**
 * @brief Constructor.
 *
 * @param alpha \f$\alpha\f$
 * @param beta \f$\beta\f$
 * @param c1 \f$c_1\f$
 * @param h \f$h\f$
 */
StructuralParameters::StructuralParameters(const double alpha,
                                           const double beta,
                                           const double c1,
                                           const double h):
    m_alpha(alpha), m_beta(beta), m_c1(c1), m_hBold(h)
{}

/**
 * @brief Generic method to compute all the structural parameters
 * \f$V_1, var(Y_0)\f$ and \f$\theta\f$.
 *
 * @param gen Generator for the random variable.
 * @param f Function \f$f\f$ of the random variable.
 * @param model Model for the random variable.
 * @param N Number of simulations.
 */
void StructuralParameters::computeParameters(mt19937_64& gen,
                                             std::function<double(double)> f,
                                             const modelfPtr model,
                                             const unsigned int N)
{
    // First we compute the var(Y0) and V1
    computeVarY0(gen, f, model, N);
    computeV1(gen, f, model, N);

    // And then theta
    computeTheta();

    // We display the parameters
    displayParameters();
    writeParameters("parameters.txt");
}


/**
 * @brief Method to compute the parameter \f$V_1\f$. We consider the following estimator:
 * \f[V_1(h) = (1 + M_{max}^{\frac{\beta}{2}})^{-2}h^{-\beta}\|Y_h-Y_{h/M_{max}}\|_2^2\f]
 * where we set \f$M_{max} = 10\f$ and \f$h = \mathbf{h}.\f$.
 *
 * @param gen Generator for the random variable.
 * @param f Function \f$f\f$ of the random variable.
 * @param model Model for the random variable.
 * @param N Number of simulations.
 */
void StructuralParameters::computeV1(mt19937_64& gen,
                                     std::function<double(double)> f,
                                     const modelfPtr model,
                                     const unsigned int N)
{
    unsigned int M = 10;
    double T = model->getMaturity();

    // Throw error if the rest of T/m_hBold is not 0
    double epsilon = 1e-5;
    if (fabs(T/m_hBold-(int)(T/m_hBold)) > epsilon){
        cerr << "The discretization step is not an integer divisor of the maturity. Please choose a good discretization step." << endl;
        exit(1);
    }

    int n = T/m_hBold;
    double beta = m_beta;

    // We instanciate a double monte carlo
    DoubleMonteCarlo Y = DoubleMonteCarlo(gen, f, model, n, M*n);
    // We make N simulations
    Y(N);

    double L2 = Y.meanOfSquares();

    m_V1 = pow(1.0+pow((double)M, -beta*0.5), -2) * pow(m_hBold,-beta) * L2;
}

/**
 * @brief Method to compute the variance \f$var(Y_0) \sim var(Y_{\mathbf{h}})\f$.
 *
 * We take here the easiest choice h=min(1,T) i.e. 1 discretization.
 *
 * @param gen Generator for the random variable.
 * @param f Function \f$f\f$ of the random variable.
 * @param model Model for the random variable.
 * @param N Number of simulations.
 */
void StructuralParameters::computeVarY0(mt19937_64& gen,
                                        std::function<double(double)> f,
                                        const modelfPtr model,
                                        const unsigned int N)
{
    double T = model->getMaturity();

    // Throw error if the rest of T/m_hBold is not 0
    double epsilon = 1e-5;
    if (fabs(T/m_hBold-(int)(T/m_hBold)) > epsilon){
        cerr << "The discretization step is not an integer divisor of the maturity. Please choose a good discretization step." << endl;
        exit(1);
    }

    int n = T/m_hBold;

    // We instanciate a single Monte Carlo
    MonteCarlo Y = MonteCarlo(gen, f, model, n);
    // We make N simulations
    Y(N);
    m_varY0 = Y.var();
}

/**
 * @brief Method to compute \f$\theta = \sqrt{\frac{V_1}{var(Y_0)}}\f$.
 * */
void StructuralParameters::computeTheta()
{
    m_theta = sqrt(m_V1/m_varY0);
}

/**
 * @brief Method to display the parameters.
 */
void StructuralParameters::displayParameters()
{
    cout << "Alpha : "      << m_alpha  << endl;
    cout << "Beta : "       << m_beta   << endl;
    cout << "c_1 : "        << m_c1     << endl;
    cout << "h bold : "     << m_hBold  << endl;
    cout << "V_1 : "        << m_V1     << endl;
    cout << "var(Y_0) : "   << m_varY0  << endl;
    cout << "theta: "       << m_theta  << endl;

    cout << endl;
}

/**
 * @brief Method to write the parameters in a file.
 */
void StructuralParameters::writeParameters(const string fileName)
{
    // Open the file where we want to write the results
    ofstream file_out(fileName.c_str(), ios::out | ios::binary);

    file_out.precision(6);

    if (!file_out)  { // if the opening succeeded
        cerr << endl << " Write parameters in file ERROR! " << endl;
        cerr << " Cannot open solution file (" << fileName << ") to write the structural parameters." << endl << endl;
        exit(1);
    }


    file_out << "Structural parameters for multilevel estimator" << endl;

    file_out << endl;

    file_out << "Alpha : "      << m_alpha  << endl;
    file_out << "Beta : "       << m_beta   << endl;
    file_out << "c_1 : "        << m_c1     << endl;
    file_out << "h bold : "     << m_hBold  << endl;
    file_out << "V_1 : "        << m_V1     << endl;
    file_out << "var(Y_0) : "   << m_varY0  << endl;
    file_out << "theta: "       << m_theta  << endl;

    file_out << endl;

    file_out.close();
}
