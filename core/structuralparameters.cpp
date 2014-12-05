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


    file_out << "# STRUCTURAL PARAMETERS : " << endl;

    file_out << endl;

    file_out << "# Alpha : "      << m_alpha  << endl;
    file_out << "# Beta : "       << m_beta   << endl;
    file_out << "# c_1 : "        << m_c1     << endl;
    file_out << "# h bold : "     << m_hBold  << endl;
    file_out << "# V_1 : "        << m_V1     << endl;
    file_out << "# var(Y_0) : "   << m_varY0  << endl;
    file_out << "# theta: "       << m_theta  << endl;

    file_out << endl;
    file_out << "# -------------------------------" << endl << endl;

    file_out.close();
}
