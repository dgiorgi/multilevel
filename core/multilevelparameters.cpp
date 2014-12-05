#include "multilevelparameters.hpp"

using namespace std;


void Refiners::clear()
{
    m_data.clear();
    m_M = 2;
}

/**
 * @brief Constructor.
 *
 * @param epsilon Precision.
 * @param structParam Structural parameters.
 * @param type Estimator type, can be MC or RR.
 * @param M Root.
 * @param R Order.
 * @param h Biais.
 * @param N Estimator size.
 */
MultilevelParameters::MultilevelParameters(double epsilon,
                                           StructuralParameters& structParam,
                                           estimator_type type,
                                           int M,
                                           int R,
                                           double h,
                                           double N):
    m_epsilon(epsilon), m_structParam(structParam), m_type(type), m_M(M),  m_R(R), m_h(h), m_N(N), m_flagDoneComputation(false)
{
    if (M == -1)
        m_flagGivenRoot = false;
    else
        m_flagGivenRoot = true;
    if (R == -1)
        m_flagGivenOrder = false;
    else
        m_flagGivenOrder = true;
    if (fabs(h+1) < 1e-5)
        m_flagGivenBiais = false;
    else
        m_flagGivenBiais = true;
    if (fabs(N+1) < 1e-5)
        m_flagGivenEstimatorSize = false;
    else
        m_flagGivenEstimatorSize = true;

    computeOptimalParameters();
}

void MultilevelParameters::initialize()
{
    m_flagDoneComputation = false;

    m_W.clear();
    m_n.clear();
    m_q.clear();
}

/**
 * @brief Generic method to compute all the optimal parameters.
 *
 * This computes the optimal parameters \f$R, h, q, N\f$, choosing the \f$M\f$ which
 * optimizes the cost.
 * It computes also the refiners and (if the estimator type is RR)
 * the weights of the allocation matrix.
 */
void MultilevelParameters::computeOptimalParameters()
{
    if (!m_flagGivenRoot){
        double bestCost = MY_INF;
        int bestM = 2;

        for (int M=2; M<10; ++M){

            computeOptimalParametersForM(M);

            double tempCost = computeCost();

            bestM = tempCost < bestCost ? M : bestM;
            bestCost = min(tempCost, bestCost);
        }

        m_M = bestM;
        m_cost = bestCost;
        computeOptimalParametersForM(bestM);
    }
    else{
        computeOptimalParametersForM(m_M);
        computeCost();
    }
}

/**
 * @brief Generic method to compute all the optimal parameters.
 *
 * This computes the optimal parameters \f$R, h, q, N\f$ for a given \f$M\f$.
 * It computes also the refiners and (if the estimator type is RR)
 * the weights of the allocation matrix.
 */
void MultilevelParameters::computeOptimalParametersForM(const int M)
{
    m_M = M;

    initialize();

    // First we compute the order R
    if (!m_flagGivenOrder)
        computeOrder();

    // Once we have R, we can build the refiners
    m_n = Refiners(m_R,m_M);

    // We compute the biais h
    if (!m_flagGivenBiais)
        computeBiais();
    else
        m_hInverse = 1./m_h;

    // If it's a Richardson-Romberg estimator we compute the weights of the allocation matrix
    if (m_type == RR){
        computeWeights();
    }

    // We compute the stratification
    computeStratification();

    // We compute the estimator size
    if (!m_flagGivenEstimatorSize)
        computeN();

    m_flagDoneComputation = true;
}


/**
 * @brief Method to compute the order \f$R\f$.
 */
void MultilevelParameters::computeOrder()
{
    double alpha = m_structParam.getAlpha();

    switch (m_type){
    case RR:{
        double c_tilde = 1.;
        double H = m_structParam.getHBold();

        // We round to the nearest integer floor(0.5+R)
        double R = 0.5 + log(pow(c_tilde, 1/alpha)*H)/log((double)m_M)
                + sqrt(pow(0.5+log(pow(c_tilde, 1./alpha)*H)/log((double)m_M),2) + 2.*log(1./m_epsilon)/(alpha*log((double)m_M)));
        m_R = ceil(R); //floor(0.5 + R);
        break;
    }
    case MC:{
        double c1 = m_structParam.getC1();
        double R = 1. + log(fabs(c1)/m_epsilon)/(alpha*log((double)m_M));

        m_R = ceil(R); //floor(0.5 + R);
        break;
    }
    }

    if (m_R<2)
        m_R = 2;
}

/**
 * @brief Method to compute the biais \f$h\f$.
 */
void MultilevelParameters::computeBiais()
{
    double alpha = m_structParam.getAlpha();

    switch(m_type){
    case RR:
    {
        double f1 = pow( 1.0+2.*alpha*(double)m_R , 1.0/(2.*alpha*(double)m_R) );
        double f2 = pow( m_epsilon , -1.0/(alpha*(double)m_R) );
        double f3 = pow( (double)m_M, -((double)m_R-1.)*0.5 );

        m_hInverse = ceil(f1*f2*f3);
        break;
    }
    case MC:
    {
        double c1 = m_structParam.getC1();

        double f1 = pow( 1.+2.*alpha , 1.0/(2.*alpha) );
        double f2 = pow( m_epsilon/fabs(c1) , -1.0/(alpha) );
        double f3 = pow( (double)m_M, -((double)m_R-1) );

        m_hInverse = ceil(f1*f2*f3);
        break;
    }
    }

    m_h = 1/(double)m_hInverse;
}

/**
 * @brief Method to compute the allocation matrix weights \f$W_i\f$.
 *
 * This method is called only for the RR multilevel estimator.
 * \f$W_1 = 1, \; W_j = \sum_{k=j}^R w_k\f$ for \f$j=2,\ldots,R\f$.
 */
void MultilevelParameters::computeWeights()
{
    double alpha = m_structParam.getAlpha();

    // First we compute the weight vector as the unique solution of the Vandermonde system
    vector<double> w = vector<double>();

    for (int i=0; i<m_R; ++i){
        double num = pow(-1, (double)m_R-(double)(i+1)) * pow((double)m_n[i], alpha*((double)m_R-1.));
        double prodInf = 1;
        for(int j=0; j<i; ++j)
            prodInf *= pow((double)m_n[i], alpha) - pow((double)m_n[j], alpha);
        double prodSup = 1;
        for(int j=i+1; j<m_R; ++j)
            prodSup *= pow((double)m_n[j], alpha) - pow((double)m_n[i], alpha);

        w.push_back(num / (prodInf*prodSup));
    }

    // Then we compute the weights of the allocation matrix as partial sums of the weight
    // vector elements.
    m_W.push_back(1.0);
    for (int j=1; j<m_R; ++j){
        double sum = 0.0;
        for (int k=j; k<m_R; ++k){
            sum += w[k];
        }
        m_W.push_back(sum);
    }
}

/**
 * @brief Method to compute the stratification \f$q=(q_1, \ldots, q_R)\f$.
 *
 * \f$q_i>0\f$ and \f$\sum_jq_j=1\f$.
 * \f$N_j = q_j *N\f$ for \f$j=1,\ldots,R\f$.
 */
void MultilevelParameters::computeStratification()
{
    double theta = m_structParam.getTheta();
    double beta = m_structParam.getBeta();

    double num = 0.;
    double denum = 0.;

    // First we compute the unnormalized elements
    vector<double> temp_q = vector<double>();

    temp_q.push_back(1 + theta*pow(m_h, 0.5*beta));

    for (int j=1; j<m_R; ++j){
        num = pow((double)m_n[j-1], -0.5*beta) + pow((double)m_n[j], -0.5*beta);
        denum = sqrt(m_n[j-1] + m_n[j]);

        temp_q.push_back(theta*pow(m_h, 0.5*beta)*num/denum);

        if (m_type == RR)
            temp_q[j] *= fabs(m_W[j]);
    }

    // Then we normalize
    double sum = 0;
    for (int j=0; j<m_R; ++j)
        sum += temp_q[j];

    for (int j=0; j<m_R; ++j)
        m_q.push_back(temp_q[j]/sum);
}

/**
 * @brief Method to compute the estimator size \f$N\f$.
 */
void MultilevelParameters::computeN()
{
    double alpha = m_structParam.getAlpha();
    double beta = m_structParam.getBeta();
    double theta = m_structParam.getTheta();
    double varY0 = m_structParam.getVarY0();

    double sum = 0.;

    // We compute the first element of the sum
    // We remember n_0 = n_0^-1 = 0
    double value = pow((double)m_n[0], -0.5*beta) * sqrt((double)m_n[0]);
    if(m_type == RR)
        value *= m_W[0];
    sum += value;

    for (int j=1; j<m_R; ++j){
        value = (pow((double)m_n[j-1], -0.5*beta) + pow((double)m_n[j],-0.5*beta)) * sqrt((double)m_n[j-1]+(double)m_n[j]);
        if (m_type == RR)
            value *= m_W[j];
        sum += value;
    }

    double num = varY0*pow(1+theta*pow(m_h, 0.5*beta)*sum,2);

    sum = 0;
    sum += m_q[0]*m_n[0];
    for (int j=1; j<m_R; ++j)
        sum += m_q[j]*(m_n[j-1]+m_n[j]);

    double denum = pow(m_epsilon,2)*sum;

    double coeff = 0;
    switch (m_type){
    case MC:
        coeff = (1. + 0.5/alpha);
        break;
    case RR:
        coeff = (1. + 0.5/(alpha*(double)m_R));
        break;
    }

    m_N = coeff*num/denum;
}

/**
 *
 * @return Computational cost of the estimator.
 */
double MultilevelParameters::computeCost()
{
    if (!m_flagDoneComputation){
        cerr << "Cannot compute the computation cost before computing the multilevel parameters." << endl;
        cerr << "Please be sure you computed the multilevel parameters first." << endl;
        exit(1);
    }

    double sum = m_q[0]*m_n[0];
    for (int j=1; j<m_R; ++j)
        sum += m_q[j]*(m_n[j-1]+m_n[j]);

    sum /= m_h;

    return m_N*sum;
}


/**
 * @brief Method to display the parameters.
 */
void MultilevelParameters::displayParameters()
{
    cout << "Epsilon : " << m_epsilon << endl << endl;
    if (m_type == 0)
        cout << "MLMC MULTILEVEL PARAMETERS : " << endl;
    else
        cout << "ML2R MULTILEVEL PARAMETERS : " << endl;
    cout << "Order R : "        << m_R << endl;
    cout << "Root M : "         << m_M  << endl;
    cout << "Biais inverse h^-1 : " << m_hInverse << endl;
    cout << "Number of simulations N : " << m_N << endl;
    cout << "Estimator cost : " << m_cost << endl;
    cout << "Stratification strategy q : (";
    for (int i=0; i<m_q.size()-1; ++i)
        cout << m_q[i] << ", ";
    cout << m_q[m_q.size()-1] << ")" << endl;

    cout << endl;
}

/**
 * @brief Method to write the parameters in a file.
 */
void MultilevelParameters::writeParameters(const string fileName)
{
    // Open the file where we want to write the results
    ofstream file_out(fileName.c_str(), ios::out | ios::binary | ios::app);

    file_out.precision(6);

    if (!file_out)  { // if the opening succeeded
        cerr << endl << " Write parameters in file error! " << endl;
        cerr << " Cannot open solution file (" << fileName << ") to write the multilevel parameters." << endl << endl;
        exit(1);
    }

    file_out << "# Epsilon : "        << m_epsilon << endl << endl;

    file_out << "# MULTILEVEL PARAMETERS :" << endl << endl;

    file_out << "# Order R : "        << m_R << endl;
    file_out << "# Root M : "         << m_M  << endl;
    file_out << "# Biais inverse h^-1 : " << m_hInverse << endl;
    file_out << "# Stratification strategy q : (";
    for (int i=0; i<m_q.size()-1; ++i)
        file_out << m_q[i] << ", ";
    file_out << m_q[m_q.size()-1] << ")" << endl;

    file_out << "# Number of simulations N : " << m_N << endl;
    file_out << "# Estimator cost : " << m_cost << endl;
    file_out << endl;

    file_out.close();
}
