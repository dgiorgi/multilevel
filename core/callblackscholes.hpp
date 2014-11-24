#include <iostream>
#include <fstream>
#include <random>
#include <chrono>
#include <vector>
#include <iomanip>  // for setprecision

#include "functions.hpp"
#include "model.hpp"
#include "montecarlo.hpp"
#include "structuralparameters.hpp"
#include "multilevelparameters.hpp"
#include "estimator.hpp"
#include "scheme.hpp"

std::vector<double> callBlackScholes(double x0=100, double r=0.05, double sigma=0.2, double K=100, double T=5, unsigned int seed = 0);
