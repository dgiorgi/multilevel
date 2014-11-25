/* File: callblackscholes.i */
%module callblackscholes

/* This is mandatory for the vector<type> use in Python */  
%include "std_vector.i"

// Instantiate templates 
namespace std {
   %template(IntVector) vector<int>;
   %template(DoubleVector) vector<double>;
}

%header %{
#include "callblackscholes.hpp"
%}

// Include the header file with above prototypes
%include "callblackscholes.hpp"
