#include <iostream>
#include <cmath>
#include <fstream>

using std::ostream;
using std::ofstream;
using std::ifstream;
using std::cout;
using std::endl;
using std::ios;

#include "include.h"

int xSPM::M;
int xSPM::N;

/**
 * initialize the static variables and allocate and initialize the static lists
 * @param M_in nr of spatial sp orbs
 * @param N_in nr of particles
 */
void xSPM::init(int M_in,int N_in){

   M = M_in;
   N = N_in;

}

/**
 * standard constructor
 */
xSPM::xSPM() : BlockMatrix(2) {

   //set the dimension and degeneracy of the two blocks:
   this->setMatrixDim(0,M,1);
   this->setMatrixDim(1,M,3);

}

/**
 * copy constructor: constructs Matrix object and fills it with the content of matrix xspm_c
 * @param xspm_c object that will be copied into this.
 */
xSPM::xSPM(const xSPM &xspm_c) : BlockMatrix(xspm_c){ }

/**
 * destructor
 */
xSPM::~xSPM(){ }

/**
 * @return number of particles
 */
int xSPM::gN() const{

   return N;

}

/**
 * @return number of sp orbitals
 */
int xSPM::gM() const{

   return M;

}
