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

int rSPM::M;
int rSPM::N;

/**
 * initialize the static variables and allocate and initialize the static lists
 * @param M_in nr of spatial sp orbs
 * @param N_in nr of particles
 */
void rSPM::init(int M_in,int N_in){

   M = M_in;
   N = N_in;

}

/**
 * standard constructor
 */
rSPM::rSPM() : BlockMatrix(2) {

   //set the dimension and degeneracy of the two blocks:
   this->setMatrixDim(0,M,1);
   this->setMatrixDim(1,M,3);

}

/**
 * copy constructor: constructs Matrix object and fills it with the content of matrix rspm_c
 * @param rspm_c object that will be copied into this.
 */
rSPM::rSPM(const rSPM &rspm_c) : BlockMatrix(rspm_c){ }

/**
 * destructor
 */
rSPM::~rSPM(){ }

/**
 * @return number of particles
 */
int rSPM::gN() const{

   return N;

}

/**
 * @return number of sp orbitals
 */
int rSPM::gM() const{

   return M;

}
