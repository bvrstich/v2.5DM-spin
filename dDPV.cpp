#include <iostream>
#include <cmath>
#include <fstream>

using std::ostream;
using std::ofstream;
using std::ifstream;
using std::endl;
using std::ios;

#include "include.h"

int dDPV::M;
int dDPV::N;

/**
 * initialize the static variables
 * @param M_in nr of sp orbs
 * @param N_in nr of particles
 */
void dDPV::init(int M_in,int N_in){

   M = M_in;
   N = N_in;

}

/**
 * standard constructor on a dDPM input, diagonalizes a dDPM and puts the eigenvalues in the allocated dDPV object
 */
dDPV::dDPV(dDPM &ddpm) {

   ddpv = new BlockVector<rxTPM> * [M];

#pragma omp parallel for
   for(int l = 0;l < M;++l)
      ddpv[l] = new BlockVector<rxTPM> (ddpm[l]);
   
}

/**
 * copy constructor:
 * @param ddpv_c object that will be copied into this.
 */
dDPV::dDPV(const dDPV &ddpv_c) { 
   
   ddpv = new BlockVector<rxTPM> * [M];

#pragma omp parallel for
   for(int l = 0;l < M;++l)
      ddpv[l] = new BlockVector<rxTPM> (ddpv_c[l]);
}

/**
 * destructor
 */
dDPV::~dDPV(){

#pragma omp parallel for
   for(int l = 0;l < M;++l)
      delete ddpv[l];

   delete [] ddpv;
   
}


ostream &operator<<(ostream &output,const dDPV &ddpv_p){

   for(int l = 0;l < ddpv_p.gM();++l){

      output << endl;
      output << l << endl;
      output << endl;

      output << ddpv_p[l] << endl;

   }

   return output;

}

/**
 * @return number of particles
 */
int dDPV::gN() const{

   return N;

}

/**
 * @return number of sp orbitals
 */
int dDPV::gM() const{

   return M;

}

/**
 * access to the individual Vector<rxTPM> objects
 * @param l the index of the object you want
 */
BlockVector<rxTPM> &dDPV::operator[](int l){

   return *ddpv[l];

}

/**
 * access to the individual Vector<rxTPM> objects: the const version
 * @param l the index of the object you want
 */
const BlockVector<rxTPM> &dDPV::operator[](int l) const{

   return *ddpv[l];

}

/**
 * overload equality operator
 * @param ddpv_c dDPV object to be copied into this
 */
dDPV &dDPV::operator=(const dDPV &ddpv_c){

#pragma omp parallel for
   for(int l = 0;l < M;++l)
      *ddpv[l] = ddpv_c[l];

   return *this;

}

/**
 * overload equality operator on a double, set all the numbers in the vector equal to parameter a
 * @param a the number
 */
dDPV &dDPV::operator=(double a){

#pragma omp parallel for
   for(int l = 0;l < M;++l)
      *ddpv[l] = a;

   return *this;

}

/**
 * overload += operator
 * @param ddpv_p the dDPV object to be added to this
 */
dDPV &dDPV::operator+=(const dDPV &ddpv_p){

#pragma omp parallel for
   for(int l = 0;l < M;++l)
      *ddpv[l] += ddpv_p[l];

   return *this;

}

/**
 * overload -= operator
 * @param ddpv_m the dDPV object to be deducted from this
 */
dDPV &dDPV::operator-=(const dDPV &ddpv_m){

#pragma omp parallel for
   for(int l = 0;l < M;++l)
      *ddpv[l] -= ddpv_m[l];

   return *this;

}

/**
 * add the ddpv ddpv_p times the constant alpha to this
 * @param alpha the constant to multiply the ddpv_p with
 * @param ddpv_p the dDPV to be multiplied by alpha and added to this
 */
dDPV &dDPV::daxpy(double alpha,const dDPV &ddpv_p){

#pragma omp parallel for
   for(int l = 0;l < M;++l)
      ddpv[l]->daxpy(alpha,ddpv_p[l]);

   return *this;

}

/**
 * @return the sum of all the elements in the vector
 */
double dDPV::sum() const{

   double ward = 0.0;

#pragma omp parallel for reduction(+:ward)
   for(int l = 0;l < M;++l)
      ward += ddpv[l]->sum();

   return ward;

}

/**
 * @return the logarithm of the product of all the elements in the vector (so the sum of all the logarithms)
 */
double dDPV::log_product() const {

   double ward = 0.0;

#pragma omp parallel for reduction(+:ward)
   for(int l = 0;l < M;++l)
      ward += ddpv[l]->log_product();

   return ward;

}

/**
 * @return inproduct of (*this) dDPV with ddpv_i
 * @param ddpv_i input dDPV
 */
double dDPV::ddot(const dDPV &ddpv_i) const{

   double ward = 0.0;

#pragma omp parallel for reduction(+:ward)
   for(int l = 0;l < M;++l)
      ward += ddpv[l]->ddot(ddpv_i[l]);

   return ward;

}

/**
 * Scale the ddpv (*this) with parameter alpha
 * @param alpha scalefactor
 */
void dDPV::dscal(double alpha){

#pragma omp parallel for
   for(int l = 0;l < M;++l)
      ddpv[l]->dscal(alpha);

}

/**
 * @return the minimal element present in this dDPV object.
 * watch out, only works when dDPV is filled with the eigenvalues of a diagonalized Matrix object
 */
double dDPV::min() const {

   double ward = ddpv[0]->min();

   for(int l = 1;l < M;++l)
      if(ddpv[l]->min() < ward)
         ward = ddpv[l]->min();

   return ward;

}

/**
 * @return the maximal element present in this dDPV object.
 * watch out, only works when dDPV is filled with the eigenvalues of a diagonalized Matrix object
 */
double dDPV::max() const {

   double ward = ddpv[0]->max();

   for(int l = 1;l < M;++l)
      if(ddpv[l]->max() > ward)
         ward = ddpv[l]->max();

   return ward;

}

/**
 * Function used by the EIG::lsfunc function, used in the line search algorithm.
 * @param alpha steplenght in the search direction.
 */
double dDPV::lsfunc(double alpha) const {

   double ward = 0;

#pragma omp parallel for reduction(+:ward)
   for(int l = 0;l < M;++l)
      ward += ddpv[l]->lsfunc(alpha);

   return ward;

}
