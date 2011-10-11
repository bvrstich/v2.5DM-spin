#include <iostream>
#include <cmath>
#include <fstream>

using std::ostream;
using std::ofstream;
using std::ifstream;
using std::endl;
using std::ios;

#include "include.h"

int dPPHV::M;
int dPPHV::N;

/**
 * initialize the static variables
 * @param M_in nr of sp orbs
 * @param N_in nr of particles
 */
void dPPHV::init(int M_in,int N_in){

   M = M_in;
   N = N_in;

}

/**
 * standard constructor on a dPPHM input, diagonalizes a dPPHM and puts the eigenvalues in the allocated dPPHV object
 */
dPPHV::dPPHV(dPPHM &dpphm) {

   dpphv = new BlockVector<xTPM> * [M];

#pragma omp parallel for
   for(int l = 0;l < M;++l)
      dpphv[l] = new BlockVector<xTPM> (dpphm[l]);
   
}

/**
 * copy constructor:
 * @param dpphv_c object that will be copied into this.
 */
dPPHV::dPPHV(const dPPHV &dpphv_c) { 
   
   dpphv = new BlockVector<xTPM> * [M];

#pragma omp parallel for
   for(int l = 0;l < M;++l)
      dpphv[l] = new BlockVector<xTPM> (dpphv_c[l]);
}

/**
 * destructor
 */
dPPHV::~dPPHV(){

#pragma omp parallel for
   for(int l = 0;l < M;++l)
      delete dpphv[l];

   delete [] dpphv;
   
}


ostream &operator<<(ostream &output,const dPPHV &dpphv_p){

   for(int l = 0;l < dpphv_p.gM();++l){

      output << endl;
      output << l << endl;
      output << endl;

      output << dpphv_p[l] << endl;

   }

   return output;

}

/**
 * @return number of particles
 */
int dPPHV::gN() const{

   return N;

}

/**
 * @return number of sp orbitals
 */
int dPPHV::gM() const{

   return M;

}

/**
 * access to the individual Vector<xTPM> objects
 * @param l the index of the object you want
 */
BlockVector<xTPM> &dPPHV::operator[](int l){

   return *dpphv[l];

}

/**
 * access to the individual Vector<xTPM> objects: the const version
 * @param l the index of the object you want
 */
const BlockVector<xTPM> &dPPHV::operator[](int l) const{

   return *dpphv[l];

}

/**
 * overload equality operator
 * @param dpphv_c dPPHV object to be copied into this
 */
dPPHV &dPPHV::operator=(const dPPHV &dpphv_c){

#pragma omp parallel for
   for(int l = 0;l < M;++l)
      *dpphv[l] = dpphv_c[l];

   return *this;

}

/**
 * overload equality operator on a double, set all the numbers in the vector equal to parameter a
 * @param a the number
 */
dPPHV &dPPHV::operator=(double a){

#pragma omp parallel for
   for(int l = 0;l < M;++l)
      *dpphv[l] = a;

   return *this;

}

/**
 * overload += operator
 * @param dpphv_p the dPPHV object to be added to this
 */
dPPHV &dPPHV::operator+=(const dPPHV &dpphv_p){

#pragma omp parallel for
   for(int l = 0;l < M;++l)
      *dpphv[l] += dpphv_p[l];

   return *this;

}

/**
 * overload -= operator
 * @param dpphv_m the dPPHV object to be deducted from this
 */
dPPHV &dPPHV::operator-=(const dPPHV &dpphv_m){

#pragma omp parallel for
   for(int l = 0;l < M;++l)
      *dpphv[l] -= dpphv_m[l];

   return *this;

}

/**
 * add the dpphv dpphv_p times the constant alpha to this
 * @param alpha the constant to multiply the dpphv_p with
 * @param dpphv_p the dPPHV to be multiplied by alpha and added to this
 */
dPPHV &dPPHV::daxpy(double alpha,const dPPHV &dpphv_p){

#pragma omp parallel for
   for(int l = 0;l < M;++l)
      dpphv[l]->daxpy(alpha,dpphv_p[l]);

   return *this;

}

/**
 * @return the sum of all the elements in the vector
 */
double dPPHV::sum() const{

   double ward = 0.0;

#pragma omp parallel for reduction(+:ward)
   for(int l = 0;l < M;++l)
      ward += dpphv[l]->sum();

   return ward;

}

/**
 * @return the logarithm of the product of all the elements in the vector (so the sum of all the logarithms)
 */
double dPPHV::log_product() const {

   double ward = 0.0;

#pragma omp parallel for reduction(+:ward)
   for(int l = 0;l < M;++l)
      ward += dpphv[l]->log_product();

   return ward;

}

/**
 * @return inproduct of (*this) dPPHV with dpphv_i
 * @param dpphv_i input dPPHV
 */
double dPPHV::ddot(const dPPHV &dpphv_i) const{

   double ward = 0.0;

#pragma omp parallel for reduction(+:ward)
   for(int l = 0;l < M;++l)
      ward += dpphv[l]->ddot(dpphv_i[l]);

   return ward;

}

/**
 * Scale the dpphv (*this) with parameter alpha
 * @param alpha scalefactor
 */
void dPPHV::dscal(double alpha){

#pragma omp parallel for
   for(int l = 0;l < M;++l)
      dpphv[l]->dscal(alpha);

}

/**
 * @return the minimal element present in this dPPHV object.
 * watch out, only works when dPPHV is filled with the eigenvalues of a diagonalized Matrix object
 */
double dPPHV::min() const {

   double ward = dpphv[0]->min();

   for(int l = 1;l < M;++l)
      if(dpphv[l]->min() < ward)
         ward = dpphv[l]->min();

   return ward;

}

/**
 * @return the maximal element present in this dPPHV object.
 * watch out, only works when dPPHV is filled with the eigenvalues of a diagonalized Matrix object
 */
double dPPHV::max() const {

   double ward = dpphv[0]->max();

   for(int l = 1;l < M;++l)
      if(dpphv[l]->max() > ward)
         ward = dpphv[l]->max();

   return ward;

}

/**
 * Function used by the EIG::lsfunc function, used in the line search algorithm.
 * @param alpha steplenght in the search direction.
 */
double dPPHV::lsfunc(double alpha) const {

   double ward = 0;

#pragma omp parallel for reduction(+:ward)
   for(int l = 0;l < M;++l)
      ward += dpphv[l]->lsfunc(alpha);

   return ward;

}
