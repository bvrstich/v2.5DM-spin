#include <iostream>
#include <cmath>
#include <fstream>

using std::ostream;
using std::ofstream;
using std::ifstream;
using std::endl;
using std::ios;

#include "include.h"

int dPHHV::M;
int dPHHV::N;

/**
 * initialize the static variables
 * @param M_in nr of sp orbs
 * @param N_in nr of particles
 */
void dPHHV::init(int M_in,int N_in){

   M = M_in;
   N = N_in;

}

/**
 * standard constructor on a dPHHM input, diagonalizes a dPHHM and puts the eigenvalues in the allocated dPHHV object
 */
dPHHV::dPHHV(dPHHM &dphhm) {

   dphhv = new BlockVector<rxPHM> * [M];

#pragma omp parallel for
   for(int l = 0;l < M;++l)
      dphhv[l] = new BlockVector<rxPHM> (dphhm[l]);
   
}

/**
 * copy constructor:
 * @param dphhv_c object that will be copied into this.
 */
dPHHV::dPHHV(const dPHHV &dphhv_c) { 
   
   dphhv = new BlockVector<rxPHM> * [M];

#pragma omp parallel for
   for(int l = 0;l < M;++l)
      dphhv[l] = new BlockVector<rxPHM> (dphhv_c[l]);
}

/**
 * destructor
 */
dPHHV::~dPHHV(){

   for(int l = 0;l < M;++l)
      delete dphhv[l];

   delete [] dphhv;
   
}


ostream &operator<<(ostream &output,const dPHHV &dphhv_p){

   for(int l = 0;l < dphhv_p.gM();++l){

      output << endl;
      output << l << endl;
      output << endl;

      output << dphhv_p[l] << endl;

   }

   return output;

}

/**
 * @return number of particles
 */
int dPHHV::gN() const{

   return N;

}

/**
 * @return number of sp orbitals
 */
int dPHHV::gM() const{

   return M;

}

/**
 * access to the individual Vector<rxPHM> objects
 * @param l the index of the object you want
 */
BlockVector<rxPHM> &dPHHV::operator[](int l){

   return *dphhv[l];

}

/**
 * access to the individual Vector<rxPHM> objects: the const version
 * @param l the index of the object you want
 */
const BlockVector<rxPHM> &dPHHV::operator[](int l) const{

   return *dphhv[l];

}

/**
 * overload equality operator
 * @param dphhv_c dPHHV object to be copied into this
 */
dPHHV &dPHHV::operator=(const dPHHV &dphhv_c){

#pragma omp parallel for
   for(int l = 0;l < M;++l)
      *dphhv[l] = dphhv_c[l];

   return *this;

}

/**
 * overload equality operator on a double, set all the numbers in the vector equal to parameter a
 * @param a the number
 */
dPHHV &dPHHV::operator=(double a){

#pragma omp parallel for
   for(int l = 0;l < M;++l)
      *dphhv[l] = a;

   return *this;

}

/**
 * overload += operator
 * @param dphhv_p the dPHHV object to be added to this
 */
dPHHV &dPHHV::operator+=(const dPHHV &dphhv_p){

#pragma omp parallel for
   for(int l = 0;l < M;++l)
      *dphhv[l] += dphhv_p[l];

   return *this;

}

/**
 * overload -= operator
 * @param dphhv_m the dPHHV object to be deducted from this
 */
dPHHV &dPHHV::operator-=(const dPHHV &dphhv_m){

#pragma omp parallel for
   for(int l = 0;l < M;++l)
      *dphhv[l] -= dphhv_m[l];

   return *this;

}

/**
 * add the dphhv dphhv_p times the constant alpha to this
 * @param alpha the constant to multiply the dphhv_p with
 * @param dphhv_p the dPHHV to be multiplied by alpha and added to this
 */
dPHHV &dPHHV::daxpy(double alpha,const dPHHV &dphhv_p){

#pragma omp parallel for
   for(int l = 0;l < M;++l)
      dphhv[l]->daxpy(alpha,dphhv_p[l]);

   return *this;

}

/**
 * @return the sum of all the elements in the vector
 */
double dPHHV::sum() const{

   double ward = 0.0;

#pragma omp parallel for reduction(+:ward)
   for(int l = 0;l < M;++l)
      ward += dphhv[l]->sum();

   return ward;

}

/**
 * @return the logarithm of the product of all the elements in the vector (so the sum of all the logarithms)
 */
double dPHHV::log_product() const {

   double ward = 0.0;

#pragma omp parallel for reduction(+:ward)
   for(int l = 0;l < M;++l)
      ward += dphhv[l]->log_product();

   return ward;

}

/**
 * @return inproduct of (*this) dPHHV with dphhv_i
 * @param dphhv_i input dPHHV
 */
double dPHHV::ddot(const dPHHV &dphhv_i) const{

   double ward = 0.0;

#pragma omp parallel for reduction(+:ward)
   for(int l = 0;l < M;++l)
      ward += dphhv[l]->ddot(dphhv_i[l]);

   return ward;

}

/**
 * Scale the dphhv (*this) with parameter alpha
 * @param alpha scalefactor
 */
void dPHHV::dscal(double alpha){

#pragma omp parallel for
   for(int l = 0;l < M;++l)
      dphhv[l]->dscal(alpha);

}

/**
 * @return the minimal element present in this dPHHV object.
 * watch out, only works when dPHHV is filled with the eigenvalues of a diagonalized Matrix object
 */
double dPHHV::min() const {

   double ward = dphhv[0]->min();

   for(int l = 1;l < M;++l)
      if(dphhv[l]->min() < ward)
         ward = dphhv[l]->min();

   return ward;

}

/**
 * @return the maximal element present in this dPHHV object.
 * watch out, only works when dPHHV is filled with the eigenvalues of a diagonalized Matrix object
 */
double dPHHV::max() const {

   double ward = dphhv[0]->max();

   for(int l = 1;l < M;++l)
      if(dphhv[l]->max() > ward)
         ward = dphhv[l]->max();

   return ward;

}

/**
 * Function used by the EIG::lsfunc function, used in the line search algorithm.
 * @param alpha steplenght in the search direction.
 */
double dPHHV::lsfunc(double alpha) const {

   double ward = 0;

#pragma omp parallel for reduction(+:ward)
   for(int l = 0;l < M;++l)
      ward += dphhv[l]->lsfunc(alpha);

   return ward;

}
