#include <iostream>
#include <fstream>
#include <cmath>

using std::ostream;
using std::endl;

#include "include.h"

int EIG::M;
int EIG::N;
int EIG::dim;

/**
 * initialize the statics
 * @param M_in nr of sp orbs
 * @param N_in nr of particles
 */
void EIG::init(int M_in,int N_in){

   M = M_in;
   N = N_in;

   dim = M*(2*M - 2)*(2*M - 1);

}

/**
 * standard constructor with initialization on the eigenvalues of a SUP object.
 * @param X input SUP object that will be destroyed after this function is called. The eigenvectors
 * of the matrix will be stored in the columns of the original SUP matrix.
 */
EIG::EIG(SUP &X){

   v_I1 = new dDPV(X.gI1());

}

/**
 * Copy constructor\n
 * allocates the memory for the eigenvalues of a SUP object and copies the content of eig_c into it.
 * @param eig_c The input EIG that will be copied into this.
 */
EIG::EIG(const EIG &eig_c){

   v_I1 = new dDPV(eig_c.gv_I1());

}

/**
 * overload equality operator
 * @param eig_c object that will be copied into this.
 */
EIG &EIG::operator=(const EIG &eig_c){

   *v_I1 = eig_c.gv_I1();

   return *this;

}

/**
 * Destructor, deallocation of the memory
 */
EIG::~EIG(){

   delete v_I1;

}

ostream &operator<<(ostream &output,const EIG &eig_p){

   std::cout << eig_p.gv_I1() << std::endl;

   return output;

}

/**
 * @return nr of particles
 */
int EIG::gN() const{

   return N;

}

/**
 * @return nr of spatial states
 */
int EIG::gM() const{

   return M;

}

/** 
 * get the dDPV object containing the eigenvalues of the dDPM block I1
 * @return a dDPV object containing the desired eigenvalues
 */
dDPV &EIG::gv_I1(){

   return *v_I1;

}

/** 
 * the const version
 * get the dDPV object containing the eigenvalues of the dDPM block I1
 * @return a dDPV object containing the desired eigenvalues
 */
const dDPV &EIG::gv_I1() const{

   return *v_I1;

}

/**
 * @return total dimension of the EIG object
 */
int EIG::gdim() const{

   return dim;

}

/**
 * @return the minimal element present in this EIG object.
 * watch out, only works when EIG is filled with the eigenvalues of a diagonalized SUP matrix
 */
double EIG::min() const{

   //lowest eigenvalue of P block
   double ward = v_I1->min();

   return ward;

}

/**
 * @return the maximum element present in this EIG object.
 * watch out, only works when EIG is filled with the eigenvalues of a diagonalized SUP matrix
 */
double EIG::max() const{

   //highest eigenvalue of P block
   double ward = v_I1->max();

   return ward;

}

/**
 * @param alpha step length along the Newton direction
 * @return The line search function, gradient of the potential in the Newton direction as a function of the step length alpha
 */
double EIG::lsfunc(double alpha) const{

   double ward = v_I1->lsfunc(alpha);

   return ward;

}