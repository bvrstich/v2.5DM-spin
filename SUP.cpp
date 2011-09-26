#include <iostream>
#include <cmath>
#include <fstream>

using std::ostream;

#include "include.h"

int SUP::M;
int SUP::N;
int SUP::dim;

/**
 * initialize the statics
 * @param M_in nr of sp orbs
 * @param N_in nr of particles
 */
void SUP::init(int M_in,int N_in){

   M = M_in;
   N = N_in;

   dim = M*(2*M - 2)*(2*M - 1);

}

/**
 * standard constructor\n
 * Allocates two dDPM matrices, I1 and Q2 and one dPPHM matrix I2
 */
SUP::SUP(){

   I1 = new dDPM();

}

/**
 * standard constructor\n
 * Allocates two dDPM's and a dPPHM, then copies the content of input SUP X into it.
 * @param X_c input SUP
 */
SUP::SUP(const SUP &X_c){

   I1 = new dDPM(X_c.gI1());

}

/**
 * Destructor
 */
SUP::~SUP(){

   delete I1;

}

/**
 * Overload += operator
 * @param X_pl The SUP matrix that has to be added to this
 */
SUP &SUP::operator+=(const SUP &X_pl){

   *I1 += X_pl.gI1();

   return *this;

}

/**
 * Overload -= operator
 * @param X_pl The SUP that will be deducted from this
 */
SUP &SUP::operator-=(const SUP &X_pl){

   *I1 -= X_pl.gI1();

   return *this;

}

/**
 * Overload equality operator, copy X_c into this
 * @param X_c SUP to be copied into this
 */
SUP &SUP::operator=(const SUP &X_c){

   *I1 = X_c.gI1();

   return *this;

}

/**
 * overload operator = number, all the blockmatrices in SUP are put equal to the number a.
 * e.g. X = 0 makes all the Matrix elements zero.
 * @param a the number
 */
SUP &SUP::operator=(double &a){

   *I1 = a;

   return *this;

}

/**
 * @return pointer to the dDPM object containing the I1 block
 */
dDPM &SUP::gI1(){

   return *I1;

}

/**
 * the const version
 * @return pointer to the dDPM object containing the I1 block
 */
const dDPM &SUP::gI1() const{

   return *I1;

}

ostream &operator<<(ostream &output,const SUP &X_p){

   output << (X_p.gI1()) << std::endl;

   return output;

}

/**
 * Fill the SUP matrix with random elements, watch out, not necessarily positive definite
 */
void SUP::fill_Random(){

   I1->fill_Random();

}

/**
 * @return number of particles
 */
int SUP::gN() const {

   return N;

}

/**
 * @return nr of spatial orbs
 */
int SUP::gM() const{

   return M;

}

/**
 * @return total dimension of SUP (carrier) space
 */
int SUP::gdim() const{

   return dim;

}

/**
 * @param X_i input SUP_PQ X_i
 * @return inproduct between this and input matrix X_i, defined as Tr(this X_i)
 */
double SUP::ddot(const SUP &X_i) const{

   double ward = 0.0;

   ward += I1->ddot(X_i.gI1());

   return ward;

}

/**
 * Invert all the Matrices in SUP and put it in this, watch out, destroys original matrices.
 * Makes use of cholesky decomposition, so only positive definite matrices can be used as input!
 */
void SUP::invert(){

   I1->pseudo_invert();

}

/**
 * Scale all the elements in the blockmatrix with parameter alpha
 * @param alpha the scalefactor
 */
void SUP::dscal(double alpha){

   I1->dscal(alpha);

}

/**
 * Take the square root out of the positive semidefinite SUP matrix this and put it in this, watch out, original matrix is destroyed.
 * @param option = 1 positive square root, = -1 negative square root.
 */
void SUP::sqrt(int option){

   I1->pseudo_sqrt(option);

}

/**
 * Multiply symmetric SUP blockmatrix object left en right with symmetric SUP blockmatrix map to 
 * form another symmetric matrix and put it in (*this): this = map*object*map
 * @param map SUP that will be multiplied to the left en to the right of matrix object
 * @param object central SUP
 */
void SUP::L_map(const SUP &map,const SUP &object){

   I1->L_map(map.gI1(),object.gI1());

}

/**
 * add the SUP X_p times the constant alpha to this
 * @param alpha the constant to multiply the X_p with
 * @param X_p the SUP to be multiplied by alpha and added to (*this)
 */
void SUP::daxpy(double alpha,const SUP &X_p){

   I1->daxpy(alpha,X_p.gI1());

}

/**
 * @return trace of the SUP matrix, defined as sum of the traces of the separate carrierspace matrices
 */
double SUP::trace() const{

   double ward = 0.0;

   ward += I1->trace();

   return ward;

}

/**
 * General matrixproduct between two SUP matrices, act with Matrix::mprod on every block
 * @param A left hand matrix
 * @param B right hand matrix
 * @return The product AB
 */
SUP &SUP::mprod(const SUP &A,const SUP &B){

   I1->mprod(A.gI1(),B.gI1());

   return *this;

}

/**
 * Fill the SUP matrix (*this) with two dDPM matrices like: this = diag[I^l(W)  Q^l(W)]
 * @param W input dDPM
 */
void SUP::fill(const dDPM &W){

   *I1 = W;

}
