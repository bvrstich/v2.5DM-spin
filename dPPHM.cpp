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

int dPPHM::M;
int dPPHM::N;

/**
 * initialize the static variables
 * @param M_in dimension of spatial sp space
 * @param N_in nr of particles
 */
void dPPHM::init(int M_in,int N_in){

   M = M_in;
   N = N_in;

}

/**
 * standard constructor: constructs M xTPM object with parameter l = 0 -> M-1
 */
dPPHM::dPPHM() {

   dpphm = new xTPM * [M];

   for(int l = 0;l < M;++l)
      dpphm[l] = new xTPM();

}

/**
 * copy constructor: constructs an array and copies the content of W in it.
 * @param W object that will be copied into this.
 */
dPPHM::dPPHM(const dPPHM &W) { 

   dpphm = new xTPM * [M];

   for(int l = 0;l < M;++l)
      dpphm[l] = new xTPM(W[l]);

}

/**
 * destructor
 */
dPPHM::~dPPHM(){

   for(int l = 0;l < M;++l)
      delete dpphm[l];

   delete [] dpphm;

}

/**
 * @return nr of particles
 */
int dPPHM::gN() const {

   return N;

}

/**
 * @return dimension of sp space
 */
int dPPHM::gM() const {

   return M;

}

/**
 * acces to the individual xTPM objects
 * @param l the specific xTPM object you want
 * @return the xTPM object with parameter l
 */
xTPM &dPPHM::operator[](int l){

   return *dpphm[l];

}

/**
 * acces to the individual xTPM objects: the const version
 * @param l the specific xTPM object you want
 * @return the xTPM object with parameter l
 */
const xTPM &dPPHM::operator[](int l) const{

   return *dpphm[l];

}

/**
 * access the numbers in sp mode
 * @param l blockindex
 * @param S dp spin
 * @param S_ab intermediate spin for a and b
 * @param a first spatial index for the row
 * @param b second spatial index for the row
 * @param S_cd intermediate spin for c and d
 * @param c first spatial index for the column
 * @param d second spatial index for the column
 */
double dPPHM::operator()(int l,int S,int S_ab,int a,int b,int S_cd,int c,int d) const {

   int phase_i = get_inco(S,S_ab,a,b);

   if(phase_i == 0)
      return 0.0;

   int phase_j = get_inco(S,S_cd,c,d);

   if(phase_j == 0)
      return 0.0;

   return phase_i * phase_j * (*this)[l](S,xTPM::gs2t(S,S_ab,a,b),xTPM::gs2t(S,S_cd,c,d));

}

/**
 * @param S dp spin
 * @param S_ab intermediate spin
 * @param a first sp index
 * @param b second sp index
 * @return the right phase for this order of sp indices as a function of my basis.
 */
int dPPHM::get_inco(int S,int S_ab,int a,int b){

   if(S == 0){//S == 1/2

      if(a == b && S_ab == 1)
         return 0;

      if(a > b)
         return 1 - 2*S_ab;
      else
         return 1;

   }
   else{//S == 3/2

      //intermediate spin can never be 0
      if(S_ab == 0)
         return 0;

      //totally antisymmetrical
      if(a == b)
         return 0;

      if(a > b)
         return -1;
      else
         return 1;

   }

}

ostream &operator<<(ostream &output,const dPPHM &dpphm_p){

   for(int l = 0;l < dpphm_p.gM();++l){

      output << std::endl;
      output << "l = \t" << l << std::endl;
      output << std::endl;

      output << dpphm_p[l] << endl;

   }

   return output;

}

/**
 * @return the trace of the dPPHM object
 */
double dPPHM::trace() const{

   double ward = 0.0;

   for(int l = 0;l < M;++l)
      ward += dpphm[l]->trace();

   return ward;

}

/**
 * Scale the dPPHM with parameter alpha
 * @param alpha scalefactor
 */
void dPPHM::dscal(double alpha){

   for(int l = 0;l < M;++l)
      dpphm[l]->dscal(alpha);

}

/**
 * overload the equality operator
 * @param dpphm_c The dPPHM you want to be copied into this
 */
dPPHM &dPPHM::operator=(const dPPHM &dpphm_c){

   for(int l = 0;l < M;++l)
      *dpphm[l] = dpphm_c[l];

   return *this;

}

/**
 * Make all the number in your dPPHM equal to the number a, e.g. usefull for initialization (dPPHM W = 0)
 * @param a the number
 */
dPPHM &dPPHM::operator=(double a){

   for(int l = 0;l < M;++l)
      *dpphm[l] = a;

   return *this;

}

/**
 * overload the += operator for dPPHM's
 * @param dpphm_pl The dPPHM you want to add to this
 */
dPPHM &dPPHM::operator+=(const dPPHM &dpphm_pl){

   for(int l = 0;l < M;++l)
      *dpphm[l] += dpphm_pl[l];

   return *this;

}

/**
 * overload the -= operator for dPPHM's
 * @param dpphm_m The dPPHM you want to deduct from this
 */
dPPHM &dPPHM::operator-=(const dPPHM &dpphm_m){

   for(int l = 0;l < M;++l)
      *dpphm[l] -= dpphm_m[l];

   return *this;

}

/**
 * add the dpphm dpphm_pl times the constant alpha to *this
 * @param alpha the constant to multiply the dpphm_pl with
 * @param dpphm_pl the dPPHM to be multiplied by alpha and added to this
 */
dPPHM &dPPHM::daxpy(double alpha,const dPPHM &dpphm_pl){

   for(int l = 0;l < M;++l)
      dpphm[l]->daxpy(alpha,dpphm_pl[l]);

   return *this;

}

/**
 * Multiply symmetric dpphm object left en right with symmetric dpphm map to 
 * form another symmetric dpphm and put it in (*this): this = map*object*map
 * @param map dpphm that will be multiplied to the left en to the right of dpphm object
 * @param object central dpphm
 */
void dPPHM::L_map(const dPPHM &map,const dPPHM &object){

   for(int l = 0;l < M;++l)
      dpphm[l]->L_map(map[l],object[l]);

}

/**
 * dPPHM product of two general matrices A en B, put result in this
 * @param A left dpphm
 * @param B right dpphm
 */
dPPHM &dPPHM::mprod(const dPPHM &A,const dPPHM &B){

   for(int l = 0;l < M;++l)
      dpphm[l]->mprod(A[l],B[l]);

   return *this;

}

/**
 * Take the square root out of the positive semidefinite dpphm, destroys original dpphm, square root will be put in (*this)
 * @param option = 1, positive square root, = -1, negative square root.
 */
void dPPHM::sqrt(int option){

   for(int l = 0;l < M;++l)
      dpphm[l]->sqrt(option);

}

/**
 * @return inproduct of (*this) dpphm with dpphm_i, defined as Tr (A B)
 * @param dpphm_i input dpphm
 */
double dPPHM::ddot(const dPPHM &dpphm_i) const{

   double ward = 0.0;

   for(int l = 0;l < M;++l)
      ward += dpphm[l]->ddot(dpphm_i[l]);

   return ward;

}

/**
 * Invert positive semidefinite symmetric dpphm is stored in (*this), original dpphm (*this) is destroyed
 */
void dPPHM::invert(){

   for(int l = 0;l < M;++l)
      dpphm[l]->invert();

}

/**
 * copy upper in lower part of dPPHM object
 */
void dPPHM::symmetrize(){

   for(int l = 0;l < M;++l)
      dpphm[l]->symmetrize();

}

/**
 * Fill the matrix with random numbers.
 */
void dPPHM::fill_Random(){

   for(int l = 0;l < M;++l)
      dpphm[l]->fill_Random();

}