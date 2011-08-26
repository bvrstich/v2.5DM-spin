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

int dDPM::M;
int dDPM::N;

/**
 * initialize the static variables
 * @param M_in dimension of spatial sp space
 * @param N_in nr of particles
 */
void dDPM::init(int M_in,int N_in){

   M = M_in;
   N = N_in;

}

/**
 * standard constructor: constructs M rTPM object with parameter l = 0 -> M-1
 */
dDPM::dDPM() {

   ddpm = new rTPM * [M];

   for(int l = 0;l < M;++l)
      ddpm[l] = new rTPM(l);

}

/**
 * copy constructor: constructs an array and copies the content of W in it.
 * @param W object that will be copied into this.
 */
dDPM::dDPM(const dDPM &W) { 

   ddpm = new rTPM * [M];

   for(int l = 0;l < M;++l)
      ddpm[l] = new rTPM(W[l]);
   
}

/**
 * destructor
 */
dDPM::~dDPM(){

   for(int l = 0;l < M;++l)
      delete ddpm[l];

   delete [] ddpm;
   
}

/**
 * @return nr of particles
 */
int dDPM::gN() const {

   return N;

}

/**
 * @return dimension of sp space
 */
int dDPM::gM() const {

   return M;

}

/**
 * acces to the individual rTPM objects
 * @param l the specific rTPM object you want
 * @return the rTPM object with parameter l
 */
rTPM &dDPM::operator[](int l){

   return *ddpm[l];

}

/**
 * acces to the individual rTPM objects: the const version
 * @param l the specific rTPM object you want
 * @return the rTPM object with parameter l
 */
const rTPM &dDPM::operator[](int l) const{

   return *ddpm[l];

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
double operator()(int l,int S,int S_ab,int a,int b,int S_cd,int c,int d) const {

}

ostream &operator<<(ostream &output,const dDPM &ddpm_p){

   for(int l = 0;l < ddpm_p.gM();++l){

      output << std::endl;
      output << "l = \t" << l << std::endl;
      output << std::endl;

      output << ddpm_p[l] << endl;

   }

   return output;

}

/**
 * @return the trace of the dDPM object
 */
double dDPM::trace() const{

   double ward = 0.0;

   for(int l = 0;l < M;++l)
      ward += ddpm[l]->trace();

   return ward;

}

/**
 * Scale the dDPM with parameter alpha
 * @param alpha scalefactor
 */
void dDPM::dscal(double alpha){

   for(int l = 0;l < M;++l)
      ddpm[l]->dscal(alpha);

}

/**
 * overload the equality operator
 * @param ddpm_c The dDPM you want to be copied into this
 */
dDPM &dDPM::operator=(const dDPM &ddpm_c){

   for(int l = 0;l < M;++l)
      *ddpm[l] = ddpm_c[l];

   return *this;

}

/**
 * Make all the number in your dDPM equal to the number a, e.g. usefull for initialization (dDPM W = 0)
 * @param a the number
 */
dDPM &dDPM::operator=(double a){

   for(int l = 0;l < M;++l)
      *ddpm[l] = a;

   return *this;

}

/**
 * overload the += operator for dDPM's
 * @param ddpm_pl The dDPM you want to add to this
 */
dDPM &dDPM::operator+=(const dDPM &ddpm_pl){

   for(int l = 0;l < M;++l)
      *ddpm[l] += ddpm_pl[l];

   return *this;

}

/**
 * overload the -= operator for dDPM's
 * @param ddpm_m The dDPM you want to deduct from this
 */
dDPM &dDPM::operator-=(const dDPM &ddpm_m){

   for(int l = 0;l < M;++l)
      *ddpm[l] -= ddpm_m[l];

   return *this;

}

/**
 * add the ddpm ddpm_pl times the constant alpha to *this
 * @param alpha the constant to multiply the ddpm_pl with
 * @param ddpm_pl the dDPM to be multiplied by alpha and added to this
 */
dDPM &dDPM::daxpy(double alpha,const dDPM &ddpm_pl){

   for(int l = 0;l < M;++l)
      ddpm[l]->daxpy(alpha,ddpm_pl[l]);

   return *this;

}

/**
 * Multiply symmetric ddpm object left en right with symmetric ddpm map to 
 * form another symmetric ddpm and put it in (*this): this = map*object*map
 * @param map ddpm that will be multiplied to the left en to the right of ddpm object
 * @param object central ddpm
 */
void dDPM::L_map(const dDPM &map,const dDPM &object){

   for(int l = 0;l < M;++l)
      ddpm[l]->L_map(map[l],object[l]);

}

/**
 * dDPM product of two general matrices A en B, put result in this
 * @param A left ddpm
 * @param B right ddpm
 */
dDPM &dDPM::mprod(const dDPM &A,const dDPM &B){

   for(int l = 0;l < M;++l)
      ddpm[l]->mprod(A[l],B[l]);

   return *this;

}

/**
 * Take the square root out of the positive semidefinite ddpm, destroys original ddpm, square root will be put in (*this)
 * @param option = 1, positive square root, = -1, negative square root.
 */
void dDPM::sqrt(int option){

   for(int l = 0;l < M;++l)
      ddpm[l]->sqrt(option);

}

/**
 * @return inproduct of (*this) ddpm with ddpm_i, defined as Tr (A B)
 * @param ddpm_i input ddpm
 */
double dDPM::ddot(const dDPM &ddpm_i) const{

   double ward = 0.0;

   for(int l = 0;l < M;++l)
      ward += ddpm[l]->ddot(ddpm_i[l]);

   return ward;

}

/**
 * Invert positive semidefinite symmetric ddpm is stored in (*this), original ddpm (*this) is destroyed
 */
void dDPM::invert(){

   for(int l = 0;l < M;++l)
      ddpm[l]->invert();

}

/**
 * copy upper in lower part of dDPM object
 */
void dDPM::symmetrize(){

   for(int l = 0;l < M;++l)
      ddpm[l]->symmetrize();

}

/**
 * Fill the matrix with random numbers.
 */
void dDPM::fill_Random(){

   for(int l = 0;l < M;++l)
      ddpm[l]->fill_Random();

}
