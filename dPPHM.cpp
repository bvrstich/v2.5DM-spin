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

/**
 * map a dDPM on a dPPHM using the I2 map
 * @param ddpm input dDPM
 */
void dPPHM::I(const dDPM &ddpm){

   int a,b,c,d;

   int S_ab,S_cd;

   TPM tpm;
   tpm.bar(1.0/(N - 2.0),ddpm);

   for(int l = 0;l < M;++l){

      for(int S = 0;S < 2;++S){

         for(int i = 0;i < dpphm[l]->gdim(S);++i){

            S_ab = xTPM::gt2s(S,i,0);

            a = xTPM::gt2s(S,i,1);
            b = xTPM::gt2s(S,i,2);

            for(int j = i;j < dpphm[l]->gdim(S);++j){

               S_cd = xTPM::gt2s(S,j,0);

               c = xTPM::gt2s(S,j,1);
               d = xTPM::gt2s(S,j,2);

               (*this)[l](S,i,j) = 0.0;

               for(int S_ = 0;S_ < 2;++S_)
                  (*this)[l](S,i,j) -= (2* (S_ + 0.5) + 1.0) * Tools::gC(S,S_,S_ab,S_cd) * ddpm(l,S_,S_ab,a,b,S_cd,c,d);

               if(S_ab == S_cd)
                  (*this)[l](S,i,j) += tpm(S_ab,a,b,c,d);

            }
         }
      }

   }

   this->symmetrize();

}

/**
 * map a dDPM on a dPPHM using the Q1 map
 * @param ddpm input dDPM
 */
void dPPHM::Q(const dDPM &ddpm){

   int a,b,c,d;

   int S_ab,S_cd;

   int sign_ab,sign_cd;

   double norm_ab,norm_cd;

   TPM tpm;
   tpm.bar(1.0/(N - 2.0),ddpm);

   SPM spm;
   spm.bar(1.0/(N - 1.0),tpm);

   double ward = 2.0 * ddpm.trace() / (N*(N - 1.0)*(N - 2.0));

   for(int l = 0;l < M;++l){

      //first S = 1/2
      for(int i = 0;i < dpphm[l]->gdim(0);++i){

         S_ab = xTPM::gt2s(0,i,0);

         a = xTPM::gt2s(0,i,1);
         b = xTPM::gt2s(0,i,2);

         norm_ab = 1.0;

         if(a == b)
            norm_ab /= std::sqrt(2.0);

         sign_ab = 1 - 2*S_ab;

         for(int j = i;j < dpphm[l]->gdim(0);++j){

            S_cd = xTPM::gt2s(0,j,0);

            c = xTPM::gt2s(0,j,1);
            d = xTPM::gt2s(0,j,2);

            sign_cd = 1 - 2*S_cd;

            norm_cd = 1.0;

            if(c == d)
               norm_cd /= std::sqrt(2.0);

            //dp term
            (*this)[l](0,i,j) = 0.0;

            for(int S_ = 0;S_ < 2;++S_)
               (*this)[l](0,i,j) += (2* (S_ + 0.5) + 1.0) * Tools::gC(0,S_,S_ab,S_cd) * ddpm(l,S_,S_ab,a,b,S_cd,c,d);

            if(b == d){

               if(a == c){

                  //np_a
                  if(a == l)
                     (*this)[l](0,i,j) += 0.5 * std::sqrt( (2*S_ab + 1.0) * (2.0*S_cd + 1.0) ) * norm_ab * norm_cd * ward;

                  //np_d
                  if(b == l)
                     (*this)[l](0,i,j) += 0.5 * sign_ab * sign_cd * std::sqrt( (2*S_ab + 1.0) * (2.0*S_cd + 1.0) ) * norm_ab * norm_cd * ward;

                  //sp(2)_a
                  if(S_ab == S_cd)
                     (*this)[l](0,i,j) += norm_ab * norm_cd * spm(l,l);

               }

               //sp(1)_d
               if(b == l)
                  (*this)[l](0,i,j) -= norm_ab * norm_cd * sign_ab * sign_cd * 0.5 * std::sqrt( (2*S_ab + 1.0) * (2.0*S_cd + 1.0) ) * spm(a,c);

               //sp(3)_a first part
               if(c == l)
                  (*this)[l](0,i,j) -= norm_ab * norm_cd * 0.5 * std::sqrt( (2*S_ab + 1.0) * (2.0*S_cd + 1.0) ) * spm(a,l);
               
               //sp(3)_c first part
               if(a == l)
                  (*this)[l](0,i,j) -= norm_ab * norm_cd * 0.5 * std::sqrt( (2*S_ab + 1.0) * (2.0*S_cd + 1.0) ) * spm(c,l);

               //tp(2)_d
               double hard = 0.0;

               for(int Z = 0;Z < 2;++Z)
                  hard += (2*Z + 1.0) * Tools::gD(0,Z,S_ab,S_cd) * tpm(Z,a,l,c,l);

               if(a == l)
                  hard *= std::sqrt(2.0);

               if(c == l)
                  hard *= std::sqrt(2.0);

               (*this)[l](0,i,j) += sign_ab * sign_cd * std::sqrt( (2*S_ab + 1.0) * (2.0*S_cd + 1.0) ) * norm_ab * norm_cd * hard;

            }

            if(a == d){

               if(b == c){

                  //np_b
                  if(a == l)
                     (*this)[l](0,i,j) += sign_cd * 0.5 * std::sqrt( (2*S_ab + 1.0) * (2.0*S_cd + 1.0) ) * norm_ab * norm_cd * ward;

                  //np_c
                  if(b == l)
                     (*this)[l](0,i,j) += sign_ab * 0.5 * std::sqrt( (2*S_ab + 1.0) * (2.0*S_cd + 1.0) ) * norm_ab * norm_cd * ward;

                  //sp(2)_b
                  if(S_ab == S_cd)
                     (*this)[l](0,i,j) += sign_ab * norm_ab * norm_cd * spm(l,l);

               }

               //sp(1)_b
               if(a == l)
                  (*this)[l](0,i,j) -= norm_ab * norm_cd * sign_cd * 0.5 * std::sqrt( (2*S_ab + 1.0) * (2.0*S_cd + 1.0) ) * spm(b,c);

               //sp(3)_b second part
               if(c == l)
                  (*this)[l](0,i,j) -= norm_ab * norm_cd * sign_ab * 0.5 * std::sqrt( (2*S_ab + 1.0) * (2.0*S_cd + 1.0) ) * spm(l,b);

               //sp(3)_c second part
               if(b == l)
                  (*this)[l](0,i,j) -= norm_ab * norm_cd * sign_ab * 0.5 * std::sqrt( (2*S_ab + 1.0) * (2.0*S_cd + 1.0) ) * spm(l,c);

               //tp(2)_b
               double hard = 0.0;

               for(int Z = 0;Z < 2;++Z)
                  hard += (2*Z + 1.0) * Tools::gD(0,Z,S_ab,S_cd) * tpm(Z,b,l,c,l);

               if(b == l)
                  hard *= std::sqrt(2.0);

               if(c == l)
                  hard *= std::sqrt(2.0);

               (*this)[l](0,i,j) += sign_cd * std::sqrt( (2*S_ab + 1.0) * (2.0*S_cd + 1.0) ) * norm_ab * norm_cd * hard;

            }

            if(b == c){

               //sp(1)_c
               if(b == l)
                  (*this)[l](0,i,j) -= norm_ab * norm_cd * sign_ab * 0.5 * std::sqrt( (2*S_ab + 1.0) * (2.0*S_cd + 1.0) ) * spm(a,d);

               //sp(3)_a second part
               if(d == l)
                  (*this)[l](0,i,j) -= norm_ab * norm_cd * sign_cd * 0.5 * std::sqrt( (2*S_ab + 1.0) * (2.0*S_cd + 1.0) ) * spm(a,l);

               //sp(3)_d second part
               if(a == l)
                  (*this)[l](0,i,j) -= norm_ab * norm_cd * sign_cd * 0.5 * std::sqrt( (2*S_ab + 1.0) * (2.0*S_cd + 1.0) ) * spm(d,l);

               //tp(2)_c
               double hard = 0.0;

               for(int Z = 0;Z < 2;++Z)
                  hard += (2*Z + 1.0) * Tools::gD(0,Z,S_ab,S_cd) * tpm(Z,a,l,d,l);

               if(a == l)
                  hard *= std::sqrt(2.0);

               if(d == l)
                  hard *= std::sqrt(2.0);

               (*this)[l](0,i,j) += sign_ab * std::sqrt( (2*S_ab + 1.0) * (2.0*S_cd + 1.0) ) * norm_ab * norm_cd * hard;


            }

            if(a == c){

               //sp(1)_a
               if(a == l)
                  (*this)[l](0,i,j) -= norm_ab * norm_cd * 0.5 * std::sqrt( (2*S_ab + 1.0) * (2.0*S_cd + 1.0) ) * spm(b,d);

               //sp(3)_b first part
               if(d == l)
                  (*this)[l](0,i,j) -= norm_ab * norm_cd * sign_ab * sign_cd * 0.5 * std::sqrt( (2*S_ab + 1.0) * (2.0*S_cd + 1.0) ) * spm(b,l);

               //sp(3)_d first part
               if(b == l)
                  (*this)[l](0,i,j) -= norm_ab * norm_cd * sign_ab * sign_cd * 0.5 * std::sqrt( (2*S_ab + 1.0) * (2.0*S_cd + 1.0) ) * spm(d,l);

               //tp(2)_a
               double hard = 0.0;

               for(int Z = 0;Z < 2;++Z)
                  hard += (2*Z + 1.0) * Tools::gD(0,Z,S_ab,S_cd) * tpm(Z,b,l,d,l);

               if(b == l)
                  hard *= std::sqrt(2.0);

               if(d == l)
                  hard *= std::sqrt(2.0);

               (*this)[l](0,i,j) += std::sqrt( (2*S_ab + 1.0) * (2.0*S_cd + 1.0) ) * norm_ab * norm_cd * hard;

            }

            //tp(1)_a
            if(a == l){

               if(l == b)
                  (*this)[l](0,i,j) += 0.5 * norm_ab * std::sqrt( 2.0 * (2*S_ab + 1.0) * (2.0*S_cd + 1.0) ) * tpm(S_cd,l,b,c,d);
               else
                  (*this)[l](0,i,j) += 0.5 * norm_ab * std::sqrt( (2*S_ab + 1.0) * (2.0*S_cd + 1.0) ) * tpm(S_cd,l,b,c,d);

            }

            //tp(1)_b
            if(b == l){

               if(a == l)
                  (*this)[l](0,i,j) += 0.5 * sign_ab * sign_cd * norm_ab * std::sqrt( 2.0 * (2*S_ab + 1.0) * (2.0*S_cd + 1.0) ) * tpm(S_cd,a,l,c,d);
               else
                  (*this)[l](0,i,j) += 0.5 * sign_ab * sign_cd * norm_ab * std::sqrt( (2*S_ab + 1.0) * (2.0*S_cd + 1.0) ) * tpm(S_cd,a,l,c,d);

            }

            //tp(1)_c
            if(c == l){

               if(d == l)
                  (*this)[l](0,i,j) += 0.5 * norm_cd * std::sqrt( 2.0 * (2*S_ab + 1.0) * (2.0*S_cd + 1.0) ) * tpm(S_ab,a,b,l,d);
               else
                  (*this)[l](0,i,j) += 0.5 * norm_cd * std::sqrt( (2*S_ab + 1.0) * (2.0*S_cd + 1.0) ) * tpm(S_ab,a,b,l,d);

            }

            //tp(1)_d
            if(d == l){

               if(c == l)
                  (*this)[l](0,i,j) += 0.5 * norm_cd * sign_ab * sign_cd * std::sqrt( 2.0 * (2*S_ab + 1.0) * (2.0*S_cd + 1.0) ) * tpm(S_ab,a,b,c,l);
               else
                  (*this)[l](0,i,j) += 0.5 * norm_cd * sign_ab * sign_cd * std::sqrt( (2*S_ab + 1.0) * (2.0*S_cd + 1.0) ) * tpm(S_ab,a,b,c,l);

            }

         }
      }

      //S = 3/2 is a lot easier
      for(int i = 0;i < dpphm[l]->gdim(1);++i){

         S_ab = xTPM::gt2s(1,i,0);

         a = xTPM::gt2s(1,i,1);
         b = xTPM::gt2s(1,i,2);

         for(int j = i;j < dpphm[l]->gdim(1);++j){

            S_cd = xTPM::gt2s(1,j,0);

            c = xTPM::gt2s(1,j,1);
            d = xTPM::gt2s(1,j,2);

            //dp term
            (*this)[l](1,i,j) = 0.0;

            for(int S_ = 0;S_ < 2;++S_)
               (*this)[l](1,i,j) += (2* (S_ + 0.5) + 1.0) * Tools::gC(1,S_,S_ab,S_cd) * ddpm(l,S_,S_ab,a,b,S_cd,c,d);

            //sp(2) full
            if(i == j)
               (*this)[l](1,i,j) += spm(l,l);

            if(b == d){

               //tp(2)_d
               double hard = 0.0;

               for(int Z = 0;Z < 2;++Z)
                  hard += (2*Z + 1.0) * Tools::gD(1,Z,S_ab,S_cd) * tpm(Z,a,l,c,l);

               if(a == l)
                  hard *= std::sqrt(2.0);

               if(c == l)
                  hard *= std::sqrt(2.0);

               (*this)[l](1,i,j) += 3.0 * hard;

            }

            if(a == d){

               //tp(2)_b
               double hard = 0.0;

               for(int Z = 0;Z < 2;++Z)
                  hard += (2*Z + 1.0) * Tools::gD(1,Z,S_ab,S_cd) * tpm(Z,b,l,c,l);

               if(b == l)
                  hard *= std::sqrt(2.0);

               if(c == l)
                  hard *= std::sqrt(2.0);

               (*this)[l](1,i,j) -= 3.0 * hard;

            }

            if(b == c){

               //tp(2)_c
               double hard = 0.0;

               for(int Z = 0;Z < 2;++Z)
                  hard += (2*Z + 1.0) * Tools::gD(1,Z,S_ab,S_cd) * tpm(Z,a,l,d,l);

               if(a == l)
                  hard *= std::sqrt(2.0);

               if(d == l)
                  hard *= std::sqrt(2.0);

               (*this)[l](1,i,j) -= 3.0 * hard;


            }

            if(a == c){

               //tp(2)_a
               double hard = 0.0;

               for(int Z = 0;Z < 2;++Z)
                  hard += (2*Z + 1.0) * Tools::gD(1,Z,S_ab,S_cd) * tpm(Z,b,l,d,l);

               if(b == l)
                  hard *= std::sqrt(2.0);

               if(d == l)
                  hard *= std::sqrt(2.0);

               (*this)[l](1,i,j) += 3.0 * hard;

            }

         }
      }

   }

   this->symmetrize();

}
