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

int dPPHM_unc::M;
int dPPHM_unc::N;

/**
 * initialize the static variables
 * @param M_in dimension of sp space
 * @param N_in nr of particles
 */
void dPPHM_unc::init(int M_in,int N_in){

   M = M_in;
   N = N_in;

}

/**
 * standard constructor: constructs M xTPM_unc object with parameter l = 0 -> M-1
 */
dPPHM_unc::dPPHM_unc() {

   dpphm = new xTPM_unc * [M];

   for(int l = 0;l < M;++l)
      dpphm[l] = new xTPM_unc();

}

/**
 * copy constructor: constructs an array and copies the content of dpphm_c in it.
 * @param dpphm_c object that will be copied into this.
 */
dPPHM_unc::dPPHM_unc(const dPPHM_unc &dpphm_c) { 

   dpphm = new xTPM_unc * [M];

   for(int l = 0;l < M;++l)
      dpphm[l] = new xTPM_unc(dpphm_c[l]);
   
}

/**
 * destructor
 */
dPPHM_unc::~dPPHM_unc(){

   for(int l = 0;l < M;++l)
      delete dpphm[l];

   delete [] dpphm;
   
}

/**
 * @return nr of particles
 */
int dPPHM_unc::gN() const {

   return N;

}

/**
 * @return dimension of sp space
 */
int dPPHM_unc::gM() const {

   return M;

}

/**
 * acces to the individual xTPM_unc objects
 * @param l the specific xTPM_unc object you want
 * @return the xTPM_unc object with parameter l
 */
xTPM_unc &dPPHM_unc::operator[](int l){

   return *dpphm[l];

}

/**
 * acces to the individual xTPM_unc objects: the const version
 * @param l the specific xTPM_unc object you want
 * @return the xTPM_unc object with parameter l
 */
const xTPM_unc &dPPHM_unc::operator[](int l) const{

   return *dpphm[l];

}

ostream &operator<<(ostream &output,const dPPHM_unc &dpphm_p){

   for(int l = 0;l < dpphm_p.gM();++l){

      output << std::endl;
      output << l << std::endl;
      output << std::endl;

      output << dpphm_p[l] << endl;

   }

   return output;

}

/**
 * @return the trace of the dPPHM_unc object
 */
double dPPHM_unc::trace() const{

   double ward = 0.0;

   for(int l = 0;l < M;++l)
      ward += dpphm[l]->trace();

   return ward;

}

/**
 * Scale the dPPHM_unc with parameter alpha
 * @param alpha scalefactor
 */
void dPPHM_unc::dscal(double alpha){

   for(int l = 0;l < M;++l)
      dpphm[l]->dscal(alpha);

}

/**
 * overload the equality operator
 * @param dpphm_c The dPPHM_unc you want to be copied into this
 */
dPPHM_unc &dPPHM_unc::operator=(const dPPHM_unc &dpphm_c){

   for(int l = 0;l < M;++l)
      *dpphm[l] = dpphm_c[l];

   return *this;

}

/**
 * Make all the number in your dPPHM_unc equal to the number a, e.g. usefull for initialization (dPPHM_unc W = 0)
 * @param a the number
 */
dPPHM_unc &dPPHM_unc::operator=(double a){

   for(int l = 0;l < M;++l)
      *dpphm[l] = a;

   return *this;

}

/**
 * overload the += operator for dPPHM_unc's
 * @param dpphm_pl The dPPHM_unc you want to add to this
 */
dPPHM_unc &dPPHM_unc::operator+=(const dPPHM_unc &dpphm_pl){

   for(int l = 0;l < M;++l)
      *dpphm[l] += dpphm_pl[l];

   return *this;

}

/**
 * overload the -= operator for dPPHM_unc's
 * @param dpphm_m The dPPHM_unc you want to deduct from this
 */
dPPHM_unc &dPPHM_unc::operator-=(const dPPHM_unc &dpphm_m){

   for(int l = 0;l < M;++l)
      *dpphm[l] -= dpphm_m[l];

   return *this;

}

/**
 * add the dpphm dpphm_pl times the constant alpha to *this
 * @param alpha the constant to multiply the dpphm_pl with
 * @param dpphm_pl the dPPHM_unc to be multiplied by alpha and added to this
 */
dPPHM_unc &dPPHM_unc::daxpy(double alpha,const dPPHM_unc &dpphm_pl){

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
void dPPHM_unc::L_map(const dPPHM_unc &map,const dPPHM_unc &object){

   for(int l = 0;l < M;++l)
      dpphm[l]->L_map(map[l],object[l]);

}

/**
 * dPPHM_unc product of two general matrices A en B, put result in this
 * @param A left dpphm
 * @param B right dpphm
 */
dPPHM_unc &dPPHM_unc::mprod(const dPPHM_unc &A,const dPPHM_unc &B){

   for(int l = 0;l < M;++l)
      dpphm[l]->mprod(A[l],B[l]);

   return *this;

}

/**
 * Take the square root out of the positive semidefinite dpphm, destroys original dpphm, square root will be put in (*this)
 * @param option = 1, positive square root, = -1, negative square root.
 */
void dPPHM_unc::sqrt(int option){

   for(int l = 0;l < M;++l)
      dpphm[l]->sqrt(option);

}

/**
 * @return inproduct of (*this) dpphm with dpphm_i, defined as Tr (A B)
 * @param dpphm_i input dpphm
 */
double dPPHM_unc::ddot(const dPPHM_unc &dpphm_i) const{

   double ward = 0.0;

   for(int l = 0;l < M;++l)
      ward += dpphm[l]->ddot(dpphm_i[l]);

   return ward;

}

/**
 * Invert positive semidefinite symmetric dpphm is stored in (*this), original dpphm (*this) is destroyed
 */
void dPPHM_unc::invert(){

   for(int l = 0;l < M;++l)
      dpphm[l]->invert();

}

/**
 * copy upper in lower part of dPPHM_unc object
 */
void dPPHM_unc::symmetrize(){

   for(int l = 0;l < M;++l)
      dpphm[l]->symmetrize();

}

/**
 * Fill the matrix with random numbers.
 */
void dPPHM_unc::fill_Random(){

   for(int l = 0;l < M;++l)
      dpphm[l]->fill_Random();

}

/**
 * print the eigenvalues, dirty fix for the fact that my Vector class isn't a Template class here.
 */
void dPPHM_unc::print_eig(){

   for(int l = 0;l < M;++l){

      Vector v(*dpphm[l]);

      cout << endl;
      cout << l << endl;
      cout << endl;

      cout << v << endl;

   }

}

/**
 * map a coupled dPPHM onto an uncoupled one, using CG coefficients (actually 3j's)
 * @param dpphm_c the coupled dPPHM
 */
void dPPHM_unc::uncouple(const dPPHM &dpphm_c){

   int a,b,c,d;
   int s_a,s_b,s_c,s_d;
   int s_l,s_l_;

   int sign_ab,sign_cd,sign_l;

   for(int l = 0;l < M;++l){

      for(int i = 0;i < dpphm[l]->gn();++i){

         s_l = xTPM_unc::gt2s(i,0);

         a = xTPM_unc::gt2s(i,1)/2;
         b = xTPM_unc::gt2s(i,2)/2;

         s_a = xTPM_unc::gt2s(i,1)%2;
         s_b = xTPM_unc::gt2s(i,2)%2;

         for(int j = i;j < dpphm[l]->gn();++j){

            s_l_ = xTPM_unc::gt2s(j,0);

            c = xTPM_unc::gt2s(j,1)/2;
            d = xTPM_unc::gt2s(j,2)/2;

            s_c = xTPM_unc::gt2s(j,1)%2;
            s_d = xTPM_unc::gt2s(j,2)%2;

            if(s_l == s_l_)
               sign_l = 1;
            else
               sign_l = -1;

            (*this)[l](i,j) = 0.0;

            for(int S = 1;S <= 3;S+=2)
               for(int M_ = -S;M_ <= S;M_+=2)
                  for(int S_ab = 0;S_ab <= 2;S_ab+=2)
                     for(int M_ab = -S_ab;M_ab <= S_ab;M_ab+=2)
                        for(int S_cd = 0;S_cd <= 2;S_cd+=2)
                           for(int M_cd = -S_cd;M_cd <= S_cd;M_cd+=2){

                              if(M_ab == 0)
                                 sign_ab = 1;
                              else
                                 sign_ab = -1;

                              if(M_cd == 0)
                                 sign_cd = 1;
                              else
                                 sign_cd = -1;

                              (*this)[l](i,j) += sign_l * (1 - S_ab) * (1 - S_cd) * sign_ab * sign_cd * (S + 1.0) * std::sqrt( (S_ab + 1.0) * (S_cd + 1.0) )

                                 * Tools::g3j(1,1,S_ab,2*s_a - 1,2*s_b - 1,-M_ab) * Tools::g3j(S_ab,1,S,M_ab,-(2*s_l - 1),-M_)

                                 * Tools::g3j(1,1,S_cd,2*s_c - 1,2*s_d - 1,-M_cd) * Tools::g3j(S_cd,1,S,M_cd,-(2*s_l_ - 1),-M_) * dpphm_c(l,S/2,S_ab/2,a,b,S_cd/2,c,d);

                           }

            if(a == b)
               (*this)[l](i,j) *= std::sqrt(2.0);

            if(c == d)
               (*this)[l](i,j) *= std::sqrt(2.0);

         }

      }

   }

   this->symmetrize();

}

/**
 * output using sp indices, for input in the regular v2.5DM program
 */
void dPPHM_unc::out_sp(const char *filename) const{

   ofstream out(filename);
   out.precision(15);

   int a,b,c,d;
   int s_l,s_l_;

   for(int l = 0;l < M;++l){

      for(int i = 0;i < dpphm[l]->gn();++i){

         s_l = xTPM_unc::gt2s(i,0);

         a = xTPM_unc::gt2s(i,1);
         b = xTPM_unc::gt2s(i,2);

         for(int j = i;j < dpphm[l]->gn();++j){

            s_l_ = xTPM_unc::gt2s(j,0);

            c = xTPM_unc::gt2s(j,1);
            d = xTPM_unc::gt2s(j,2);

            if(s_l == s_l_)
               out << 2*l + s_l << "\t" << a << "\t" << b << "\t" << c << "\t" << d << "\t" << (*this)[l](i,j) << endl;

         }
      }
   }

}

/**
 * I2 uncoupled
 */
void dPPHM_unc::I(const dDPM_unc &ddpm) {

   *this = 0;

   int a,b,c,d;
   int s_l,s_l_;

   for(int l = 0;l < M;++l){

      for(int i = 0;i < dpphm[l]->gn();++i){

         s_l = xTPM_unc::gt2s(i,0);

         a = xTPM_unc::gt2s(i,1);
         b = xTPM_unc::gt2s(i,2);

         for(int j = i;j < dpphm[l]->gn();++j){

            s_l_ = xTPM_unc::gt2s(j,0);

            c = xTPM_unc::gt2s(j,1);
            d = xTPM_unc::gt2s(j,2);

            (*this)[l](i,j) = -ddpm[l](s_l_,a,b,s_l,c,d);

            if(s_l == s_l_)
               for(int k = 0;k < M;++k)
                  for(int s_k = 0;s_k < 2;++s_k)
                     (*this)[l](i,j) += 1.0/(N - 2.0) * ddpm[k](s_k,a,b,s_k,c,d);

         }
      }
   }

   this->symmetrize();

}
