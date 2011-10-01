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

int dDPM_unc::M;
int dDPM_unc::N;

/**
 * initialize the static variables
 * @param M_in dimension of sp space
 * @param N_in nr of particles
 */
void dDPM_unc::init(int M_in,int N_in){

   M = M_in;
   N = N_in;

}

/**
 * standard constructor: constructs M rxTPM_unc object with parameter l = 0 -> M-1
 */
dDPM_unc::dDPM_unc() {

   ddpm = new rxTPM_unc * [M];

   for(int l = 0;l < M;++l)
      ddpm[l] = new rxTPM_unc(l);

}

/**
 * copy constructor: constructs an array and copies the content of ddpm_c in it.
 * @param ddpm_c object that will be copied into this.
 */
dDPM_unc::dDPM_unc(const dDPM_unc &ddpm_c) { 

   ddpm = new rxTPM_unc * [M];

   for(int l = 0;l < M;++l)
      ddpm[l] = new rxTPM_unc(ddpm_c[l]);
   
}

/**
 * destructor
 */
dDPM_unc::~dDPM_unc(){

   for(int l = 0;l < M;++l)
      delete ddpm[l];

   delete [] ddpm;
   
}

/**
 * @return nr of particles
 */
int dDPM_unc::gN() const {

   return N;

}

/**
 * @return dimension of sp space
 */
int dDPM_unc::gM() const {

   return M;

}

/**
 * acces to the individual rxTPM_unc objects
 * @param l the specific rxTPM_unc object you want
 * @return the rxTPM_unc object with parameter l
 */
rxTPM_unc &dDPM_unc::operator[](int l){

   return *ddpm[l];

}

/**
 * acces to the individual rxTPM_unc objects: the const version
 * @param l the specific rxTPM_unc object you want
 * @return the rxTPM_unc object with parameter l
 */
const rxTPM_unc &dDPM_unc::operator[](int l) const{

   return *ddpm[l];

}

ostream &operator<<(ostream &output,const dDPM_unc &ddpm_p){

   for(int l = 0;l < ddpm_p.gM();++l){

      output << std::endl;
      output << l << std::endl;
      output << std::endl;

      output << ddpm_p[l] << endl;

   }

   return output;

}

/**
 * @return the trace of the dDPM_unc object
 */
double dDPM_unc::trace() const{

   double ward = 0.0;

   for(int l = 0;l < M;++l)
      ward += ddpm[l]->trace();

   return ward;

}

/**
 * Scale the dDPM_unc with parameter alpha
 * @param alpha scalefactor
 */
void dDPM_unc::dscal(double alpha){

   for(int l = 0;l < M;++l)
      ddpm[l]->dscal(alpha);

}

/**
 * overload the equality operator
 * @param ddpm_c The dDPM_unc you want to be copied into this
 */
dDPM_unc &dDPM_unc::operator=(const dDPM_unc &ddpm_c){

   for(int l = 0;l < M;++l)
      *ddpm[l] = ddpm_c[l];

   return *this;

}

/**
 * Make all the number in your dDPM_unc equal to the number a, e.g. usefull for initialization (dDPM_unc W = 0)
 * @param a the number
 */
dDPM_unc &dDPM_unc::operator=(double a){

   for(int l = 0;l < M;++l)
      *ddpm[l] = a;

   return *this;

}

/**
 * overload the += operator for dDPM_unc's
 * @param ddpm_pl The dDPM_unc you want to add to this
 */
dDPM_unc &dDPM_unc::operator+=(const dDPM_unc &ddpm_pl){

   for(int l = 0;l < M;++l)
      *ddpm[l] += ddpm_pl[l];

   return *this;

}

/**
 * overload the -= operator for dDPM_unc's
 * @param ddpm_m The dDPM_unc you want to deduct from this
 */
dDPM_unc &dDPM_unc::operator-=(const dDPM_unc &ddpm_m){

   for(int l = 0;l < M;++l)
      *ddpm[l] -= ddpm_m[l];

   return *this;

}

/**
 * add the ddpm ddpm_pl times the constant alpha to *this
 * @param alpha the constant to multiply the ddpm_pl with
 * @param ddpm_pl the dDPM_unc to be multiplied by alpha and added to this
 */
dDPM_unc &dDPM_unc::daxpy(double alpha,const dDPM_unc &ddpm_pl){

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
void dDPM_unc::L_map(const dDPM_unc &map,const dDPM_unc &object){

   for(int l = 0;l < M;++l)
      ddpm[l]->L_map(map[l],object[l]);

}

/**
 * dDPM_unc product of two general matrices A en B, put result in this
 * @param A left ddpm
 * @param B right ddpm
 */
dDPM_unc &dDPM_unc::mprod(const dDPM_unc &A,const dDPM_unc &B){

   for(int l = 0;l < M;++l)
      ddpm[l]->mprod(A[l],B[l]);

   return *this;

}

/**
 * Take the square root out of the positive semidefinite ddpm, destroys original ddpm, square root will be put in (*this)
 * @param option = 1, positive square root, = -1, negative square root.
 */
void dDPM_unc::sqrt(int option){

   for(int l = 0;l < M;++l)
      ddpm[l]->sqrt(option);

}

/**
 * @return inproduct of (*this) ddpm with ddpm_i, defined as Tr (A B)
 * @param ddpm_i input ddpm
 */
double dDPM_unc::ddot(const dDPM_unc &ddpm_i) const{

   double ward = 0.0;

   for(int l = 0;l < M;++l)
      ward += ddpm[l]->ddot(ddpm_i[l]);

   return ward;

}

/**
 * Invert positive semidefinite symmetric ddpm is stored in (*this), original ddpm (*this) is destroyed
 */
void dDPM_unc::invert(){

   for(int l = 0;l < M;++l)
      ddpm[l]->invert();

}

/**
 * copy upper in lower part of dDPM_unc object
 */
void dDPM_unc::symmetrize(){

   for(int l = 0;l < M;++l)
      ddpm[l]->symmetrize();

}

/**
 * Fill the matrix with random numbers.
 */
void dDPM_unc::fill_Random(){

   for(int l = 0;l < M;++l)
      ddpm[l]->fill_Random();

}

/**
 * print the eigenvalues, dirty fix for the fact that my Vector class isn't a Template class here.
 */
void dDPM_unc::print_eig(){

   for(int l = 0;l < M;++l){

      Vector v(*ddpm[l]);

      cout << endl;
      cout << l << endl;
      cout << endl;

      cout << v << endl;

   }

}

/**
 * map a coupled dDPM onto an uncoupled one, using CG coefficients (actually 3j's)
 * @param ddpm_c the coupled dDPM
 */
void dDPM_unc::uncouple(const dDPM &ddpm_c){

   int a,b,c,d;
   int s_a,s_b,s_c,s_d;
   int s_l,s_l_;

   int sign_ab,sign_cd;

   for(int l = 0;l < M;++l){

      for(int i = 0;i < ddpm[l]->gn();++i){

         s_l = rxTPM_unc::gt2s(l,i,0);

         a = rxTPM_unc::gt2s(l,i,1)/2;
         b = rxTPM_unc::gt2s(l,i,2)/2;

         s_a = rxTPM_unc::gt2s(l,i,1)%2;
         s_b = rxTPM_unc::gt2s(l,i,2)%2;

         for(int j = i;j < ddpm[l]->gn();++j){

            s_l_ = rxTPM_unc::gt2s(l,j,0);

            c = rxTPM_unc::gt2s(l,j,1)/2;
            d = rxTPM_unc::gt2s(l,j,2)/2;

            s_c = rxTPM_unc::gt2s(l,j,1)%2;
            s_d = rxTPM_unc::gt2s(l,j,2)%2;

            (*this)[l](i,j) = 0.0;

            for(int S = 1;S <= 3;S+=2)
               for(int M = -S;M <= S;M+=2)
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

                              (*this)[l](i,j) += (1 - S_ab) * (1 - S_cd) * sign_ab * sign_cd * (S + 1.0) * std::sqrt( (S_ab + 1.0) * (S_cd + 1.0) )

                              * Tools::g3j(1,1,S_ab,2*s_a - 1,2*s_b - 1,-M_ab) * Tools::g3j(S_ab,1,S,M_ab,2*s_l - 1,-M)

                              * Tools::g3j(1,1,S_cd,2*s_c - 1,2*s_d - 1,-M_cd) * Tools::g3j(S_cd,1,S,M_cd,2*s_l_ - 1,-M) * ddpm_c(l,S/2,S_ab/2,a,b,S_cd/2,c,d);

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
