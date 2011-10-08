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

int dPHHM_unc::M;
int dPHHM_unc::N;

/**
 * initialize the static variables
 * @param M_in dimension of sp space
 * @param N_in nr of particles
 */
void dPHHM_unc::init(int M_in,int N_in){

   M = M_in;
   N = N_in;

}

/**
 * standard constructor: constructs M rxPHM_unc object with parameter l = 0 -> M-1
 */
dPHHM_unc::dPHHM_unc() {

   dphhm = new rxPHM_unc * [M];

   for(int l = 0;l < M;++l)
      dphhm[l] = new rxPHM_unc(l);

}

/**
 * copy constructor: constructs an array and copies the content of dphhm_c in it.
 * @param dphhm_c object that will be copied into this.
 */
dPHHM_unc::dPHHM_unc(const dPHHM_unc &dphhm_c) { 

   dphhm = new rxPHM_unc * [M];

   for(int l = 0;l < M;++l)
      dphhm[l] = new rxPHM_unc(dphhm_c[l]);
   
}

/**
 * destructor
 */
dPHHM_unc::~dPHHM_unc(){

   for(int l = 0;l < M;++l)
      delete dphhm[l];

   delete [] dphhm;
   
}

/**
 * @return nr of particles
 */
int dPHHM_unc::gN() const {

   return N;

}

/**
 * @return dimension of sp space
 */
int dPHHM_unc::gM() const {

   return M;

}

/**
 * acces to the individual rxPHM_unc objects
 * @param l the specific rxPHM_unc object you want
 * @return the rxPHM_unc object with parameter l
 */
rxPHM_unc &dPHHM_unc::operator[](int l){

   return *dphhm[l];

}

/**
 * acces to the individual rxPHM_unc objects: the const version
 * @param l the specific rxPHM_unc object you want
 * @return the rxPHM_unc object with parameter l
 */
const rxPHM_unc &dPHHM_unc::operator[](int l) const{

   return *dphhm[l];

}

ostream &operator<<(ostream &output,const dPHHM_unc &dphhm_p){

   for(int l = 0;l < dphhm_p.gM();++l){

      output << std::endl;
      output << l << std::endl;
      output << std::endl;

      output << dphhm_p[l] << endl;

   }

   return output;

}

/**
 * @return the trace of the dPHHM_unc object
 */
double dPHHM_unc::trace() const{

   double ward = 0.0;

   for(int l = 0;l < M;++l)
      ward += dphhm[l]->trace();

   return ward;

}

/**
 * Scale the dPHHM_unc with parameter alpha
 * @param alpha scalefactor
 */
void dPHHM_unc::dscal(double alpha){

   for(int l = 0;l < M;++l)
      dphhm[l]->dscal(alpha);

}

/**
 * overload the equality operator
 * @param dphhm_c The dPHHM_unc you want to be copied into this
 */
dPHHM_unc &dPHHM_unc::operator=(const dPHHM_unc &dphhm_c){

   for(int l = 0;l < M;++l)
      *dphhm[l] = dphhm_c[l];

   return *this;

}

/**
 * Make all the number in your dPHHM_unc equal to the number a, e.g. usefull for initialization (dPHHM_unc W = 0)
 * @param a the number
 */
dPHHM_unc &dPHHM_unc::operator=(double a){

   for(int l = 0;l < M;++l)
      *dphhm[l] = a;

   return *this;

}

/**
 * overload the += operator for dPHHM_unc's
 * @param dphhm_pl The dPHHM_unc you want to add to this
 */
dPHHM_unc &dPHHM_unc::operator+=(const dPHHM_unc &dphhm_pl){

   for(int l = 0;l < M;++l)
      *dphhm[l] += dphhm_pl[l];

   return *this;

}

/**
 * overload the -= operator for dPHHM_unc's
 * @param dphhm_m The dPHHM_unc you want to deduct from this
 */
dPHHM_unc &dPHHM_unc::operator-=(const dPHHM_unc &dphhm_m){

   for(int l = 0;l < M;++l)
      *dphhm[l] -= dphhm_m[l];

   return *this;

}

/**
 * add the dphhm dphhm_pl times the constant alpha to *this
 * @param alpha the constant to multiply the dphhm_pl with
 * @param dphhm_pl the dPHHM_unc to be multiplied by alpha and added to this
 */
dPHHM_unc &dPHHM_unc::daxpy(double alpha,const dPHHM_unc &dphhm_pl){

   for(int l = 0;l < M;++l)
      dphhm[l]->daxpy(alpha,dphhm_pl[l]);

   return *this;

}

/**
 * Multiply symmetric dphhm object left en right with symmetric dphhm map to 
 * form another symmetric dphhm and put it in (*this): this = map*object*map
 * @param map dphhm that will be multiplied to the left en to the right of dphhm object
 * @param object central dphhm
 */
void dPHHM_unc::L_map(const dPHHM_unc &map,const dPHHM_unc &object){

   for(int l = 0;l < M;++l)
      dphhm[l]->L_map(map[l],object[l]);

}

/**
 * dPHHM_unc product of two general matrices A en B, put result in this
 * @param A left dphhm
 * @param B right dphhm
 */
dPHHM_unc &dPHHM_unc::mprod(const dPHHM_unc &A,const dPHHM_unc &B){

   for(int l = 0;l < M;++l)
      dphhm[l]->mprod(A[l],B[l]);

   return *this;

}

/**
 * Take the square root out of the positive semidefinite dphhm, destroys original dphhm, square root will be put in (*this)
 * @param option = 1, positive square root, = -1, negative square root.
 */
void dPHHM_unc::sqrt(int option){

   for(int l = 0;l < M;++l)
      dphhm[l]->sqrt(option);

}

/**
 * @return inproduct of (*this) dphhm with dphhm_i, defined as Tr (A B)
 * @param dphhm_i input dphhm
 */
double dPHHM_unc::ddot(const dPHHM_unc &dphhm_i) const{

   double ward = 0.0;

   for(int l = 0;l < M;++l)
      ward += dphhm[l]->ddot(dphhm_i[l]);

   return ward;

}

/**
 * Invert positive semidefinite symmetric dphhm is stored in (*this), original dphhm (*this) is destroyed
 */
void dPHHM_unc::invert(){

   for(int l = 0;l < M;++l)
      dphhm[l]->invert();

}

/**
 * copy upper in lower part of dPHHM_unc object
 */
void dPHHM_unc::symmetrize(){

   for(int l = 0;l < M;++l)
      dphhm[l]->symmetrize();

}

/**
 * Fill the matrix with random numbers.
 */
void dPHHM_unc::fill_Random(){

   for(int l = 0;l < M;++l)
      dphhm[l]->fill_Random();

}

/**
 * print the eigenvalues, dirty fix for the fact that my Vector class isn't a Template class here.
 */
void dPHHM_unc::print_eig(){

   for(int l = 0;l < M;++l){

      Vector v(*dphhm[l]);

      cout << endl;
      cout << l << endl;
      cout << endl;

      cout << v << endl;

   }

}

/**
 * map a coupled dPHHM onto an uncoupled one, using CG coefficients (actually 3j's)
 * @param dphhm_c the coupled dPHHM
 */
void dPHHM_unc::uncouple(const dPHHM &dphhm_c){

   int a,b,c,d;
   int s_a,s_b,s_c,s_d;
   int s_l,s_l_;

   int sign_bl,sign_dl,sign_ac;

   for(int l = 0;l < M;++l){

      for(int i = 0;i < dphhm[l]->gn();++i){

         s_l = rxPHM_unc::gph2s(l,i,0);

         a = rxPHM_unc::gph2s(l,i,1)/2;
         b = rxPHM_unc::gph2s(l,i,2)/2;

         s_a = rxPHM_unc::gph2s(l,i,1)%2;
         s_b = rxPHM_unc::gph2s(l,i,2)%2;

         for(int j = i;j < dphhm[l]->gn();++j){

            s_l_ = rxPHM_unc::gph2s(l,j,0);

            c = rxPHM_unc::gph2s(l,j,1)/2;
            d = rxPHM_unc::gph2s(l,j,2)/2;

            s_c = rxPHM_unc::gph2s(l,j,1)%2;
            s_d = rxPHM_unc::gph2s(l,j,2)%2;

            if(s_a == s_c)
               sign_ac = 1;
            else
               sign_ac = -1;

            (*this)[l](i,j) = 0.0;

            for(int S = 1;S <= 3;S+=2)
               for(int M_ = -S;M_ <= S;M_+=2)
                  for(int S_bl = 0;S_bl <= 2;S_bl+=2)
                     for(int M_bl = -S_bl;M_bl <= S_bl;M_bl+=2)
                        for(int S_dl = 0;S_dl <= 2;S_dl+=2)
                           for(int M_dl = -S_dl;M_dl <= S_dl;M_dl+=2){

                              if(M_bl == 0)
                                 sign_bl = 1;
                              else
                                 sign_bl = -1;

                              if(M_dl == 0)
                                 sign_dl = 1;
                              else
                                 sign_dl = -1;

                              (*this)[l](i,j) += sign_ac * (1 - S_bl) * (1 - S_dl) * sign_bl * sign_dl * (S + 1.0) * std::sqrt( (S_bl + 1.0) * (S_dl + 1.0) )

                                 * Tools::g3j(1,1,S_bl,2*s_b - 1,2*s_l - 1,-M_bl) * Tools::g3j(1,S_bl,S,-(2*s_a - 1),M_bl,-M_)

                                 * Tools::g3j(1,1,S_dl,2*s_d - 1,2*s_l_ - 1,-M_dl) * Tools::g3j(1,S_dl,S,-(2*s_c - 1),M_dl,-M_) * dphhm_c(l,S/2,S_bl/2,a,b,S_dl/2,c,d);

                           }

         }

      }

   }

   this->symmetrize();

}

/**
 * output using sp indices, for input in the regular v2.5DM program
 */
void dPHHM_unc::out_sp(const char *filename) const{

   ofstream out(filename);
   out.precision(15);

   int a,b,c,d;
   int s_l,s_l_;

   for(int l = 0;l < M;++l){

      for(int i = 0;i < dphhm[l]->gn();++i){

         s_l = rxPHM_unc::gph2s(l,i,0);

         a = rxPHM_unc::gph2s(l,i,1);
         b = rxPHM_unc::gph2s(l,i,2);

         for(int j = i;j < dphhm[l]->gn();++j){

            s_l_ = rxPHM_unc::gph2s(l,j,0);

            c = rxPHM_unc::gph2s(l,j,1);
            d = rxPHM_unc::gph2s(l,j,2);

            if(s_l == s_l_)
               out << 2*l + s_l << "\t" << a << "\t" << b << "\t" << c << "\t" << d << "\t" << (*this)[l](i,j) << endl;

         }
      }
   }

}
