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

int dPHHM::M;
int dPHHM::N;

/**
 * initialize the static variables
 * @param M_in dimension of spatial sp space
 * @param N_in nr of particles
 */
void dPHHM::init(int M_in,int N_in){

   M = M_in;
   N = N_in;

}

/**
 * standard constructor: constructs M rxPHM object with parameter l = 0 -> M-1
 */
dPHHM::dPHHM() {

   dphhm = new rxPHM * [M];

   for(int l = 0;l < M;++l)
      dphhm[l] = new rxPHM(l);

}

/**
 * copy constructor: constructs an array and copies the content of W in it.
 * @param W object that will be copied into this.
 */
dPHHM::dPHHM(const dPHHM &dphhm_c) { 

   dphhm = new rxPHM * [M];

   for(int l = 0;l < M;++l)
      dphhm[l] = new rxPHM(dphhm_c[l]);

}

/**
 * destructor
 */
dPHHM::~dPHHM(){

   for(int l = 0;l < M;++l)
      delete dphhm[l];

   delete [] dphhm;

}

/**
 * @return nr of particles
 */
int dPHHM::gN() const {

   return N;

}

/**
 * @return dimension of sp space
 */
int dPHHM::gM() const {

   return M;

}

/**
 * acces to the individual rxPHM objects
 * @param l the specific rxPHM object you want
 * @return the rxPHM object with parameter l
 */
rxPHM &dPHHM::operator[](int l){

   return *dphhm[l];

}

/**
 * acces to the individual rxPHM objects: the const version
 * @param l the specific rxPHM object you want
 * @return the rxPHM object with parameter l
 */
const rxPHM &dPHHM::operator[](int l) const{

   return *dphhm[l];

}

/**
 * access the numbers in sp mode
 * @param l blockindex
 * @param S pph spin
 * @param S_bl intermediate spin for b and l
 * @param a first spatial index for the row
 * @param b second spatial index for the row
 * @param S_dl intermediate spin for d and l
 * @param c first spatial index for the column
 * @param d second spatial index for the column
 */
double dPHHM::operator()(int l,int S,int S_bl,int a,int b,int S_dl,int c,int d) const {

   if( S == 1 && ( S_bl == 0 || S_dl == 0) )
      return 0.0;

   if(S_bl == 1 && b == l)
      return 0.0;

   if(S_dl == 1 && d == l)
      return 0.0;

   return (*this)[l](S,rxPHM::gs2ph(l,S,S_bl,a,b),rxPHM::gs2ph(l,S,S_dl,c,d));

}

ostream &operator<<(ostream &output,const dPHHM &dphhm_p){

   for(int l = 0;l < dphhm_p.gM();++l){

      output << std::endl;
      output << "l = \t" << l << std::endl;
      output << std::endl;

      output << dphhm_p[l] << endl;

   }

   return output;

}

/**
 * @return the trace of the dPHHM object
 */
double dPHHM::trace() const{

   double ward = 0.0;

   for(int l = 0;l < M;++l)
      ward += dphhm[l]->trace();

   return ward;

}

/**
 * Scale the dPHHM with parameter alpha
 * @param alpha scalefactor
 */
void dPHHM::dscal(double alpha){

   for(int l = 0;l < M;++l)
      dphhm[l]->dscal(alpha);

}

/**
 * overload the equality operator
 * @param dphhm_c The dPHHM you want to be copied into this
 */
dPHHM &dPHHM::operator=(const dPHHM &dphhm_c){

   for(int l = 0;l < M;++l)
      *dphhm[l] = dphhm_c[l];

   return *this;

}

/**
 * Make all the number in your dPHHM equal to the number a, e.g. usefull for initialization (dPHHM W = 0)
 * @param a the number
 */
dPHHM &dPHHM::operator=(double a){

   for(int l = 0;l < M;++l)
      *dphhm[l] = a;

   return *this;

}

/**
 * overload the += operator for dPHHM's
 * @param dphhm_pl The dPHHM you want to add to this
 */
dPHHM &dPHHM::operator+=(const dPHHM &dphhm_pl){

   for(int l = 0;l < M;++l)
      *dphhm[l] += dphhm_pl[l];

   return *this;

}

/**
 * overload the -= operator for dPHHM's
 * @param dphhm_m The dPHHM you want to deduct from this
 */
dPHHM &dPHHM::operator-=(const dPHHM &dphhm_m){

   for(int l = 0;l < M;++l)
      *dphhm[l] -= dphhm_m[l];

   return *this;

}

/**
 * add the dphhm dphhm_pl times the constant alpha to *this
 * @param alpha the constant to multiply the dphhm_pl with
 * @param dphhm_pl the dPHHM to be multiplied by alpha and added to this
 */
dPHHM &dPHHM::daxpy(double alpha,const dPHHM &dphhm_pl){

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
void dPHHM::L_map(const dPHHM &map,const dPHHM &object){

   for(int l = 0;l < M;++l)
      dphhm[l]->L_map(map[l],object[l]);

}

/**
 * dPHHM product of two general matrices A en B, put result in this
 * @param A left dphhm
 * @param B right dphhm
 */
dPHHM &dPHHM::mprod(const dPHHM &A,const dPHHM &B){

   for(int l = 0;l < M;++l)
      dphhm[l]->mprod(A[l],B[l]);

   return *this;

}

/**
 * Take the square root out of the positive semidefinite dphhm, destroys original dphhm, square root will be put in (*this)
 * @param option = 1, positive square root, = -1, negative square root.
 */
void dPHHM::sqrt(int option){

   for(int l = 0;l < M;++l)
      dphhm[l]->sqrt(option);

}

/**
 * @return inproduct of (*this) dphhm with dphhm_i, defined as Tr (A B)
 * @param dphhm_i input dphhm
 */
double dPHHM::ddot(const dPHHM &dphhm_i) const{

   double ward = 0.0;

   for(int l = 0;l < M;++l)
      ward += dphhm[l]->ddot(dphhm_i[l]);

   return ward;

}

/**
 * Invert positive semidefinite symmetric dphhm is stored in (*this), original dphhm (*this) is destroyed
 */
void dPHHM::invert(){

   for(int l = 0;l < M;++l)
      dphhm[l]->invert();

}

/**
 * copy upper in lower part of dPHHM object
 */
void dPHHM::symmetrize(){

   for(int l = 0;l < M;++l)
      dphhm[l]->symmetrize();

}

/**
 * Fill the matrix with random numbers.
 */
void dPHHM::fill_Random(){

   for(int l = 0;l < M;++l)
      dphhm[l]->fill_Random();

}

/**
 * Map a dDPM on a dPHHM using the G1 map
 * @param ddpm input dDPM
 */
void dPHHM::G1(const dDPM &ddpm){

   int a,b,c,d;
   int S_bl,S_dl;

   int sign_bl,sign_dl;

   TPM tpm;
   tpm.bar(1.0/(N - 2.0),ddpm);

   SPM spm;
   spm.bar(1.0/(N - 1.0),tpm);

   for(int l = 0;l < M;++l){

      //first S = 1/2
      for(int i = 0;i < dphhm[l]->gdim(0);++i){

         S_bl = rxPHM::gph2s(l,0,i,0);

         a = rxPHM::gph2s(l,0,i,1);
         b = rxPHM::gph2s(l,0,i,2);

         sign_bl = 1 - 2*S_bl;

         for(int j = 0;j < dphhm[l]->gdim(0);++j){

            S_dl = rxPHM::gph2s(l,0,j,0);

            c = rxPHM::gph2s(l,0,j,1);
            d = rxPHM::gph2s(l,0,j,2);

            sign_dl = 1 - 2*S_dl;

            //let the games begin: first the dp
            (*this)[l](0,i,j) = 0.0;
/*
            for(int S_ = 0;S_ < 2;S_++)
               for(int S_ab = 0;S_ab < 2;++S_ab)
                  for(int S_cd = 0;S_cd < 2;++S_cd){

                     (*this)[l](0,i,j) += (2*(S_ + 0.5) + 1.0) * std::sqrt( (2.0*S_ab + 1.0) * (2.0*S_cd + 1.0) * (2.0*S_bl + 1.0) * (2.0*S_dl + 1.0) )
                     
                        * Tools::g6j(2*S_ + 1,1,2*S_dl,1,1,2*S_ab) * Tools::g6j(2*S_ + 1,1,2*S_bl,1,1,2*S_cd) * Tools::g6j(2*S_ + 1,2*S_bl,1,1,2*S_dl,1)

                        * ddpm(l,S_,S_ab,a,d,S_cd,c,b);

                  }

            if(a == d)
               (*this)[l](0,i,j) *= std::sqrt(2.0);

            if(c == b)
               (*this)[l](0,i,j) *= std::sqrt(2.0);
*/
            //tp(1)
            if(S_bl == S_dl)
               if(a == c){

                  double hard = tpm(S_bl,b,l,d,l);

                  if(b == l)
                     hard *= std::sqrt(2.0);

                  if(d == l)
                     hard *= std::sqrt(2.0);

                  (*this)[l](0,i,j) += hard;

               }

            //part only in S = 1/2 blocks
            if(a == b){

               //sp_d
               if(c == d)
                  (*this)[l](0,i,j) += 0.5 * std::sqrt( (2.0*S_bl + 1.0) * (2.0*S_dl + 1.0) ) * sign_bl * sign_dl * spm(l,l);

               //sp_b
               if(c == l)
                  (*this)[l](0,i,j) += 0.5 * std::sqrt( (2.0*S_bl + 1.0) * (2.0*S_dl + 1.0) ) * sign_bl * spm(d,l);
/*
               //tp(2)_a
               double hard = tpm(S_dl,c,l,d,l);

               if(c == l)
                  hard *= std::sqrt(2.0);

               if(d == l)
                  hard *= std::sqrt(2.0);

               (*this)[l](0,i,j) -= 0.5 * std::sqrt( (2.0*S_bl + 1.0) * (2.0*S_dl + 1.0) ) * sign_bl * sign_dl * hard;
*/
            }

            if(c == d){

               //sp_c
               if(a == l)
                  (*this)[l](0,i,j) += 0.5 * std::sqrt( (2.0*S_bl + 1.0) * (2.0*S_dl + 1.0) ) * sign_dl * spm(b,l);
/*
               //tp(2)_d
               double hard = tpm(S_bl,a,l,b,l);

               if(a == l)
                  hard *= std::sqrt(2.0);

               if(b == l)
                  hard *= std::sqrt(2.0);

               (*this)[l](0,i,j) -= 0.5 * std::sqrt( (2.0*S_bl + 1.0) * (2.0*S_dl + 1.0) ) * sign_bl * sign_dl * hard;
*/
            }

            if(a == l){

               //sp_a
               if(c == l)
                  (*this)[l](0,i,j) += 0.5 * std::sqrt( (2.0*S_bl + 1.0) * (2.0*S_dl + 1.0) ) * spm(b,d);
/*
               //tp(2)_b
               double hard = tpm(S_dl,c,b,d,l);

               if(c == b)
                  hard *= std::sqrt(2.0);

               if(d == l)
                  hard *= std::sqrt(2.0);

               (*this)[l](0,i,j) -= 0.5 * std::sqrt( (2.0*S_bl + 1.0) * (2.0*S_dl + 1.0) ) * sign_dl * hard;
*/
            }
/*
            if(c == l){

               //tp(2)_c
               double hard = tpm(S_bl,a,d,b,l);

               if(a == d)
                  hard *= std::sqrt(2.0);

               if(b == l)
                  hard *= std::sqrt(2.0);

               (*this)[l](0,i,j) -= 0.5 * std::sqrt( (2.0*S_bl + 1.0) * (2.0*S_dl + 1.0) ) * sign_bl * hard;

            }
*/
         }
      }

      //then S = 3/2
      for(int i = 0;i < dphhm[l]->gdim(1);++i){

         S_bl = rxPHM::gph2s(l,1,i,0);

         a = rxPHM::gph2s(l,1,i,1);
         b = rxPHM::gph2s(l,1,i,2);

         for(int j = 0;j < dphhm[l]->gdim(1);++j){

            S_dl = rxPHM::gph2s(l,1,j,0);

            c = rxPHM::gph2s(l,1,j,1);
            d = rxPHM::gph2s(l,1,j,2);

            //dp
            (*this)[l](1,i,j) = 0.0;
/*
            for(int S_ = 0;S_ < 2;S_++)
               for(int S_ab = 0;S_ab < 2;++S_ab)
                  for(int S_cd = 0;S_cd < 2;++S_cd){

                     (*this)[l](1,i,j) += (2*(S_ + 0.5) + 1.0) * std::sqrt( (2.0*S_ab + 1.0) * (2.0*S_cd + 1.0) * (2.0*S_bl + 1.0) * (2.0*S_dl + 1.0) )
                     
                        * Tools::g6j(2*S_ + 1,1,2*S_dl,1,1,2*S_ab) * Tools::g6j(2*S_ + 1,1,2*S_bl,1,1,2*S_cd) * Tools::g6j(2*S_ + 1,2*S_bl,1,3,2*S_dl,1)

                        * ddpm(l,S_,S_ab,a,d,S_cd,c,b);

                  }

            if(a == d)
               (*this)[l](1,i,j) *= std::sqrt(2.0);

            if(c == b)
               (*this)[l](1,i,j) *= std::sqrt(2.0);
*/
            //tp(1)
            if(a == c)
               (*this)[l](1,i,j) += tpm(1,b,l,d,l);

         }
      }

   }

   //this->symmetrize();

}

/**
 * test the symmetry after the map
 */
void dPHHM::test_sym() const {

   for(int l = 0;l < M;++l){

      cout << "l = " << l << endl;
      cout << endl;

      for(int S = 0;S < 2;++S)
         for(int S_bl = 0;S_bl < 2;++S_bl)
            for(int S_dl = 0;S_dl < 2;++S_dl)
               for(int a = 0;a < M;++a)
                  for(int b = 0;b < M;++b)
                     for(int c = 0;c < M;++c){

                        if(fabs((*this)(l,S,S_bl,a,b,S_dl,c,b)- (1 - 2*S_bl) * (1 - 2*S_dl) * (*this)(b,S,S_bl,a,l,S_dl,c,l)) > 1.0e-12){

                           cout << 2*S + 1 << "/2\t(" << S_bl << ";" << S_dl << ")\t|\t" << a << "\t" << b << "\t" << c << "\t|\t"
                        
                           << (*this)(l,S,S_bl,a,b,S_dl,c,b) << "\t" << (1 - 2*S_bl) * (1 - 2*S_dl) * (*this)(b,S,S_bl,a,l,S_dl,c,l) << endl;

                        }

                     }

   }

}
