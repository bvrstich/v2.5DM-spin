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

double **dDPM::_6j;

/**
 * initialize the static variables
 * @param M_in dimension of spatial sp space
 * @param N_in nr of particles
 */
void dDPM::init(int M_in,int N_in){

   M = M_in;
   N = N_in;

   //allocate
   _6j = new double * [2];

   for(int S = 0;S < 2;++S)
      _6j[S] = new double [2];

   //initialize
   _6j[0][0] = -0.5;
   _6j[0][1] = 0.5;
   _6j[1][0] = 0.5;
   _6j[1][1] = 1.0/6.0;

}

/**
 * function that deallocates the static lists
 */
void dDPM::clear(){

   for(int S = 0;S < 2;++S)
      delete [] _6j[S];

   delete [] _6j;

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
double dDPM::operator()(int l,int S,int S_ab,int a,int b,int S_cd,int c,int d) const {

   int l_copy = l;

   double phasefac_i = get_inco(l_copy,S,S_ab,a,b);

   if(phasefac_i == 0)
      return 0.0;

   double phasefac_j = get_inco(l,S,S_cd,c,d);

   if(phasefac_j == 0)
      return 0.0;

   if(l_copy != l)
      return 0.0;

   return phasefac_i * phasefac_j * (*this)[l](S,rTPM::gs2t(l,S,S_ab,a,b),rTPM::gs2t(l,S,S_cd,c,d));

}

/**
 * reorder the sp indices and return the right phase and factor for the new vector
 * @param l third (block) index
 * @param S dp spin
 * @param S_ab intermediate spin
 * @param a first sp index
 * @param b second sp index
 */
double dDPM::get_inco(int &l,int S,int &S_ab,int &a,int &b){

   if(l == a && l == b)
      return 0.0;

   if(S == 0){//S == 1/2

      if(a == l){

         if(S_ab == 0){

            l = b;
            b = a;

            return -1.0/std::sqrt(2.0);

         }
         else{

            l = b;
            b = a;

            S_ab = 0;

            return -3.0/std::sqrt(2.0);

         }

      }
      else if(b == l){

         if(S_ab == 0){

            l = a;
            a = b;

            return -1.0/std::sqrt(2.0);

         }
         else{

            l = a;
            a = b;

            S_ab = 0;

            return 3.0/std::sqrt(2.0);

         }

      }
      else{//a != l and b != l

         if(a > b)
            return 1.0 - 2.0*S_ab;
         else
            return 1.0;

      }

   }
   else{//S == 3/2

      if(a == b || a == l || b == l)
         return 0.0;

      if(S_ab == 0)
         return 0.0;

      if(a > b)
         return -1.0;
      else
         return 1.0;

   }

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

/**
 * project this on the a dDPM object with correct symmetry with the third index
 */
void dDPM::proj_W(){

   int i,j;

   for(int l = 0;l < M;++l){

      for(int a = 0;a < M;++a){

         if(a == l)
            ++a;

         if(a == M)
            break;

         for(int b = a + 1;b < M;++b){

            if(b == l)
               ++b;

            if(b == M)
               break;

            //first set the 4 l block elements right!
            for(int S_ab = 0;S_ab < 2;++S_ab)
               for(int S_cd = S_ab;S_cd < 2;++S_cd){

                  //1
                  double ward = (*this)(l,0,S_ab,a,b,S_cd,a,b);

                  //2
                  for(int S_al = 0;S_al < 2;++S_al)
                     for(int S_cl = 0;S_cl < 2;++S_cl){
                         
                         ward += std::sqrt( (2.0*S_ab + 1.0) * (2.0*S_cd + 1.0) * (2.0*S_al + 1.0) * (2.0*S_cl + 1.0) ) * (1 - 2*S_ab) * (1 - 2*S_cd) * (1 - 2*S_al)

                            * (1 - 2*S_cl) * _6j[S_al][S_ab] * _6j[S_cl][S_cd] * (*this)(b,0,S_al,a,l,S_cl,a,l);

                     }

                  //3
                  for(int S_lb = 0;S_lb < 2;++S_lb)
                     for(int S_ld = 0;S_ld < 2;++S_ld){
                         
                         ward += std::sqrt( (2.0*S_ab + 1.0) * (2.0*S_cd + 1.0) * (2.0*S_lb + 1.0) * (2.0*S_ld + 1.0) )

                            * _6j[S_lb][S_ab] * _6j[S_ld][S_cd] * (*this)(b,0,S_lb,l,b,S_ld,l,b);

                     }

                  i = rTPM::gs2t(l,0,S_ab,a,b);
                  j = rTPM::gs2t(l,0,S_cd,a,b);

                  (*this)[l](0,i,j) = (*this)[l](0,j,i) = ward/3.0;

               }

               //then the 4 "b" block elements
               for(int S_al = 0;S_al < 2;++S_al)
                  for(int S_cl = S_al;S_cl < 2;++S_cl){

                     i = rTPM::gs2t(b,0,S_al,a,l);
                     j = rTPM::gs2t(b,0,S_cl,a,l);

                     (*this)[b](0,i,j) = 0.0;

                     for(int S_ab = 0;S_ab < 2;++S_ab)
                        for(int S_cd = 0;S_cd < 2;++S_cd){

                           (*this)[b](0,i,j) += std::sqrt( (2.0*S_ab + 1.0) * (2.0*S_cd + 1.0) * (2.0*S_al + 1.0) * (2.0*S_cl + 1.0) ) 

                              * (1 - 2*S_ab) * (1 - 2*S_cd) * (1 - 2*S_al) * (1 - 2*S_cl) * _6j[S_al][S_ab] * _6j[S_cl][S_cd] * (*this)(l,0,S_ab,a,b,S_cd,a,b);

                        }

                     (*this)[b](0,j,i) = (*this)[b](0,i,j);

                  }

               //finally the 4 "a" block elements
               for(int S_lb = 0;S_lb < 2;++S_lb)
                  for(int S_ld = S_lb;S_ld < 2;++S_ld){

                     i = rTPM::gs2t(a,0,S_lb,l,b);
                     j = rTPM::gs2t(a,0,S_ld,l,b);

                     (*this)[a](0,i,j) = 0.0;

                     for(int S_ab = 0;S_ab < 2;++S_ab)
                        for(int S_cd = 0;S_cd < 2;++S_cd){

                           (*this)[a](0,i,j) += std::sqrt( (2.0*S_ab + 1.0) * (2.0*S_cd + 1.0) * (2.0*S_lb + 1.0) * (2.0*S_ld + 1.0) ) 

                              * _6j[S_lb][S_ab] * _6j[S_ld][S_cd] * (*this)(l,0,S_ab,a,b,S_cd,a,b);

                        }

                     (*this)[a](0,j,i) = (*this)[a](0,i,j);

                  }


         }

      }

   }

}

/**
 * test if the projection is correct
 */
void dDPM::test_proj() const {

   for(int l = 0;l < M;++l){

      cout << l << endl;
      cout << endl;

      for(int a = 0;a < M;++a)
         for(int b = 0;b < M;++b){

            for(int S_ab = 0;S_ab < 2;++S_ab)
               for(int S_cd = 0;S_cd < 2;++S_cd){

                  cout << a << "\t" << b << "\t" << S_ab << "\t" << S_cd << "\t" << (*this)(l,0,S_ab,a,b,S_cd,a,b) << "\t";

                  double ward = 0.0;

                  for(int S_al = 0;S_al < 2;++S_al)
                     for(int S_cl = 0;S_cl < 2;++S_cl)
                        ward += std::sqrt( ( 2.0*S_ab + 1.0) * (2*S_cd + 1.0) * (2*S_al + 1.0) * (2*S_cl + 1.0) ) * (1 - 2*S_ab) * (1 - 2*S_cd) * (1 - 2*S_al) * (1 - 2*S_cl)

                        * _6j[S_ab][S_al] * _6j[S_cd][S_cl] * (*this)(b,0,S_al,a,l,S_cl,a,l);

                  cout << ward << "\t";

                  ward = 0.0;

                  for(int S_lb = 0;S_lb < 2;++S_lb)
                     for(int S_ld = 0;S_ld < 2;++S_ld)
                        ward += std::sqrt( ( 2.0*S_ab + 1.0) * (2*S_cd + 1.0) * (2*S_lb + 1.0) * (2*S_ld + 1.0) ) 

                           * _6j[S_ab][S_lb] * _6j[S_cd][S_ld] * (*this)(a,0,S_lb,l,b,S_ld,l,b);

                  cout << ward << endl;

               }

         }

   }

}
