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
 * standard constructor: constructs M rxTPM object with parameter l = 0 -> M-1
 */
dDPM::dDPM() {

   ddpm = new rxTPM * [M];

   for(int l = 0;l < M;++l)
      ddpm[l] = new rxTPM(l);

}

/**
 * copy constructor: constructs an array and copies the content of W in it.
 * @param W object that will be copied into this.
 */
dDPM::dDPM(const dDPM &W) { 

   ddpm = new rxTPM * [M];

   for(int l = 0;l < M;++l)
      ddpm[l] = new rxTPM(W[l]);

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
 * acces to the individual rxTPM objects
 * @param l the specific rxTPM object you want
 * @return the rxTPM object with parameter l
 */
rxTPM &dDPM::operator[](int l){

   return *ddpm[l];

}

/**
 * acces to the individual rxTPM objects: the const version
 * @param l the specific rxTPM object you want
 * @return the rxTPM object with parameter l
 */
const rxTPM &dDPM::operator[](int l) const{

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

   int phase_i = get_inco(l,S,S_ab,a,b);

   if(phase_i == 0)
      return 0.0;

   int phase_j = get_inco(l,S,S_cd,c,d);

   if(phase_j == 0)
      return 0.0;

   return phase_i * phase_j * (*this)[l](S,rxTPM::gs2t(l,S,S_ab,a,b),rxTPM::gs2t(l,S,S_cd,c,d));

}

/**
 * @param l third (block) index
 * @param S dp spin
 * @param S_ab intermediate spin
 * @param a first sp index
 * @param b second sp index
 * @return the right phase for this order of sp indices as a function of my basis.
 */
int dDPM::get_inco(int l,int S,int S_ab,int a,int b){

   //3 indices can never be equal at the same time!
   if(l == a && l == b)
      return 0;

   if(S == 0){//S == 1/2

      if(a == b && S_ab == 1)
         return 0;

      if(a > b)
         return 1 - 2*S_ab;
      else
         return 1;

   }
   else{//S == 3/2

      //here totally antisymmetrical
      if(a == b || a == l || b == l)
         return 0;

      //intermediate spin can never be 0
      if(S_ab == 0)
         return 0;

      if(a > b)
         return -1;
      else
         return 1;

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
 * Pseudo - Invert positive semidefinite symmetric ddpm is stored in (*this), original ddpm (*this) is destroyed
 */
void dDPM::pseudo_invert(){

   for(int l = 0;l < M;++l)
      ddpm[l]->pseudo_invert();

}

/**
 * Pseudo - sqrt of the ddpm object, calls the rxTPM::pseudo_sqrt() function
 * @param option == 1 positive sqrt , if == -1 negative sqrt
 */
void dDPM::pseudo_sqrt(int option){

   for(int l = 0;l < M;++l)
      ddpm[l]->pseudo_sqrt(option);

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

   int phase_i,phase_j;

   //storage for "in block" symmetries
   double mat[2][2];
   double vec[2];

   for(int l = 0;l < M;++l){

      //first the S = 1/2 , the difficult part

      //all equal: W^l|aa;aa: three equalities!
      for(int a = 0;a < M;++a){

         if(a == l)
            a++;

         if(a == M)
            break;

         for(int S_ab = 0;S_ab < 2;++S_ab)
            for(int S_cd = 0;S_cd < 2;++S_cd)
               mat[S_ab][S_cd] = (*this)(l,0,S_ab,l,a,S_cd,l,a);

         for(int S_ab = 0;S_ab < 2;++S_ab)
            for(int S_cd = S_ab;S_cd < 2;++S_cd){

               //1
               double ward = mat[S_ab][S_cd];

               //2
               for(int S_lb = 0;S_lb < 2;++S_lb)
                  ward += std::sqrt( (2.0*S_ab + 1.0) * (2.0*S_lb + 1.0) ) * _6j[S_ab][S_lb] * mat[S_lb][S_cd];

               //3
               for(int S_ld = 0;S_ld < 2;++S_ld)
                  ward += std::sqrt( (2.0*S_cd + 1.0) * (2.0*S_ld + 1.0) ) * _6j[S_cd][S_ld] * mat[S_ab][S_ld];

               //4
               for(int S_lb = 0;S_lb < 2;++S_lb)
                  for(int S_ld = 0;S_ld < 2;++S_ld){

                     ward += std::sqrt( (2.0*S_cd + 1.0) * (2.0*S_ld + 1.0) * (2.0*S_ab + 1.0) * (2.0*S_lb + 1.0) ) 

                        * _6j[S_ab][S_lb] * _6j[S_cd][S_ld] * mat[S_lb][S_ld];

                  }

               //5
               ward += 2.0*std::sqrt( (2.0*S_ab + 1.0) * (2.0*S_cd + 1.0) ) * _6j[S_ab][0] * _6j[S_cd][0] * (1 - 2*S_ab) * (1 - 2*S_cd) * (*this)(a,0,0,l,l,0,l,l);

               i = rxTPM::gs2t(l,0,S_ab,l,a);
               j = rxTPM::gs2t(l,0,S_cd,l,a);

               if(l > a){

                  phase_i = 1 - 2*S_ab;
                  phase_j = 1 - 2*S_cd;

               }
               else{

                  phase_i = 1;
                  phase_j = 1;

               }

               (*this)[l](0,i,j) = 0.2 * phase_i * phase_j * ward;
               (*this)[l](0,j,i) = (*this)[l](0,i,j);

            }

         i = rxTPM::gs2t(a,0,0,l,l);

         (*this)[a](0,i,i) = 0.0;

         for(int S_ab = 0;S_ab < 2;++S_ab)
            for(int S_cd = 0;S_cd < 2;++S_cd){

               (*this)[a](0,i,i) += 0.5*std::sqrt( (2.0*S_ab + 1.0) * (2.0*S_cd + 1.0) ) * _6j[S_ab][0] * _6j[S_cd][0] 

                  * (1 - 2*S_ab) * (1 - 2*S_cd) * (*this)(l,0,S_ab,l,a,S_cd,l,a);

            }

      }

      //other with 3 equalities: a=b, b=c and l = d
      for(int a = 0;a < M;++a){

         if(a == l)
            a++;

         if(a == M)
            break;

         for(int S_cd = 0;S_cd < 2;++S_cd)
            vec[S_cd] = (*this)(l,0,0,a,a,S_cd,a,l);

         //first the "in block symmetries"
         for(int S_cd = 0;S_cd < 2;++S_cd){

            //1)
            double ward = vec[S_cd];

            //2)
            for(int S_cl = 0;S_cl < 2;++S_cl)
               ward += std::sqrt( (2.0*S_cl + 1.0) * (2.0*S_cd + 1.0) ) * (1 - 2*S_cl) * (1 - 2*S_cd) * _6j[S_cd][S_cl] * vec[S_cl];

            //3)
            for(int S_lb = 0;S_lb < 2;++S_lb)
               ward += 2.0*std::sqrt( (2.0*S_cd + 1.0) * (2.0*S_lb + 1.0) ) * _6j[0][S_lb] * _6j[S_cd][0] * (*this)(a,0,S_lb,l,a,0,l,l);

            i = rxTPM::gs2t(l,0,0,a,a);
            j = rxTPM::gs2t(l,0,S_cd,a,l);

            if(a > l)
               phase_j = 1 - 2*S_cd;
            else
               phase_j = 1;

            (*this)[l](0,i,j) = 0.25 * phase_j * ward;
            (*this)[l](0,j,i) = (*this)[l](0,i,j);

         }

         //then the "a" block
         for(int S_lb = 0;S_lb < 2;++S_lb){

            i = rxTPM::gs2t(a,0,S_lb,l,a);
            j = rxTPM::gs2t(a,0,0,l,l);

            if(l > a)
               phase_i = 1 - 2*S_lb;
            else
               phase_i = 1;

            (*this)[a](0,i,j) = 0.0;

            for(int S_cd = 0;S_cd < 2;++S_cd)
               (*this)[a](0,i,j) += phase_i * std::sqrt( (2.0*S_cd + 1.0) * (2.0*S_lb + 1.0) ) * _6j[0][S_lb] * _6j[S_cd][0] * (*this)(l,0,0,a,a,S_cd,a,l);

            (*this)[a](0,j,i) =  (*this)[a](0,i,j);

         }

      }

      //then 2 equalities, first ( a = c ; b = d ) W^l_{ab;ab} <--> W^a_{lb;lb} <--> W^b_{al;al}
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

            //first calculate the avarage:
            for(int S_ab = 0;S_ab < 2;++S_ab)
               for(int S_cd = S_ab;S_cd < 2;++S_cd){

                  double ward = (*this)(l,0,S_ab,a,b,S_cd,a,b);

                  for(int S_lb = 0;S_lb < 2;++S_lb)
                     for(int S_ld = 0;S_ld < 2;++S_ld){

                        ward += std::sqrt( (2.0*S_lb + 1.0) * (2.0*S_ld + 1.0) * (2.0*S_ab + 1.0) * (2.0*S_cd + 1.0) ) * _6j[S_ab][S_lb] * _6j[S_cd][S_ld] 

                           * (*this)(a,0,S_lb,l,b,S_ld,l,b);

                     }

                  for(int S_al = 0;S_al < 2;++S_al)
                     for(int S_cl = 0;S_cl < 2;++S_cl){

                        ward += std::sqrt( (2.0*S_al + 1.0) * (2.0*S_cl + 1.0) * (2.0*S_ab + 1.0) * (2.0*S_cd + 1.0) ) * _6j[S_ab][S_al] * _6j[S_cd][S_cl] 

                           * (1 - 2*S_al) * (1 - 2*S_cl) * (1 - 2*S_ab) * (1 - 2*S_cd) * (*this)(b,0,S_al,a,l,S_cl,a,l);

                     }

                  i = rxTPM::gs2t(l,0,S_ab,a,b);
                  j = rxTPM::gs2t(l,0,S_cd,a,b);

                  (*this)[l](0,i,j) = ward/3.0;
                  (*this)[l](0,j,i) = (*this)[l](0,i,j);

               }

            //then make the rest symmetric
            for(int S_lb = 0;S_lb < 2;++S_lb)
               for(int S_ld = S_lb;S_ld < 2;++S_ld){

                  i = rxTPM::gs2t(a,0,S_lb,l,b);
                  j = rxTPM::gs2t(a,0,S_ld,l,b);

                  if(l > b){

                     phase_i = 1 - 2*S_lb;
                     phase_j = 1 - 2*S_ld;

                  }
                  else{

                     phase_i = 1;
                     phase_j = 1;

                  }

                  (*this)[a](0,i,j) = 0.0;

                  for(int S_ab = 0;S_ab < 2;++S_ab)
                     for(int S_cd = 0;S_cd < 2;++S_cd){

                        (*this)[a](0,i,j) += phase_i * phase_j * std::sqrt( (2.0*S_ab + 1.0) * (2.0*S_cd + 1.0) * (2.0*S_lb + 1.0) * (2.0*S_ld + 1.0) ) 

                           * _6j[S_ab][S_lb] * _6j[S_cd][S_ld] * (*this)(l,0,S_ab,a,b,S_cd,a,b);

                     }

                  (*this)[a](0,j,i) = (*this)[a](0,i,j);

               }

            for(int S_al = 0;S_al < 2;++S_al)
               for(int S_cl = S_al;S_cl < 2;++S_cl){

                  i = rxTPM::gs2t(b,0,S_al,a,l);
                  j = rxTPM::gs2t(b,0,S_cl,a,l);

                  if(a > l){

                     phase_i = 1 - 2*S_al;
                     phase_j = 1 - 2*S_cl;

                  }
                  else{

                     phase_i = 1;
                     phase_j = 1;

                  }

                  (*this)[b](0,i,j) = 0.0;

                  for(int S_ab = 0;S_ab < 2;++S_ab)
                     for(int S_cd = 0;S_cd < 2;++S_cd){

                        (*this)[b](0,i,j) += phase_i * phase_j * std::sqrt( (2.0*S_ab + 1.0) * (2.0*S_cd + 1.0) * (2.0*S_al + 1.0) * (2.0*S_cl + 1.0) ) * 

                           (1 - 2*S_ab) * (1 - 2*S_cd) * (1 - 2*S_al) * (1 - 2*S_cl) * _6j[S_ab][S_al] * _6j[S_cd][S_cl] * (*this)(l,0,S_ab,a,b,S_cd,a,b);

                     }

                  (*this)[b](0,j,i) = (*this)[b](0,i,j);

               }

         }
      }

      //next with two equalities:
      for(int a = 0;a < M;++a){

         if(a == l)
            ++a;

         if(a == M)
            break;

         //a = d and  b = l: W^l_{al;ca}
         for(int c = 0;c < a;++c){

            if(c == l)
               c++;

            if(c == a)
               break;

            for(int S_ab = 0;S_ab < 2;++S_ab)
               for(int S_cd = 0;S_cd < 2;++S_cd)
                  mat[S_ab][S_cd] = (*this)(l,0,S_ab,a,l,S_cd,c,a);

            //first set the "l" block right
            for(int S_ab = 0;S_ab < 2;++S_ab)
               for(int S_cd = 0;S_cd < 2;++S_cd){

                  double ward = mat[S_ab][S_cd];

                  for(int S_al = 0;S_al < 2;++S_al)
                     ward += std::sqrt( (2.0*S_ab + 1.0) * (2.0*S_al + 1.0) ) * (1 - 2*S_al) * (1 - 2*S_ab) * _6j[S_al][S_ab] * mat[S_al][S_cd];

                  for(int S_cl = 0;S_cl < 2;++S_cl){

                     ward += std::sqrt( 2.0 * (2.0*S_ab + 1.0) * (2.0*S_cd + 1.0) * (2.0*S_cl + 1.0) ) * (1 - 2*S_cd) * (1 - 2*S_cl) 

                        * _6j[S_ab][0] * _6j[S_cd][S_cl] * (*this)(a,0,0,l,l,S_cl,c,l);

                  }

                  i = rxTPM::gs2t(l,0,S_ab,a,l);
                  j = rxTPM::gs2t(l,0,S_cd,c,a);

                  if(a > l)
                     phase_i = 1 - 2*S_ab;
                  else
                     phase_i = 1;

                  (*this)[l](0,i,j) = phase_i * ward / 3.0;
                  (*this)[l](0,j,i) = (*this)[l](0,i,j);

               }

            //then the "a" block
            for(int S_cl = 0;S_cl < 2;++S_cl){

               i = rxTPM::gs2t(a,0,0,l,l);
               j = rxTPM::gs2t(a,0,S_cl,c,l);

               if(c > l)
                  phase_j = 1 - 2*S_cl;
               else
                  phase_j = 1;

               (*this)[a](0,i,j) = 0.0;

               for(int S_ab = 0;S_ab < 2;++S_ab)
                  for(int S_cd = 0;S_cd < 2;++S_cd){

                     (*this)[a](0,i,j) += phase_j * std::sqrt( 0.5 * (2.0*S_ab + 1.0) * (2.0*S_cd + 1.0) * (2.0*S_cl + 1.0) ) * (1 - 2*S_cl) * (1 - 2*S_cd)

                        * _6j[S_ab][0] * _6j[S_cl][S_cd] * (*this)(l,0,S_ab,a,l,S_cd,c,a);

                  }

               (*this)[a](0,j,i) = (*this)[a](0,i,j);

            }

         }

         //a = c and  b = l: W^l_{al;ad}
         for(int d = a + 1;d < M;++d){

            if(d == l)
               d++;

            if(d == M)
               break;

            for(int S_ab = 0;S_ab < 2;++S_ab)
               for(int S_cd = 0;S_cd < 2;++S_cd)
                  mat[S_ab][S_cd] = (*this)(l,0,S_ab,a,l,S_cd,a,d);

            //first set the "l" block right
            for(int S_ab = 0;S_ab < 2;++S_ab)
               for(int S_cd = 0;S_cd < 2;++S_cd){

                  double ward = mat[S_ab][S_cd];

                  for(int S_al = 0;S_al < 2;++S_al)
                     ward += std::sqrt( (2.0*S_ab + 1.0) * (2.0*S_al + 1.0) ) * (1 - 2*S_al) * (1 - 2*S_ab) * _6j[S_al][S_ab] * mat[S_al][S_cd];

                  for(int S_ld = 0;S_ld < 2;++S_ld){

                     ward += std::sqrt( 2.0 * (2.0*S_ab + 1.0) * (2.0*S_cd + 1.0) * (2.0*S_ld + 1.0) )

                        * _6j[S_ab][0] * _6j[S_cd][S_ld] * (*this)(a,0,0,l,l,S_ld,l,d);

                  }

                  i = rxTPM::gs2t(l,0,S_ab,a,l);
                  j = rxTPM::gs2t(l,0,S_cd,a,d);

                  if(a > l)
                     phase_i = 1 - 2*S_ab;
                  else
                     phase_i = 1;

                  (*this)[l](0,i,j) = phase_i * ward / 3.0;
                  (*this)[l](0,j,i) = (*this)[l](0,i,j);

               }

            //then the "a" block
            for(int S_ld = 0;S_ld < 2;++S_ld){

               i = rxTPM::gs2t(a,0,0,l,l);
               j = rxTPM::gs2t(a,0,S_ld,l,d);

               if(l > d)
                  phase_j = 1 - 2*S_ld;
               else
                  phase_j = 1;

               (*this)[a](0,i,j) = 0.0;

               for(int S_ab = 0;S_ab < 2;++S_ab)
                  for(int S_cd = 0;S_cd < 2;++S_cd){

                     (*this)[a](0,i,j) += phase_j * std::sqrt( 0.5 * (2.0*S_ab + 1.0) * (2.0*S_cd + 1.0) * (2.0*S_ld + 1.0) ) 

                        * _6j[S_ab][0] * _6j[S_ld][S_cd] * (*this)(l,0,S_ab,a,l,S_cd,a,d);

                  }

               (*this)[a](0,j,i) = (*this)[a](0,i,j);

            }

         }

      }

      //another one with 2 equalities: a = c = l
      for(int b = l + 1;b < M;++b)
         for(int d = b + 1;d < M;++d){

            //save the numbers first in the matrix mat
            for(int S_ab = 0;S_ab < 2;++S_ab)
               for(int S_cd = 0;S_cd < 2;++S_cd)
                  mat[S_ab][S_cd] = (*this)(l,0,S_ab,l,b,S_cd,l,d);

            for(int S_ab = 0;S_ab < 2;++S_ab)
               for(int S_cd = 0;S_cd < 2;++S_cd){

                  i = rxTPM::gs2t(l,0,S_ab,l,b);
                  j = rxTPM::gs2t(l,0,S_cd,l,d);

                  double ward = mat[S_ab][S_cd];

                  for(int S_lb = 0;S_lb < 2;++S_lb)
                     ward += std::sqrt( (2.0*S_ab + 1.0) * (2.0*S_lb + 1.0) ) * _6j[S_ab][S_lb] * mat[S_lb][S_cd];

                  for(int S_ld = 0;S_ld < 2;++S_ld)
                     ward += std::sqrt( (2.0*S_cd + 1.0) * (2.0*S_ld + 1.0) ) * _6j[S_cd][S_ld] * mat[S_ab][S_ld];

                  for(int S_lb = 0;S_lb < 2;++S_lb)
                     for(int S_ld = 0;S_ld < 2;++S_ld){

                        ward += std::sqrt( (2.0*S_ab + 1.0) * (2.0*S_cd + 1.0) * (2.0*S_lb + 1.0) * (2.0*S_ld + 1.0) )

                           * _6j[S_ab][S_lb] * _6j[S_cd][S_ld] * mat[S_lb][S_ld];

                     }

                  (*this)[l](0,i,j) = 0.25 * ward;
                  (*this)[l](0,j,i) = (*this)[l](0,i,j);

               }

         }

      //next one with 2 equalities: b = c = l
      for(int a = 0;a < l;++a)
         for(int d = l + 1;d < M;++d){

            //save the numbers first in the matrix mat
            for(int S_ab = 0;S_ab < 2;++S_ab)
               for(int S_cd = 0;S_cd < 2;++S_cd)
                  mat[S_ab][S_cd] = (*this)(l,0,S_ab,a,l,S_cd,l,d);

            for(int S_ab = 0;S_ab < 2;++S_ab)
               for(int S_cd = 0;S_cd < 2;++S_cd){

                  i = rxTPM::gs2t(l,0,S_ab,a,l);
                  j = rxTPM::gs2t(l,0,S_cd,l,d);

                  double ward = mat[S_ab][S_cd];

                  for(int S_al = 0;S_al < 2;++S_al)
                     ward += std::sqrt( (2.0*S_ab + 1.0) * (2.0*S_al + 1.0) ) * (1 - 2*S_al) * (1 - 2*S_ab) * _6j[S_ab][S_al] * mat[S_al][S_cd];

                  for(int S_ld = 0;S_ld < 2;++S_ld)
                     ward += std::sqrt( (2.0*S_cd + 1.0) * (2.0*S_ld + 1.0) ) * _6j[S_cd][S_ld] * mat[S_ab][S_ld];

                  for(int S_al = 0;S_al < 2;++S_al)
                     for(int S_ld = 0;S_ld < 2;++S_ld){

                        ward += std::sqrt( (2.0*S_ab + 1.0) * (2.0*S_cd + 1.0) * (2.0*S_al + 1.0) * (2.0*S_ld + 1.0) ) * (1 - 2*S_al) * (1 - 2*S_ab)

                           * _6j[S_ab][S_al] * _6j[S_cd][S_ld] * mat[S_al][S_ld];

                     }

                  (*this)[l](0,i,j) = 0.25*ward;
                  (*this)[l](0,j,i) = (*this)[l](0,i,j);

               }

         }

      //next one with 2 equalities: b = d = l
      for(int a = 0;a < l;++a)
         for(int c = a + 1;c < l;++c){

            //save the numbers first in the matrix mat
            for(int S_ab = 0;S_ab < 2;++S_ab)
               for(int S_cd = 0;S_cd < 2;++S_cd)
                  mat[S_ab][S_cd] = (*this)(l,0,S_ab,a,l,S_cd,c,l);

            for(int S_ab = 0;S_ab < 2;++S_ab)
               for(int S_cd = 0;S_cd < 2;++S_cd){

                  i = rxTPM::gs2t(l,0,S_ab,a,l);
                  j = rxTPM::gs2t(l,0,S_cd,c,l);

                  double ward = mat[S_ab][S_cd];

                  for(int S_al = 0;S_al < 2;++S_al)
                     ward += std::sqrt( (2.0*S_ab + 1.0) * (2.0*S_al + 1.0) ) * (1 - 2*S_al) * (1 - 2*S_ab) * _6j[S_ab][S_al] * mat[S_al][S_cd];

                  for(int S_cl = 0;S_cl < 2;++S_cl)
                     ward += std::sqrt( (2.0*S_cd + 1.0) * (2.0*S_cl + 1.0) ) * (1 - 2*S_cd) * (1 - 2*S_cl) * _6j[S_cd][S_cl] * mat[S_ab][S_cl];

                  for(int S_al = 0;S_al < 2;++S_al)
                     for(int S_cl = 0;S_cl < 2;++S_cl){

                        ward += std::sqrt( (2.0*S_ab + 1.0) * (2.0*S_cd + 1.0) * (2.0*S_al + 1.0) * (2.0*S_cl + 1.0) ) * (1 - 2*S_al) * (1 - 2*S_ab)

                           * (1 - 2*S_cd) * (1 - 2*S_cl) * _6j[S_ab][S_al] * _6j[S_cd][S_cl] * mat[S_al][S_cl];

                     }

                  (*this)[l](0,i,j) = 0.25*ward;
                  (*this)[l](0,j,i) = (*this)[l](0,i,j);

               }

         }

      //one equality: start with a == c
      for(int a = 0;a < M;++a){

         if(a == l)
            a++;

         if(a == M)
            break;

         for(int b = a + 1;b < M;++b){

            if(b == l)
               b++;

            if(b == M)
               break;

            for(int d = b + 1;d < M;++d){

               if(d == l)
                  d++;

               if(d == M)
                  break;

               for(int S_ab = 0;S_ab < 2;++S_ab)
                  for(int S_cd = 0;S_cd < 2;++S_cd){

                     //first take the average
                     double ward = (*this)(l,0,S_ab,a,b,S_cd,a,d);

                     for(int S_lb = 0;S_lb < 2;++S_lb)
                        for(int S_ld = 0;S_ld < 2;++S_ld){

                           ward += std::sqrt( (2.0*S_ab + 1.0) * (2.0*S_cd + 1.0) * (2.0*S_lb + 1.0) * (2.0*S_ld + 1.0) ) * _6j[S_ab][S_lb] * _6j[S_cd][S_ld]

                              * (*this)(a,0,S_lb,l,b,S_ld,l,d);

                        }

                     i = rxTPM::gs2t(l,0,S_ab,a,b);
                     j = rxTPM::gs2t(l,0,S_cd,a,d);

                     (*this)[l](0,i,j) = 0.5*ward;
                     (*this)[l](0,j,i) = (*this)[l](0,i,j);

                  }

               //then symmetrize the other term
               for(int S_lb = 0;S_lb < 2;++S_lb)
                  for(int S_ld = 0;S_ld < 2;++S_ld){

                     i = rxTPM::gs2t(a,0,S_lb,l,b);
                     j = rxTPM::gs2t(a,0,S_ld,l,d);

                     if(l > b)
                        phase_i = 1 - 2*S_lb;
                     else
                        phase_i = 1;

                     if(l > d)
                        phase_j = 1 - 2*S_ld;
                     else
                        phase_j = 1;

                     (*this)[a](0,i,j) = 0.0;

                     for(int S_ab = 0;S_ab < 2;++S_ab)
                        for(int S_cd = 0;S_cd < 2;++S_cd){

                           (*this)[a](0,i,j) += phase_i * phase_j * std::sqrt( (2.0*S_ab + 1.0) * (2.0*S_cd + 1.0) * (2.0*S_lb + 1.0) * (2.0*S_ld + 1.0) )

                              * _6j[S_ab][S_lb] * _6j[S_cd][S_ld] * (*this)(l,0,S_ab,a,b,S_cd,a,d);

                        }

                     (*this)[a](0,j,i) = (*this)[a](0,i,j);

                  }

            }
         }
      }


      //then b == c: ab;bd
      for(int a = 0;a < M;++a){

         if(a == l)
            a++;

         if(a == M)
            break;

         for(int b = a + 1;b < M;++b){

            if(b == l)
               b++;

            if(b == M)
               break;

            for(int d = b + 1;d < M;++d){

               if(d == l)
                  d++;

               if(d == M)
                  break;

               //first take average
               for(int S_ab = 0;S_ab < 2;++S_ab)
                  for(int S_cd = 0;S_cd < 2;++S_cd){

                     double ward = (*this)(l,0,S_ab,a,b,S_cd,b,d);

                     for(int S_al = 0;S_al < 2;++S_al)
                        for(int S_ld = 0;S_ld < 2;++S_ld){

                           ward += std::sqrt( (2.0*S_ab + 1.0) * (2.0*S_cd + 1.0) * (2.0*S_al + 1.0) * (2.0*S_ld + 1.0) ) * (1 - 2*S_ab) * (1 - 2*S_al)

                              * _6j[S_al][S_ab] * _6j[S_cd][S_ld] * (*this)(b,0,S_al,a,l,S_ld,l,d);

                        }

                     i = rxTPM::gs2t(l,0,S_ab,a,b);
                     j = rxTPM::gs2t(l,0,S_cd,b,d);

                     (*this)[l](0,i,j) = 0.5 * ward;
                     (*this)[l](0,j,i) = (*this)[l](0,i,j);

                  }

               //then symmetrize rest
               for(int S_al = 0;S_al < 2;++S_al)
                  for(int S_ld = 0;S_ld < 2;++S_ld){

                     i = rxTPM::gs2t(b,0,S_al,a,l);
                     j = rxTPM::gs2t(b,0,S_ld,l,d);

                     if(a > l)
                        phase_i = 1 - 2*S_al;
                     else
                        phase_i = 1;

                     if(l > d)
                        phase_j = 1 - 2*S_ld;
                     else
                        phase_j = 1;

                     (*this)[b](0,i,j) = 0.0;

                     for(int S_ab = 0;S_ab < 2;++S_ab)
                        for(int S_cd = 0;S_cd < 2;++S_cd){

                           (*this)[b](0,i,j) += phase_i * phase_j * std::sqrt( (2.0*S_ab + 1.0) * (2.0*S_cd + 1.0) * (2.0*S_al + 1.0) * (2.0*S_ld + 1.0) ) 

                              * (1 - 2*S_ab) * (1 - 2*S_al) * _6j[S_al][S_ab] * _6j[S_cd][S_ld] * (*this)(l,0,S_ab,a,b,S_cd,b,d);

                        }

                     (*this)[b](0,j,i) = (*this)[b](0,i,j);

                  }

            }
         }
      }

      //last regular: b == d --> W^l_{ab;cb}
      for(int a = 0;a < M;++a){

         if(a == l)
            a++;

         if(a == M)
            break;

         for(int c = a + 1;c < M;++c){

            if(c == l)
               c++;

            if(c == M)
               break;

            for(int b = c + 1;b < M;++b){

               if(b == l)
                  b++;

               if(b == M)
                  break;

               //first average out
               for(int S_ab = 0;S_ab < 2;++S_ab)
                  for(int S_cd = 0;S_cd < 2;++S_cd){

                     double ward = (*this)(l,0,S_ab,a,b,S_cd,c,b);

                     for(int S_al = 0;S_al < 2;++S_al)
                        for(int S_cl = 0;S_cl < 2;++S_cl){

                           ward += std::sqrt( (2.0*S_al + 1.0) * (2.0*S_cl + 1.0) * (2.0*S_ab + 1.0) * (2.0*S_cd + 1.0) ) * (1 - 2*S_ab) * (1 - 2*S_cd) 

                              * (1 - 2*S_al) * (1 - 2*S_cl) * _6j[S_al][S_ab] * _6j[S_cl][S_cd] * (*this)(b,0,S_al,a,l,S_cl,c,l);

                        }

                     i = rxTPM::gs2t(l,0,S_ab,a,b);
                     j = rxTPM::gs2t(l,0,S_cd,c,b);

                     (*this)[l](0,i,j) = 0.5*ward;
                     (*this)[l](0,j,i) = (*this)[l](0,i,j);

                  }

               //then make the rest symmetric
               for(int S_al = 0;S_al < 2;++S_al)
                  for(int S_cl = 0;S_cl < 2;++S_cl){

                     i = rxTPM::gs2t(b,0,S_al,a,l);
                     j = rxTPM::gs2t(b,0,S_cl,c,l);

                     if(a > l)
                        phase_i = 1 - 2*S_al;
                     else
                        phase_i = 1;

                     if(c > l)
                        phase_j = 1 - 2*S_cl;
                     else
                        phase_j = 1;

                     (*this)[b](0,i,j) = 0.0;

                     for(int S_ab = 0;S_ab < 2;++S_ab)
                        for(int S_cd = 0;S_cd < 2;++S_cd){

                           (*this)[b](0,i,j) += phase_i*phase_j * std::sqrt( (2.0*S_al + 1.0) * (2.0*S_cl + 1.0) * (2.0*S_ab + 1.0) * (2.0*S_cd + 1.0) ) 

                              * (1 - 2*S_ab) * (1 - 2*S_cd) * (1 - 2*S_al ) * (1 - 2*S_cl) * _6j[S_al][S_ab] * _6j[S_cl][S_cd] * (*this)(l,0,S_ab,a,b,S_cd,c,b);

                        }

                     (*this)[b](0,j,i) = (*this)[b](0,i,j);

                  }

            }
         }
      }

      //two more terms: b == l
      for(int a = 0;a < l;++a){

         for(int S_cd = 0;S_cd < 2;++S_cd){

            for(int c = 0;c < M;++c){

               if(c != l && c != a){

                  for(int d = c + S_cd;d < M;++d){

                     if(d != l && d != a){

                        for(int S_ab = 0;S_ab < 2;++S_ab)
                           vec[S_ab] = (*this)(l,0,S_ab,a,l,S_cd,c,d);

                        for(int S_ab = 0;S_ab < 2;++S_ab){

                           double ward = vec[S_ab];

                           for(int S_al = 0;S_al < 2;++S_al)
                              ward += std::sqrt( (2.0*S_al + 1.0) * (2.0*S_ab + 1.0) ) * (1 - 2*S_al) * (1 - 2*S_ab) * _6j[S_al][S_ab] * vec[S_al];

                           i = rxTPM::gs2t(l,0,S_ab,a,l);
                           j = rxTPM::gs2t(l,0,S_cd,c,d);

                           (*this)[l](0,i,j) = 0.5 * ward;
                           (*this)[l](0,j,i) = (*this)[l](0,i,j);

                        }

                     }

                  }

               }

            }

         }

      }

      //and at last: a == l
      for(int b = l + 1;b < M;++b){

         for(int S_cd = 0;S_cd < 2;++S_cd){

            for(int c = 0;c < M;++c){

               if(c != l && c != b){

                  for(int d = c + S_cd;d < M;++d){

                     if(d != l && d != b){

                        for(int S_ab = 0;S_ab < 2;++S_ab)
                           vec[S_ab] = (*this)(l,0,S_ab,l,b,S_cd,c,d);

                        for(int S_ab = 0;S_ab < 2;++S_ab){

                           double ward = vec[S_ab];

                           for(int S_lb = 0;S_lb < 2;++S_lb)
                              ward += std::sqrt( (2.0*S_lb + 1.0) * (2.0*S_ab + 1.0) ) * _6j[S_lb][S_ab] * vec[S_lb];

                           i = rxTPM::gs2t(l,0,S_ab,l,b);
                           j = rxTPM::gs2t(l,0,S_cd,c,d);

                           (*this)[l](0,i,j) = 0.5 * ward;
                           (*this)[l](0,j,i) = (*this)[l](0,i,j);

                        }

                     }

                  }

               }

            }

         }

      }

      //easy part: S = 3/2

      //first diagonal part: ab;ab
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

            i = rxTPM::gs2t(l,1,1,a,b);

            (*this)[l](1,i,i) = 1.0/3.0 * ( (*this)[l](1,i,i) + (*this)(a,1,1,l,b,1,l,b) + (*this)(b,1,1,a,l,1,a,l) );

            //rest is symmetric
            i = rxTPM::gs2t(a,1,1,l,b);

            (*this)[a](1,i,i) = (*this)(l,1,1,a,b,1,a,b);

            i = rxTPM::gs2t(b,1,1,a,l);

            (*this)[b](1,i,i) = (*this)(l,1,1,a,b,1,a,b);

         }
      }

      //three terms with one diagonal index:

      //1) b == d: ab;cb
      for(int a = 0;a < M;++a){

         if(a == l)
            a++;

         if(a == M)
            break;

         for(int c = a + 1;c < M;++c){

            if(c == l)
               c++;

            if(c == M)
               break;

            for(int b = c + 1;b < M;++b){

               if(b == l)
                  b++;

               if(b == M)
                  break;

               i = rxTPM::gs2t(l,1,1,a,b);
               j = rxTPM::gs2t(l,1,1,c,b);

               (*this)[l](1,i,j) = 0.5 * ( (*this)[l](1,i,j) + (*this)(b,1,1,a,l,1,c,l) );
               (*this)[l](1,j,i) = (*this)[l](1,i,j);

               i = rxTPM::gs2t(b,1,1,a,l);
               j = rxTPM::gs2t(b,1,1,c,l);

               if(a > l)
                  phase_i = -1;
               else
                  phase_i = 1;

               if(c > l)
                  phase_j = -1;
               else
                  phase_j = 1;

               (*this)[b](1,i,j) = phase_i*phase_j*(*this)(l,1,1,a,b,1,c,b);
               (*this)[b](1,j,i) = (*this)[b](1,i,j);

            }
         }
      }

      //2) b == c: ab;bc
      for(int a = 0;a < M;++a){

         if(a == l)
            a++;

         if(a == M)
            break;

         for(int b = a + 1;b < M;++b){

            if(b == l)
               b++;

            if(b == M)
               break;

            for(int c = b + 1;c < M;++c){

               if(c == l)
                  c++;

               if(c == M)
                  break;

               i = rxTPM::gs2t(l,1,1,a,b);
               j = rxTPM::gs2t(l,1,1,b,c);

               (*this)[l](1,i,j) = 0.5 * ( (*this)[l](1,i,j) + (*this)(b,1,1,a,l,1,l,c) );
               (*this)[l](1,j,i) = (*this)[l](1,i,j);

               i = rxTPM::gs2t(b,1,1,a,l);
               j = rxTPM::gs2t(b,1,1,l,c);

               if(a > l)
                  phase_i = -1;
               else
                  phase_i = 1;

               if(l > c)
                  phase_j = -1;
               else
                  phase_j = 1;

               (*this)[b](1,i,j) = phase_i*phase_j*(*this)(l,1,1,a,b,1,b,c);
               (*this)[b](1,j,i) = (*this)[b](1,i,j);

            }
         }
      }

      //3) a == c: ab;ac
      for(int a = 0;a < M;++a){

         if(a == l)
            a++;

         if(a == M)
            break;

         for(int b = a + 1;b < M;++b){

            if(b == l)
               b++;

            if(b == M)
               break;

            for(int c = b + 1;c < M;++c){

               if(c == l)
                  c++;

               if(c == M)
                  break;

               i = rxTPM::gs2t(l,1,1,a,b);
               j = rxTPM::gs2t(l,1,1,a,c);

               (*this)[l](1,i,j) = 0.5 * ( (*this)[l](1,i,j) + (*this)(a,1,1,l,b,1,l,c) );
               (*this)[l](1,j,i) = (*this)[l](1,i,j);

               i = rxTPM::gs2t(a,1,1,l,b);
               j = rxTPM::gs2t(a,1,1,l,c);

               if(l > b)
                  phase_i = -1;
               else
                  phase_i = 1;

               if(l > c)
                  phase_j = -1;
               else
                  phase_j = 1;

               (*this)[a](1,i,j) = phase_i * phase_j * (*this)(l,1,1,a,b,1,a,c);
               (*this)[a](1,j,i) = (*this)[a](1,i,j);

            }
         }
      }

   }

}

/**
 * project *this on the traceless space but staying in the right symmetry area.
 */
void dDPM::proj_Tr(){

   double ward = this->ddot(Tools::gunit())/Tools::gunit().ddot(Tools::gunit());

   this->daxpy(-ward,Tools::gunit());

}


/**
 * total projection, both on traceless space and good third index symmetry
 */
void dDPM::proj(){

   this->proj_W();
   this->proj_Tr();

}

/**
 * test if the projection is correct, number 1 -> print all
 */
void dDPM::test_proj_1() const {

   for(int l = 0;l < M;++l){

      cout << endl;
      cout << "l = \t" << l << endl;
      cout << endl;

      cout << "S = 1/2" << endl;
      cout << endl;

      for(int a = 0;a < M;++a)
         for(int b = 0;b < M;++b)
            for(int c = 0;c < M;++c){

               cout << endl;

               //1) b = d
               for(int S_ab = 0;S_ab < 2;++S_ab)
                  for(int S_cd = 0;S_cd < 2;++S_cd){

                     cout << "ab;cb\t" << a << "\t" << b << "\t" << c << "\t(" << S_ab << ")\t(" << S_cd << ")\t" << (*this)(l,0,S_ab,a,b,S_cd,c,b) << "\t";

                     double ward = 0.0;

                     for(int S_al = 0;S_al < 2;++S_al)
                        for(int S_cl = 0;S_cl < 2;++S_cl){

                           double hard = std::sqrt( (2.0*S_al + 1.0) * (2.0*S_cl + 1.0) * (2.0*S_ab + 1.0) * (2.0*S_cd + 1.0) ) * (1 - 2*S_al) * (1 - 2*S_cl)

                              * (1 - 2*S_ab) * (1 - 2*S_cd) * _6j[S_ab][S_al] * _6j[S_cd][S_cl] * (*this)(b,0,S_al,a,l,S_cl,c,l);

                           if(a == b)
                              hard /= std::sqrt(2.0);

                           if(c == b)
                              hard /= std::sqrt(2.0);

                           if(a == l)
                              hard *= std::sqrt(2.0);

                           if(c == l)
                              hard *= std::sqrt(2.0);

                           ward += hard;

                        }

                     cout << ward << endl;

                  }

               //2) b = c
               for(int S_ab = 0;S_ab < 2;++S_ab)
                  for(int S_cd = 0;S_cd < 2;++S_cd){

                     cout << "ab;bc\t" << a << "\t" << b << "\t" << c << "\t(" << S_ab << ")\t(" << S_cd << ")\t" << (*this)(l,0,S_ab,a,b,S_cd,b,c) << "\t";

                     double ward = 0.0;

                     for(int S_al = 0;S_al < 2;++S_al)
                        for(int S_ld = 0;S_ld < 2;++S_ld){

                           double hard = std::sqrt( (2.0*S_al + 1.0) * (2.0*S_ld + 1.0) * (2.0*S_ab + 1.0) * (2.0*S_cd + 1.0) ) * (1 - 2*S_al) * 

                              (1 - 2*S_ab) * _6j[S_ab][S_al] * _6j[S_cd][S_ld] * (*this)(b,0,S_al,a,l,S_ld,l,c);

                           if(a == b)
                              hard /= std::sqrt(2.0);

                           if(c == b)
                              hard /= std::sqrt(2.0);

                           if(a == l)
                              hard *= std::sqrt(2.0);

                           if(c == l)
                              hard *= std::sqrt(2.0);

                           ward += hard;

                        }

                     cout << ward << endl;

                  }

               //3) a = d
               for(int S_ab = 0;S_ab < 2;++S_ab)
                  for(int S_cd = 0;S_cd < 2;++S_cd){

                     cout << "ab;ca\t" << a << "\t" << b << "\t" << c << "\t(" << S_ab << ")\t(" << S_cd << ")\t" << (*this)(l,0,S_ab,a,b,S_cd,c,a) << "\t";

                     double ward = 0.0;

                     for(int S_lb = 0;S_lb < 2;++S_lb)
                        for(int S_cl = 0;S_cl < 2;++S_cl){

                           double hard = std::sqrt( (2.0*S_lb + 1.0) * (2.0*S_cl + 1.0) * (2.0*S_ab + 1.0) * (2.0*S_cd + 1.0) ) * (1 - 2*S_cl)

                              * (1 - 2*S_cd) * _6j[S_ab][S_lb] * _6j[S_cd][S_cl] * (*this)(a,0,S_lb,l,b,S_cl,c,l);

                           if(a == b)
                              hard /= std::sqrt(2.0);

                           if(c == a)
                              hard /= std::sqrt(2.0);

                           if(l == b)
                              hard *= std::sqrt(2.0);

                           if(c == l)
                              hard *= std::sqrt(2.0);

                           ward += hard;

                        }

                     cout << ward << endl;

                  }

               //4) a = c
               for(int S_ab = 0;S_ab < 2;++S_ab)
                  for(int S_cd = 0;S_cd < 2;++S_cd){

                     cout << "ab;ac\t" << a << "\t" << b << "\t" << c << "\t(" << S_ab << ")\t(" << S_cd << ")\t" << (*this)(l,0,S_ab,a,b,S_cd,a,c) << "\t";

                     double ward = 0.0;

                     for(int S_lb = 0;S_lb < 2;++S_lb)
                        for(int S_ld = 0;S_ld < 2;++S_ld){

                           double hard = std::sqrt( (2.0*S_lb + 1.0) * (2.0*S_ld + 1.0) * (2.0*S_ab + 1.0) * (2.0*S_cd + 1.0) ) 

                              * _6j[S_ab][S_lb] * _6j[S_cd][S_ld] * (*this)(a,0,S_lb,l,b,S_ld,l,c);

                           if(a == b)
                              hard /= std::sqrt(2.0);

                           if(c == a)
                              hard /= std::sqrt(2.0);

                           if(l == b)
                              hard *= std::sqrt(2.0);

                           if(c == l)
                              hard *= std::sqrt(2.0);

                           ward += hard;

                        }

                     cout << ward << endl;

                  }

               cout << endl;

               //)5 a = l
               for(int S_ab = 0;S_ab < 2;++S_ab)
                  for(int S_cd = 0;S_cd < 2;++S_cd){

                     cout << "lb;cd\t" << a << "\t" << b << "\t" << c << "\t(" << S_ab << ")\t(" << S_cd << ")\t" << (*this)(l,0,S_ab,l,a,S_cd,b,c) << "\t";

                     double ward = 0.0;

                     for(int S_lb = 0;S_lb < 2;++S_lb)
                        ward += std::sqrt( (2.0*S_ab + 1.0) * (2.0*S_lb + 1.0) ) * _6j[S_lb][S_ab] * (*this)(l,0,S_lb,l,a,S_cd,b,c);

                     cout << ward << endl;

                  }

               //6) b = l
               for(int S_ab = 0;S_ab < 2;++S_ab)
                  for(int S_cd = 0;S_cd < 2;++S_cd){

                     cout << "al;cd\t" << a << "\t" << b << "\t" << c << "\t(" << S_ab << ")\t(" << S_cd << ")\t" << (*this)(l,0,S_ab,a,l,S_cd,b,c) << "\t";

                     double ward = 0.0;

                     for(int S_al = 0;S_al < 2;++S_al)
                        ward += std::sqrt( (2.0*S_ab + 1.0) * (2.0*S_al + 1.0) ) * (1 - 2*S_ab) * (1 - 2*S_al) * _6j[S_al][S_ab] * (*this)(l,0,S_al,a,l,S_cd,b,c);

                     cout << ward << endl;

                  }

               //7) c = l
               for(int S_ab = 0;S_ab < 2;++S_ab)
                  for(int S_cd = 0;S_cd < 2;++S_cd){

                     cout << "ab;ld\t" << a << "\t" << b << "\t" << c << "\t(" << S_ab << ")\t(" << S_cd << ")\t" << (*this)(l,0,S_ab,a,b,S_cd,l,c) << "\t";

                     double ward = 0.0;

                     for(int S_ld = 0;S_ld < 2;++S_ld)
                        ward += std::sqrt( (2.0*S_cd + 1.0) * (2.0*S_ld + 1.0) ) * _6j[S_cd][S_ld] * (*this)(l,0,S_ab,a,b,S_ld,l,c);

                     cout << ward << endl;

                  }

               //8) d = l
               for(int S_ab = 0;S_ab < 2;++S_ab)
                  for(int S_cd = 0;S_cd < 2;++S_cd){

                     cout << "ab;cl\t" << a << "\t" << b << "\t" << c << "\t(" << S_ab << ")\t(" << S_cd << ")\t" << (*this)(l,0,S_ab,a,b,S_cd,c,l) << "\t";

                     double ward = 0.0;

                     for(int S_cl = 0;S_cl < 2;++S_cl)
                        ward += std::sqrt( (2.0*S_cd + 1.0) * (2.0*S_cl + 1.0) ) * (1 - 2*S_cd) * (1 - 2*S_cl) * _6j[S_cd][S_cl] * (*this)(l,0,S_ab,a,b,S_cl,c,l);

                     cout << ward << endl;

                  }

            }

      cout << endl;


      cout << endl;
      cout << "S = 3/2" << endl;
      cout << endl;

      for(int a = 0;a < M;++a)
         for(int b = 0;b < M;++b)
            for(int c = 0;c < M;++c){

               cout << endl;
               cout << a << "\t" << b << "\t" << c << "\t|\t" << (*this)(l,1,1,a,b,1,c,b) << "\t" << (*this)(b,1,1,a,l,1,c,l) << endl;
               cout << a << "\t" << b << "\t" << c << "\t|\t" << (*this)(l,1,1,a,b,1,b,c) << "\t" << (*this)(b,1,1,a,l,1,l,c) << endl;
               cout << a << "\t" << b << "\t" << c << "\t|\t" << (*this)(l,1,1,a,b,1,c,a) << "\t" << (*this)(a,1,1,l,b,1,c,l) << endl;
               cout << a << "\t" << b << "\t" << c << "\t|\t" << (*this)(l,1,1,a,b,1,a,c) << "\t" << (*this)(a,1,1,l,b,1,l,c) << endl;
               cout << endl;

            }

   }

}

/**
 * test if the projection is correct, number 2 -> print only the errors
 */
void dDPM::test_proj_2() const {

   for(int l = 0;l < M;++l){

      cout << endl;
      cout << "l = \t" << l << endl;
      cout << endl;

      cout << "S = 1/2" << endl;
      cout << endl;

      for(int a = 0;a < M;++a)
         for(int b = 0;b < M;++b)
            for(int c = 0;c < M;++c){

               //1) b = d
               for(int S_ab = 0;S_ab < 2;++S_ab)
                  for(int S_cd = 0;S_cd < 2;++S_cd){

                     double ward = 0.0;

                     for(int S_al = 0;S_al < 2;++S_al)
                        for(int S_cl = 0;S_cl < 2;++S_cl){

                           double hard = std::sqrt( (2.0*S_al + 1.0) * (2.0*S_cl + 1.0) * (2.0*S_ab + 1.0) * (2.0*S_cd + 1.0) ) * (1 - 2*S_al) * (1 - 2*S_cl)

                              * (1 - 2*S_ab) * (1 - 2*S_cd) * _6j[S_ab][S_al] * _6j[S_cd][S_cl] * (*this)(b,0,S_al,a,l,S_cl,c,l);

                           if(a == b)
                              hard /= std::sqrt(2.0);

                           if(c == b)
                              hard /= std::sqrt(2.0);

                           if(a == l)
                              hard *= std::sqrt(2.0);

                           if(c == l)
                              hard *= std::sqrt(2.0);

                           ward += hard;

                        }

                     if(fabs(ward - (*this)(l,0,S_ab,a,b,S_cd,c,b)) > 1.0e-12)
                        cout << "ab;cb\t" << a << "\t" << b << "\t" << c << "\t(" << S_ab << ")\t(" << S_cd << ")\t" << (*this)(l,0,S_ab,a,b,S_cd,c,b) << "\t" << ward << endl;

                  }

               //2) b = c
               for(int S_ab = 0;S_ab < 2;++S_ab)
                  for(int S_cd = 0;S_cd < 2;++S_cd){

                     double ward = 0.0;

                     for(int S_al = 0;S_al < 2;++S_al)
                        for(int S_ld = 0;S_ld < 2;++S_ld){

                           double hard = std::sqrt( (2.0*S_al + 1.0) * (2.0*S_ld + 1.0) * (2.0*S_ab + 1.0) * (2.0*S_cd + 1.0) ) * (1 - 2*S_al) * 

                              (1 - 2*S_ab) * _6j[S_ab][S_al] * _6j[S_cd][S_ld] * (*this)(b,0,S_al,a,l,S_ld,l,c);

                           if(a == b)
                              hard /= std::sqrt(2.0);

                           if(c == b)
                              hard /= std::sqrt(2.0);

                           if(a == l)
                              hard *= std::sqrt(2.0);

                           if(c == l)
                              hard *= std::sqrt(2.0);

                           ward += hard;

                        }

                     if(fabs(ward - (*this)(l,0,S_ab,a,b,S_cd,b,c)) > 1.0e-12)
                        cout << "ab;bc\t" << a << "\t" << b << "\t" << c << "\t(" << S_ab << ")\t(" << S_cd << ")\t" << (*this)(l,0,S_ab,a,b,S_cd,b,c) << "\t" << ward << endl;

                  }

               //3) a = d
               for(int S_ab = 0;S_ab < 2;++S_ab)
                  for(int S_cd = 0;S_cd < 2;++S_cd){

                     double ward = 0.0;

                     for(int S_lb = 0;S_lb < 2;++S_lb)
                        for(int S_cl = 0;S_cl < 2;++S_cl){

                           double hard = std::sqrt( (2.0*S_lb + 1.0) * (2.0*S_cl + 1.0) * (2.0*S_ab + 1.0) * (2.0*S_cd + 1.0) ) * (1 - 2*S_cl)

                              * (1 - 2*S_cd) * _6j[S_ab][S_lb] * _6j[S_cd][S_cl] * (*this)(a,0,S_lb,l,b,S_cl,c,l);

                           if(a == b)
                              hard /= std::sqrt(2.0);

                           if(c == a)
                              hard /= std::sqrt(2.0);

                           if(l == b)
                              hard *= std::sqrt(2.0);

                           if(c == l)
                              hard *= std::sqrt(2.0);

                           ward += hard;

                        }

                     if(fabs(ward - (*this)(l,0,S_ab,a,b,S_cd,c,a)) > 1.0e-12)
                        cout << "ab;ca\t" << a << "\t" << b << "\t" << c << "\t(" << S_ab << ")\t(" << S_cd << ")\t" << (*this)(l,0,S_ab,a,b,S_cd,c,a) << "\t" << ward << endl;

                  }

               //4) a = c
               for(int S_ab = 0;S_ab < 2;++S_ab)
                  for(int S_cd = 0;S_cd < 2;++S_cd){

                     double ward = 0.0;

                     for(int S_lb = 0;S_lb < 2;++S_lb)
                        for(int S_ld = 0;S_ld < 2;++S_ld){

                           double hard = std::sqrt( (2.0*S_lb + 1.0) * (2.0*S_ld + 1.0) * (2.0*S_ab + 1.0) * (2.0*S_cd + 1.0) ) 

                              * _6j[S_ab][S_lb] * _6j[S_cd][S_ld] * (*this)(a,0,S_lb,l,b,S_ld,l,c);

                           if(a == b)
                              hard /= std::sqrt(2.0);

                           if(c == a)
                              hard /= std::sqrt(2.0);

                           if(l == b)
                              hard *= std::sqrt(2.0);

                           if(c == l)
                              hard *= std::sqrt(2.0);

                           ward += hard;

                        }

                     if(fabs(ward - (*this)(l,0,S_ab,a,b,S_cd,a,c))>1.0e-12)
                        cout << "ab;ac\t" << a << "\t" << b << "\t" << c << "\t(" << S_ab << ")\t(" << S_cd << ")\t" << (*this)(l,0,S_ab,a,b,S_cd,a,c) << "\t" << ward << endl;

                  }

               //)5 a = l
               for(int S_ab = 0;S_ab < 2;++S_ab)
                  for(int S_cd = 0;S_cd < 2;++S_cd){

                     double ward = 0.0;

                     for(int S_lb = 0;S_lb < 2;++S_lb)
                        ward += std::sqrt( (2.0*S_ab + 1.0) * (2.0*S_lb + 1.0) ) * _6j[S_lb][S_ab] * (*this)(l,0,S_lb,l,a,S_cd,b,c);

                     if(fabs(ward - (*this)(l,0,S_ab,l,a,S_cd,b,c)) > 1.0e-12)
                        cout << "lb;cd\t" << a << "\t" << b << "\t" << c << "\t(" << S_ab << ")\t(" << S_cd << ")\t" << (*this)(l,0,S_ab,l,a,S_cd,b,c) << "\t" << ward << endl;


                  }

               //6) b = l
               for(int S_ab = 0;S_ab < 2;++S_ab)
                  for(int S_cd = 0;S_cd < 2;++S_cd){

                     double ward = 0.0;

                     for(int S_al = 0;S_al < 2;++S_al)
                        ward += std::sqrt( (2.0*S_ab + 1.0) * (2.0*S_al + 1.0) ) * (1 - 2*S_ab) * (1 - 2*S_al) * _6j[S_al][S_ab] * (*this)(l,0,S_al,a,l,S_cd,b,c);

                     if(fabs(ward - (*this)(l,0,S_ab,a,l,S_cd,b,c)) > 1.0e-12)
                        cout << "al;cd\t" << a << "\t" << b << "\t" << c << "\t(" << S_ab << ")\t(" << S_cd << ")\t" << (*this)(l,0,S_ab,a,l,S_cd,b,c) << "\t" << ward <<endl;

                  }

               //7) c = l
               for(int S_ab = 0;S_ab < 2;++S_ab)
                  for(int S_cd = 0;S_cd < 2;++S_cd){

                     double ward = 0.0;

                     for(int S_ld = 0;S_ld < 2;++S_ld)
                        ward += std::sqrt( (2.0*S_cd + 1.0) * (2.0*S_ld + 1.0) ) * _6j[S_cd][S_ld] * (*this)(l,0,S_ab,a,b,S_ld,l,c);

                     if(fabs(ward - (*this)(l,0,S_ab,a,b,S_cd,l,c)) > 1.0e-12)
                        cout << "ab;ld\t" << a << "\t" << b << "\t" << c << "\t(" << S_ab << ")\t(" << S_cd << ")\t" << (*this)(l,0,S_ab,a,b,S_cd,l,c) << "\t" << ward << endl;

                  }

               //8) d = l
               for(int S_ab = 0;S_ab < 2;++S_ab)
                  for(int S_cd = 0;S_cd < 2;++S_cd){

                     double ward = 0.0;

                     for(int S_cl = 0;S_cl < 2;++S_cl)
                        ward += std::sqrt( (2.0*S_cd + 1.0) * (2.0*S_cl + 1.0) ) * (1 - 2*S_cd) * (1 - 2*S_cl) * _6j[S_cd][S_cl] * (*this)(l,0,S_ab,a,b,S_cl,c,l);


                     if(fabs(ward - (*this)(l,0,S_ab,a,b,S_cd,c,l)) > 1.0e-12)
                        cout << "ab;cl\t" << a << "\t" << b << "\t" << c << "\t(" << S_ab << ")\t(" << S_cd << ")\t" << (*this)(l,0,S_ab,a,b,S_cd,c,l) << "\t" << ward << endl;


                  }

            }

      cout << endl;


      cout << endl;
      cout << "S = 3/2" << endl;
      cout << endl;

      for(int a = 0;a < M;++a)
         for(int b = 0;b < M;++b)
            for(int c = 0;c < M;++c){

               if(fabs((*this)(l,1,1,a,b,1,c,b) - (*this)(b,1,1,a,l,1,c,l)) > 1.0e-12)
                  cout << a << "\t" << b << "\t" << c << "\t|\t" << (*this)(l,1,1,a,b,1,c,b) << "\t" << (*this)(b,1,1,a,l,1,c,l) << endl;

               if(fabs((*this)(l,1,1,a,b,1,b,c) - (*this)(b,1,1,a,l,1,l,c)) > 1.0e-12)
                  cout << a << "\t" << b << "\t" << c << "\t|\t" << (*this)(l,1,1,a,b,1,b,c) << "\t" << (*this)(b,1,1,a,l,1,l,c) << endl;

               if(fabs((*this)(l,1,1,a,b,1,c,a) - (*this)(a,1,1,l,b,1,c,l)) > 1.0e-12)
                  cout << a << "\t" << b << "\t" << c << "\t|\t" << (*this)(l,1,1,a,b,1,c,a) << "\t" << (*this)(a,1,1,l,b,1,c,l) << endl;

               if(fabs((*this)(l,1,1,a,b,1,a,c) - (*this)(a,1,1,l,b,1,l,c)) > 1.0e-12)
                  cout << a << "\t" << b << "\t" << c << "\t|\t" << (*this)(l,1,1,a,b,1,a,c) << "\t" << (*this)(a,1,1,l,b,1,l,c) << endl;

            }

   }

}

/**
 * lift a TPM tpm object up to dDPM space (*this), so that the dotproduct of (*this) with a dDPM matrix is equal 
 * to the dotproduct of tpm with the "barred" dDPM matrix.
 */
void dDPM::up(const TPM &tpm){

   int a,b,c,d;
   int S_ab,S_cd;

   for(int l = 0;l < M;++l){

      for(int S = 0;S < 2;++S){

         for(int i = 0;i < ddpm[l]->gdim(S);++i){

            S_ab = rxTPM::gt2s(l,S,i,0);

            a = rxTPM::gt2s(l,S,i,1);
            b = rxTPM::gt2s(l,S,i,2);

            for(int j = i;j < ddpm[l]->gdim(S);++j){

               S_cd = rxTPM::gt2s(l,S,j,0);

               c = rxTPM::gt2s(l,S,j,1);
               d = rxTPM::gt2s(l,S,j,2);

               if(S_ab == S_cd)
                  (*this)[l](S,i,j) = tpm(S_ab,a,b,c,d)/(N - 2.0);
               else
                  (*this)[l](S,i,j) = 0.0;

            }
         }
      }

   }

   this->symmetrize();

}

/**
 * construct the spinsymmetrical hubbard hamiltonian with on site repulsion U
 * @param U onsite repulsion term
 */
void dDPM::hubbard(double U){

   TPM ham;
   ham.hubbard(U);

   this->up(ham);

}

/**
 * initialize (*this) on the correctly normalized unitmatrix (so that the trace is N(N-1)(N-2)/2)
 */
void dDPM::unit(){

   int S_ab,S_cd;
   int a,b,c,d;

   double norm;

   for(int l = 0;l < M;++l){

      for(int S = 0;S < 2;++S){

         for(int i = 0;i < ddpm[l]->gdim(S);++i){

            S_ab = rxTPM::gt2s(l,S,i,0);

            a = rxTPM::gt2s(l,S,i,1);
            b = rxTPM::gt2s(l,S,i,2);

            for(int j = i;j < ddpm[l]->gdim(S);++j){

               S_cd = rxTPM::gt2s(l,S,j,0);

               c = rxTPM::gt2s(l,S,j,1);
               d = rxTPM::gt2s(l,S,j,2);

               //set the norm
               norm = 1.0;

               if(a == b)
                  norm /= std::sqrt(2.0);

               if(c == d)
                  norm /= std::sqrt(2.0);

               (*this)[l](S,i,j) = 0.0;

               //set the unitmatrix
               if(a == c && b == d){

                  if(S_ab == S_cd)
                     (*this)[l](S,i,j) += 1.0;

                  if(a == l)
                     (*this)[l](S,i,j) += std::sqrt( (2*S_ab + 1.0) * (2.0*S_cd + 1.0) ) * _6j[S_ab][S_cd];

                  if(b == l)
                     (*this)[l](S,i,j) += std::sqrt( (2*S_ab + 1.0) * (2.0*S_cd + 1.0) ) * _6j[S_ab][S_cd] * (1 - 2*S_ab) * (1 - 2*S_cd);

               }

               if(a == d && b == c){

                  if(S_ab == S_cd)
                     (*this)[l](S,i,j) += (1 - 2*S_ab);

                  if(a == l)
                     (*this)[l](S,i,j) += std::sqrt( (2*S_ab + 1.0) * (2.0*S_cd + 1.0) ) * _6j[S_ab][S_cd] * (1 - 2*S_cd);

                  if(b == l)
                     (*this)[l](S,i,j) += std::sqrt( (2*S_ab + 1.0) * (2.0*S_cd + 1.0) ) * _6j[S_ab][S_cd] * (1 - 2*S_ab);

               }

               (*this)[l](S,i,j) *= norm * (N*(N - 1.0)*(N - 2.0)/(2*M*(2*M - 1.0)*(2*M - 2.0)));

               (*this)[l](S,j,i) = (*this)[l](S,i,j);

            }

         }

      }

   }

}

/**
 * set the matrix equal to the I1 part of a u^0 matrix.
 */
void dDPM::set_u_0(){

   *this = Tools::gunit();

   double ward = this->trace();

   this->dscal(N*(N - 1.0)*(N - 2.0)/(2.0*ward));

}

/**
 * The spincoupled Q2 map: maps a dDPM object onto itself
 * @param option if == 'U' up, if == 'D' down
 * @param ddpm_i input TPM
 */
void dDPM::Q(char option,const dDPM &ddpm_i){

   if(option == 'U'){

      TPM tpm;
      tpm.bar(1.0/(N - 2.0),ddpm_i);

      SPM spm(1.0/(N - 1.0),tpm);

      double ward = 2.0 * ddpm_i.trace()/( N*(N - 1.0)*(N - 2.0) );

      int a,b,c,d;
      int S_ab,S_cd;

      int sign_ab,sign_cd;

      double norm_ab,norm_cd;

      double hard;

      for(int l = 0;l < M;++l){

         //start with the S = 1/2 block, this is the most difficult one:
         for(int i = 0;i < ddpm[l]->gdim(0);++i){

            S_ab = rxTPM::gt2s(l,0,i,0);

            a = rxTPM::gt2s(l,0,i,1);
            b = rxTPM::gt2s(l,0,i,2);

            sign_ab = 1 - 2*S_ab;

            norm_ab = 1.0;

            if(a == b)
               norm_ab /= std::sqrt(2.0);

            for(int j = i;j < ddpm[l]->gdim(0);++j){

               S_cd = rxTPM::gt2s(l,0,j,0);

               c = rxTPM::gt2s(l,0,j,1);
               d = rxTPM::gt2s(l,0,j,2);

               sign_cd = 1 - 2*S_cd;

               norm_cd = 1.0;

               if(c == d)
                  norm_cd /= std::sqrt(2.0);

               hard = std::sqrt( (2*S_ab + 1.0) * (2*S_cd + 1.0) ) * _6j[S_ab][S_cd];

               //dp part
               (*this)[l](0,i,j) = -ddpm_i[l](0,i,j);

               //np(1)
               if(i == j)
                  (*this)[l](0,i,j) += ward - spm(l,l);

               //terms that contribute when the spin is diagonal:
               if(S_ab == S_cd){

                  //tp(1)
                  (*this)[l](0,i,j) += tpm(S_ab,a,b,c,d);

                  //sp(1) first term
                  if(b == d)
                     (*this)[l](0,i,j) -= norm_ab * norm_cd * spm(a,c);

                  //sp(2) first term
                  if(a == d)
                     (*this)[l](0,i,j) -= sign_ab * norm_ab * norm_cd * spm(b,c);

                  //sp(4) first term
                  if(b == c)
                     (*this)[l](0,i,j) -= sign_cd * norm_ab * norm_cd * spm(a,d);

                  //sp(5) first term
                  if(a == c)
                     (*this)[l](0,i,j) -= norm_ab * norm_cd * spm(b,d);

               }

               if(b == l){

                  //tp(2)
                  if(a == l)
                     (*this)[l](0,i,j) += std::sqrt(2.0) * norm_ab * sign_ab * sign_cd * hard * tpm(S_cd,a,l,c,d);
                  else
                     (*this)[l](0,i,j) += norm_ab * sign_ab * sign_cd * hard * tpm(S_cd,a,l,c,d);

                  //sp(1) second term
                  if(d == l)
                     (*this)[l](0,i,j) -= sign_ab * sign_cd * norm_ab * norm_cd * hard * spm(a,c);

                  //sp(3)
                  if(a == d)
                     (*this)[l](0,i,j) -= sign_ab * norm_ab * norm_cd * hard * spm(c,l);

                  //sp(4) second term
                  if(l == c)
                     (*this)[l](0,i,j) -= sign_ab * norm_ab * norm_cd * hard * spm(a,d);

                  //sp(6)
                  if(a == c)
                     (*this)[l](0,i,j) -= sign_ab * sign_cd * norm_ab * norm_cd * hard * spm(l,d);

                  //np(4)
                  if(c == l && a == d)
                     (*this)[l](0,i,j) += sign_ab * norm_ab * norm_cd * hard * ward;

                  //np(6)
                  if(d == l && a == c)
                     (*this)[l](0,i,j) += sign_ab * sign_cd * norm_ab * norm_cd * hard * ward;

               }

               if(a == l){

                  //tp(3)
                  if(b == l)
                     (*this)[l](0,i,j) += std::sqrt(2.0) * norm_ab * hard * tpm(S_cd,l,b,c,d);
                  else
                     (*this)[l](0,i,j) += norm_ab * hard * tpm(S_cd,l,b,c,d);

                  //sp(2) second term
                  if(d == l)
                     (*this)[l](0,i,j) -= sign_cd * norm_ab * norm_cd * hard * spm(b,c);

                  //sp(5) second term
                  if(c == l)
                     (*this)[l](0,i,j) -= norm_ab * norm_cd * hard * spm(b,d);

                  //sp(3) second part
                  if(b == d)
                     (*this)[l](0,i,j) -= norm_ab * norm_cd * hard * spm(c,l);

                  //sp(6) second part
                  if(b == c)
                     (*this)[l](0,i,j) -= sign_cd * norm_ab * norm_cd * hard * spm(d,l);

                  //np(3)
                  if(c == l && b == d)
                     (*this)[l](0,i,j) += norm_ab * norm_cd * hard * ward;

                  //np(5)
                  if(d == l && b == c)
                     (*this)[l](0,i,j) += sign_cd * norm_ab * norm_cd * hard * ward;

               }

               if(l == d){

                  //tp(4)
                  if(c == l)
                     (*this)[l](0,i,j) += std::sqrt(2.0) * norm_cd * sign_ab * sign_cd * hard * tpm(S_ab,a,b,c,l);
                  else
                     (*this)[l](0,i,j) += norm_cd * sign_ab * sign_cd * hard * tpm(S_ab,a,b,c,l);

                  //sp(7) first term
                  if(b == c)
                     (*this)[l](0,i,j) -= norm_ab * norm_cd * sign_cd * hard * spm(a,l);

                  //sp(8) first term
                  if(a == c)
                     (*this)[l](0,i,j) -= norm_ab * norm_cd * sign_ab * sign_cd * hard * spm(b,l);

               }

               if(b == d){

                  //tp(5)
                  double hulp = 0.0;

                  //sum over intermediate spin
                  for(int Z = 0;Z < 2;++Z)
                     hulp += (2*Z + 1.0) * _6j[Z][S_ab] * _6j[Z][S_cd] * tpm(Z,a,l,c,l);

                  //correct for norms of the tpm
                  if(a == l)
                     hulp *= std::sqrt(2.0);

                  if(c == l)
                     hulp *= std::sqrt(2.0);

                  (*this)[l](0,i,j) += norm_ab * norm_cd * sign_ab * sign_cd * std::sqrt( (2*S_ab + 1.0) * (2*S_cd + 1.0) ) * hulp;

                  //sp(7) second term
                  if(c == l)
                     (*this)[l](0,i,j) -= norm_ab * norm_cd * hard * spm(a,l);

               }

               if(a == d){

                  //tp(6)
                  double hulp = 0.0;

                  //sum over intermediate spin
                  for(int Z = 0;Z < 2;++Z)
                     hulp += (2*Z + 1.0) * _6j[Z][S_ab] * _6j[Z][S_cd] * tpm(Z,b,l,c,l);

                  if(b == l)
                     hulp *= std::sqrt(2.0);

                  if(c == l)
                     hulp *= std::sqrt(2.0);

                  (*this)[l](0,i,j) += sign_cd * std::sqrt( (2*S_ab + 1) * (2*S_cd + 1.0) ) * norm_ab * norm_cd * hulp;

                  //sp(8) second term
                  if(c == l)
                     (*this)[l](0,i,j) -= sign_ab * norm_ab * norm_cd * hard * spm(b,l);

               }

               if(c == l){

                  //tp(7)
                  if(d == l)
                     (*this)[l](0,i,j) += std::sqrt(2.0) * norm_cd * hard * tpm(S_ab,a,b,l,d);
                  else
                     (*this)[l](0,i,j) += norm_cd * hard * tpm(S_ab,a,b,l,d);

               }

               if(b == c){

                  //tp(8)
                  double hulp = 0.0;

                  //sum over intermediate spin
                  for(int Z = 0;Z < 2;++Z)
                     hulp += (2*Z + 1.0) * _6j[Z][S_ab] * _6j[Z][S_cd] * tpm(Z,a,l,d,l);

                  if(a == l)
                     hulp *= std::sqrt(2.0);

                  if(d == l)
                     hulp *= std::sqrt(2.0);

                  (*this)[l](0,i,j) += sign_ab * std::sqrt( (2*S_ab + 1) * (2*S_cd + 1.0) ) * norm_ab * norm_cd * hulp;

               }

               if(a == c){

                  //tp(8)
                  double hulp = 0.0;

                  //sum over intermediate spin
                  for(int Z = 0;Z < 2;++Z)
                     hulp += (2*Z + 1.0) * _6j[Z][S_ab] * _6j[Z][S_cd] * tpm(Z,b,l,d,l);

                  if(b == l)
                     hulp *= std::sqrt(2.0);

                  if(d == l)
                     hulp *= std::sqrt(2.0);

                  (*this)[l](0,i,j) += std::sqrt( (2*S_ab + 1) * (2*S_cd + 1.0) ) * norm_ab * norm_cd * hulp;

               }

            }
         }

         //then the S = 3/2 block, this should be easy, totally antisymmetrical 
         for(int i = 0;i < ddpm[l]->gdim(1);++i){

            a = rxTPM::gt2s(l,1,i,1);
            b = rxTPM::gt2s(l,1,i,2);

            for(int j = i;j < ddpm[l]->gdim(1);++j){

               c = rxTPM::gt2s(l,1,j,1);
               d = rxTPM::gt2s(l,1,j,2);


               (*this)[l](1,i,j) = tpm(1,a,b,c,d) - ddpm_i[l](1,i,j);

               if(i == j)
                  (*this)[l](1,i,j) += ward - spm(l,l);

               if(b == d)
                  (*this)[l](1,i,j) += tpm(1,a,l,c,l) - spm(a,c);

               if(b == c)
                  (*this)[l](1,i,j) -= tpm(1,a,l,d,l) - spm(a,d);

               if(a == c)
                  (*this)[l](1,i,j) += tpm(1,b,l,d,l) - spm(b,d);

            }
         }

      }

   }
   else{


      TPM tpm;
      tpm.bar(1.0/(N - 2.0),ddpm_i);

      SPM spm(1.0/(N - 1.0),tpm);

      dSPM dspm;
      dspm.trace(1/((N - 1.0)*(N - 2.0)),ddpm_i);

      SPM breve;
      breve.breve(1.0/((N - 1.0)*(N - 2.0)),ddpm_i);

      ssdTPM ssdtpm;
      ssdtpm.bar(1.0/((N - 1.0)*(N - 2.0)),ddpm_i);

      PHM phm;
      phm.spinsum(1.0/(N - 2.0),ddpm_i);

      dTPM dtpm;
      dtpm.bar(1.0/(N - 2.0),ddpm_i);

      double ward = 2.0 * ddpm_i.dotunit()/( N*(N - 1.0)*(N - 2.0) );

      int a,b,c,d;
      int S_ab,S_cd;

      int sign_ab,sign_cd;

      double norm_ab,norm_cd;

      for(int l = 0;l < M;++l){

         //start with the S = 1/2 block, this is the most difficult one:
         for(int i = 0;i < ddpm[l]->gdim(0);++i){

            S_ab = rxTPM::gt2s(l,0,i,0);

            a = rxTPM::gt2s(l,0,i,1);
            b = rxTPM::gt2s(l,0,i,2);

            sign_ab = 1 - 2*S_ab;

            norm_ab = 1.0;

            if(a == b)
               norm_ab /= std::sqrt(2.0);

            for(int j = i;j < ddpm[l]->gdim(0);++j){

               S_cd = rxTPM::gt2s(l,0,j,0);

               c = rxTPM::gt2s(l,0,j,1);
               d = rxTPM::gt2s(l,0,j,2);

               sign_cd = 1 - 2*S_cd;

               norm_cd = 1.0;

               if(c == d)
                  norm_cd /= std::sqrt(2.0);

               //dp part
               (*this)[l](0,i,j) = -ddpm_i[l](0,i,j);

               if(i == j)
                  (*this)[l](0,i,j) += ward - 0.5 * (dspm[a] + dspm[b]);

               if(S_ab == S_cd){

                  (*this)[l](0,i,j) += tpm(S_ab,a,b,c,d) + phm(S_ab,a,b,c,d) + (1 - 2*S_ab)*phm(S_ab,b,a,c,d) + (1 - 2*S_cd) * phm(S_cd,d,c,a,b) + phm(S_cd,c,d,a,b);

                  if(b == d)
                     (*this)[l](0,i,j) -= norm_ab * norm_cd * ( spm(a,c) + breve(a,c) + ssdtpm[a](a,c) + ssdtpm[c](a,c) - dtpm[b](S_ab,a,c) );

                  if(a == d)
                     (*this)[l](0,i,j) -= sign_ab * norm_ab * norm_cd * ( spm(b,c) + breve(b,c) + ssdtpm[b](b,c) + ssdtpm[c](b,c) - dtpm[a](S_ab,b,c) );

                  if(b == c)
                     (*this)[l](0,i,j) -= sign_cd * norm_ab * norm_cd * ( spm(a,d) + breve(a,d) + ssdtpm[a](a,d) + ssdtpm[d](a,d) - dtpm[b](S_ab,a,d));

                  if(a == c)
                     (*this)[l](0,i,j) -= norm_ab * norm_cd * ( spm(b,d) + breve(b,d) + ssdtpm[b](b,d) + ssdtpm[d](b,d) - dtpm[a](S_ab,b,d) );

               }

            }
         }

         //then the S = 3/2 block, this should be easy, totally antisymmetrical 
         for(int i = 0;i < ddpm[l]->gdim(1);++i){

            a = rxTPM::gt2s(l,1,i,1);
            b = rxTPM::gt2s(l,1,i,2);

            for(int j = i;j < ddpm[l]->gdim(1);++j){

               c = rxTPM::gt2s(l,1,j,1);
               d = rxTPM::gt2s(l,1,j,2);

               (*this)[l](1,i,j) = tpm(1,a,b,c,d) - ddpm_i[l](1,i,j) + phm(1,a,b,c,d) - phm(1,b,a,c,d) - phm(1,d,c,a,b) + phm(1,c,d,a,b);

               if(i == j)
                  (*this)[l](1,i,j) += ward - 0.5 * (dspm[a] + dspm[b]);

               if(b == d)
                  (*this)[l](1,i,j) -= spm(a,c) + breve(a,c) + ssdtpm[a](a,c) + ssdtpm[c](a,c) - dtpm[b](1,a,c);

               if(b == c)
                  (*this)[l](1,i,j) += spm(a,d) + breve(a,d) + ssdtpm[a](a,d) + ssdtpm[d](a,d) - dtpm[b](1,a,d);

               if(a == c)
                  (*this)[l](1,i,j) -= spm(b,d) + breve(b,d) + ssdtpm[b](b,d) + ssdtpm[d](b,d) - dtpm[a](1,b,d);

            }
         }

      }

   }

   this->symmetrize();

}

/**
 * @return the trace of *this with the unitmatrix in dDP space.
 */
double dDPM::dotunit() const{

   double ward = this->trace();

   for(int S_ab = 0;S_ab < 2;++S_ab)
      for(int S_cd = 0;S_cd < 2;++S_cd){

         double hard = 0.0;

         for(int l = 0;l < M;++l)
            for(int b = 0;b < M;++b){

               if(l == b)
                  hard += 2.0 * (*this)(l,0,S_ab,l,b,S_cd,l,b);
               else
                  hard += (*this)(l,0,S_ab,l,b,S_cd,l,b);

            }

         ward += 2.0 * hard * std::sqrt( (2.0*S_ab + 1.0) * (2.0*S_cd + 1.0) ) * _6j[S_ab][S_cd];

      }

   return ward;

}

/**
 * map a dPPHM on a dDPM object with the I2 down map.
 * @param dpphm input dPPHM matrix
 */
void dDPM::I(const dPPHM &dpphm){

   int a,b,c,d;

   int S_ab,S_cd;

   TPM tpm;
   tpm.bar(1.0/(N - 2.0),dpphm);

   for(int l = 0;l < M;++l){

      for(int S = 0;S < 2;++S){

         for(int i = 0;i < ddpm[l]->gdim(S);++i){

            S_ab = rxTPM::gt2s(l,S,i,0);

            a = rxTPM::gt2s(l,S,i,1);
            b = rxTPM::gt2s(l,S,i,2);

            for(int j = i;j < ddpm[l]->gdim(S);++j){

               S_cd = rxTPM::gt2s(l,S,j,0);

               c = rxTPM::gt2s(l,S,j,1);
               d = rxTPM::gt2s(l,S,j,2);

               (*this)[l](S,i,j) = 0.0;

               for(int S_ = 0;S_ < 2;++S_)
                  (*this)[l](S,i,j) -= (2* (S_ + 0.5) + 1.0) * Tools::gC(S,S_,S_ab,S_cd) * dpphm(l,S_,S_ab,a,b,S_cd,c,d);

               if(S_ab == S_cd)
                  (*this)[l](S,i,j) += tpm(S_ab,a,b,c,d);

            }
         }
      }

   }

   this->symmetrize();

}

/**
 * map a dPPHM on a dDPM object with the Q1 down map.
 * @param dpphm input dPPHM matrix
 */
void dDPM::Q(const dPPHM &dpphm){

   int a,b,c,d;

   int S_ab,S_cd;

   TPM tpm;
   tpm.bar(1.0/(N - 2.0),dpphm);

   SPM breve;
   breve.breve(0.5/((N - 1.0)*(N - 2.0)),dpphm);

   dSPM dspm;
   dspm.trace(1.0/((N - 1.0)*(N - 2.0)),dpphm);

   ssdTPM ssdtpm;
   ssdtpm.bar(0.5/((N - 1.0)*(N - 2.0)),dpphm);

   PHM phm;
   phm.spinsum(1.0/(N - 2.0),dpphm);

   dTPM dtpm;
   dtpm.bar(1.0/(N - 2.0),dpphm);

   double ward = 2.0 * dpphm.barbreve()/( N*(N - 1.0)*(N - 2.0) );

   double norm_ab,norm_cd;
   int sign_ab,sign_cd;

   for(int l = 0;l < M;++l){

      for(int S = 0;S < 2;++S){

         for(int i = 0;i < ddpm[l]->gdim(S);++i){

            S_ab = rxTPM::gt2s(l,S,i,0);

            a = rxTPM::gt2s(l,S,i,1);
            b = rxTPM::gt2s(l,S,i,2);

            norm_ab = 1.0;

            sign_ab = 1 - 2*S_ab;

            if(a == b)
               norm_ab /= std::sqrt(2.0);

            for(int j = i;j < ddpm[l]->gdim(S);++j){

               S_cd = rxTPM::gt2s(l,S,j,0);

               c = rxTPM::gt2s(l,S,j,1);
               d = rxTPM::gt2s(l,S,j,2);

               norm_cd = 1.0;

               sign_cd = 1 - 2*S_cd;

               if(c == d)
                  norm_cd /= std::sqrt(2.0);

               (*this)[l](S,i,j) = 0.0;

               for(int S_ = 0;S_ < 2;++S_)
                  (*this)[l](S,i,j) += (2* (S_ + 0.5) + 1.0) * Tools::gC(S,S_,S_ab,S_cd) * dpphm(l,S_,S_ab,a,b,S_cd,c,d);

               if(i == j)
                  (*this)[l](S,i,j) += ward + 0.5 * ( dspm[a] + dspm[b] );

               if(S_ab == S_cd){

                  (*this)[l](S,i,j) += phm(S_cd,a,b,c,d) + sign_cd * phm(S_cd,b,a,c,d) + sign_ab * phm(S_ab,d,c,a,b) + phm(S_ab,c,d,a,b) ;

                  if(a == c)
                     (*this)[l](S,i,j) -= norm_ab * norm_cd * ( breve(b,d) + ssdtpm[b](b,d) + ssdtpm[d](b,d) + dtpm[a](S_ab,b,d) );

                  if(b == c)
                     (*this)[l](S,i,j) -= norm_ab * norm_cd * sign_ab * ( breve(a,d) + ssdtpm[a](a,d) + ssdtpm[d](a,d) + dtpm[b](S_ab,a,d) );

                  if(a == d)
                     (*this)[l](S,i,j) -= norm_ab * norm_cd * sign_cd * ( breve(b,c) + ssdtpm[b](b,c) + ssdtpm[c](b,c) + dtpm[a](S_ab,b,c) );

                  if(b == d)
                     (*this)[l](S,i,j) -= norm_ab * norm_cd * ( breve(a,c) + ssdtpm[a](a,c) + ssdtpm[c](a,c) + dtpm[b](S_ab,a,c) );

               }

            }
         }
      }

   }

   this->symmetrize();

}

/**
 * map a dPHHM on a dDPM object with the G1 down map.
 * @param dphhm input dPHHM matrix
 */
void dDPM::G1(const dPHHM &dphhm){

   int a,b,c,d;

   int S_ab,S_cd;

   double norm_ab,norm_cd;
   int sign_ab,sign_cd;

   dTPM bar;
   bar.bar(1.0/(N - 2.0),dphhm);

   SPM breve;
   breve.breve(0.5/((N - 1.0)*(N - 2.0)),dphhm);

   ssdTPM ssdtpm;
   ssdtpm.skew_bar(0.5/((N - 1.0)*(N - 2.0)),dphhm);

   dSPM dspm;
   dspm.skew_trace(0.5/((N - 1.0)*(N - 2.0)),dphhm);

   dTPM skew_bar;
   skew_bar.skew_bar(1.0/(N - 2.0),dphhm);

   PHM phm;
   phm.spinsum(1.0/(N - 2.0),dphhm);

   for(int l = 0;l < M;++l){

      for(int S = 0;S < 2;++S){

         for(int i = 0;i < ddpm[l]->gdim(S);++i){

            S_ab = rxTPM::gt2s(l,S,i,0);

            a = rxTPM::gt2s(l,S,i,1);
            b = rxTPM::gt2s(l,S,i,2);

            norm_ab = 1.0;

            sign_ab = 1 - 2*S_ab;

            if(a == b)
               norm_ab /= std::sqrt(2.0);

            for(int j = i;j < ddpm[l]->gdim(S);++j){

               S_cd = rxTPM::gt2s(l,S,j,0);

               c = rxTPM::gt2s(l,S,j,1);
               d = rxTPM::gt2s(l,S,j,2);

               norm_cd = 1.0;

               sign_cd = 1 - 2*S_cd;

               if(c == d)
                  norm_cd /= std::sqrt(2.0);

               (*this)[l](S,i,j) = 0.0;

               for(int S_ = 0;S_ < 2;S_++)
                  for(int S_bl = 0;S_bl < 2;++S_bl)
                     for(int S_dl = 0;S_dl < 2;++S_dl){

                        (*this)[l](S,i,j) += norm_ab * norm_cd * (2*(S_ + 0.5) + 1.0) 
                        
                           * std::sqrt( (2.0*S_ab + 1.0) * (2.0*S_cd + 1.0) * (2.0*S_bl + 1.0) * (2.0*S_dl + 1.0) )

                           * Tools::g6j(2*S + 1,1,2*S_dl,1,1,2*S_ab) * Tools::g6j(2*S + 1,1,2*S_bl,1,1,2*S_cd) * Tools::g6j(2*S + 1,2*S_bl,1,2*S_ + 1,2*S_dl,1)

                           * ( dphhm(l,S_,S_bl,a,d,S_dl,c,b) + sign_ab * dphhm(l,S_,S_bl,b,d,S_dl,c,a) + sign_cd * dphhm(l,S_,S_bl,a,c,S_dl,d,b)
                           
                                 + sign_ab * sign_cd * dphhm(l,S_,S_bl,b,c,S_dl,d,a) );

                     }


               if(i == j)
                  (*this)[l](S,i,j) += dspm[a] + dspm[b];

               if(S_ab == S_cd){

                  (*this)[l](S,i,j) -= norm_ab * norm_cd * ( phm(S_ab,a,d,c,b) + phm(S_ab,c,b,a,d) + phm(S_ab,b,c,d,a) + phm(S_ab,d,a,b,c)

                        + sign_ab * ( phm(S_ab,a,c,d,b) + phm(S_ab,d,b,a,c) + phm(S_ab,b,d,c,a) + phm(S_ab,c,a,b,d) ) );

                  if(a == c){

                     (*this)[l](S,i,j) += norm_ab * norm_cd * ( bar[a](S_ab,b,d) + breve(b,d) + ssdtpm[b](b,d) + ssdtpm[d](d,b)

                           - skew_bar[a](S_ab,b,d) - skew_bar[a](S_ab,d,b) );

                  }

                  if(b == c){

                     (*this)[l](S,i,j) += sign_ab * norm_ab * norm_cd * ( bar[b](S_ab,a,d) + breve(a,d) + ssdtpm[a](a,d) + ssdtpm[d](d,a)

                           - skew_bar[b](S_ab,a,d) - skew_bar[b](S_ab,d,a) );

                  }

                  if(a == d){

                     (*this)[l](S,i,j) += sign_ab * norm_ab * norm_cd * ( bar[a](S_ab,b,c) + breve(b,c) + ssdtpm[b](b,c) + ssdtpm[c](c,b)

                           - skew_bar[a](S_ab,b,c) - skew_bar[a](S_ab,c,b) );

                  }

                  if(b == d){

                     (*this)[l](S,i,j) += norm_ab * norm_cd * ( bar[b](S_ab,a,c) + breve(a,c) + ssdtpm[a](a,c) + ssdtpm[c](c,a)

                           - skew_bar[b](S_ab,a,c) - skew_bar[b](S_ab,c,a) );

                  }

               }

            }
         }
      }

   }

      this->symmetrize();

}


/**
 * map a dPHHM on a dDPM object with the G1 down map.
 * @param dphhm input dPHHM matrix
 */
void dDPM::G2(const dPHHM &dphhm){

   SPM barbar;
   barbar.barbar(0.5/((N - 1.0)*(N - 2.0)),dphhm);

   SPM breve_si;
   breve_si.breve_si(0.5/((N - 1.0)*(N - 2.0)),dphhm);

   PHM phm;
   phm.bar(1.0/(N - 2.0),dphhm);

   dTPM dtpm;
   dtpm.ssbar(1.0/(N - 2.0),dphhm);

   PHM hat_si;
   hat_si.spinsum_si(1.0/(N - 2.0),dphhm);

   int a,b,c,d;

   int S_ab,S_cd;

   double norm_ab,norm_cd;
   int sign_ab,sign_cd;

   for(int l = 0;l < M;++l){

      for(int S = 0;S < 2;++S){

         for(int i = 0;i < ddpm[l]->gdim(S);++i){

            S_ab = rxTPM::gt2s(l,S,i,0);

            a = rxTPM::gt2s(l,S,i,1);
            b = rxTPM::gt2s(l,S,i,2);

            norm_ab = 1.0;

            sign_ab = 1 - 2*S_ab;

            if(a == b)
               norm_ab /= std::sqrt(2.0);

            for(int j = i;j < ddpm[l]->gdim(S);++j){

               S_cd = rxTPM::gt2s(l,S,j,0);

               c = rxTPM::gt2s(l,S,j,1);
               d = rxTPM::gt2s(l,S,j,2);

               norm_cd = 1.0;

               sign_cd = 1 - 2*S_cd;

               if(c == d)
                  norm_cd /= std::sqrt(2.0);

               (*this)[l](S,i,j) = 0.0;

               for(int S_ = 0;S_ < 2;S_++)
                  for(int S_bl = 0;S_bl < 2;++S_bl)
                     for(int S_dl = 0;S_dl < 2;++S_dl){

                        (*this)[l](S,i,j) -= norm_ab * norm_cd * (2*(S_ + 0.5) + 1.0) 
                        
                           * std::sqrt( (2.0*S_ab + 1.0) * (2.0*S_cd + 1.0) * (2.0*S_bl + 1.0) * (2.0*S_dl + 1.0) )

                           * Tools::g6j(2*S + 1,1,2*S_dl,1,1,2*S_ab) * Tools::g6j(2*S + 1,1,2*S_bl,1,1,2*S_cd) * Tools::g6j(2*S + 1,2*S_bl,1,2*S_ + 1,2*S_dl,1)

                           * ( dphhm(l,S_,S_bl,a,d,S_dl,c,b) + sign_ab * dphhm(l,S_,S_bl,b,d,S_dl,c,a) + sign_cd * dphhm(l,S_,S_bl,a,c,S_dl,d,b)
                           
                                 + sign_ab * sign_cd * dphhm(l,S_,S_bl,b,c,S_dl,d,a) );

                     }

               if(S_ab == S_cd){

                  //tp_a 
                  (*this)[l](S,i,j) -= norm_ab * norm_cd * ( phm(S_ab,a,d,c,b) + sign_ab * phm(S_ab,a,c,d,b) + sign_ab * phm(S_ab,b,d,c,a) + phm(S_ab,b,c,d,a) );

                  //tp_b and tp_c
                  (*this)[l](S,i,j) -= norm_ab * norm_cd * ( hat_si(S_ab,a,d,c,b) + sign_ab * hat_si(S_ab,b,d,c,a) + sign_ab * hat_si(S_ab,a,c,d,b) + hat_si(S_ab,b,c,d,a)

                        + hat_si(S_ab,c,b,a,d) + sign_ab * hat_si(S_ab,c,a,b,d) + sign_ab * hat_si(S_ab,d,b,a,c) + hat_si(S_ab,d,a,b,c) );

                  if(b == d)
                     (*this)[l](S,i,j) += norm_ab * norm_cd * ( barbar(a,c) + breve_si(a,c) - dtpm[b](S_ab,a,c) );

                  if(a == d)
                     (*this)[l](S,i,j) += sign_ab * norm_ab * norm_cd * ( barbar(b,c) + breve_si(b,c) - dtpm[a](S_ab,b,c) );

                  if(b == c)
                     (*this)[l](S,i,j) += sign_ab * norm_ab * norm_cd * ( barbar(a,d) + breve_si(a,d) - dtpm[b](S_ab,a,d) );

                  if(a == c)
                     (*this)[l](S,i,j) += norm_ab * norm_cd * ( barbar(b,d) + breve_si(b,d) - dtpm[a](S_ab,b,d) );

               }


            }
         }
      }

   }

   this->symmetrize();

}

/**
 * Collaps a SUP matrix X onto a dDPM matrix like this:\n\n
 * sum_i Tr (X u^i)f^i = this
 * @param X input SUP
 */
void dDPM::collaps(const SUP &X){

   dDPM hulp(X.gI1());

   *this = hulp;

#ifdef __Q2_CON
   hulp.Q('D',X.gQ2());

   *this += hulp;
#endif

#ifdef __I2_CON
   hulp.I(X.gI2());

   *this += hulp;
#endif

#ifdef __Q1_CON
   hulp.Q(X.gQ1());

   *this += hulp;
#endif

#ifdef __G1_CON
   hulp.G1(X.gG1());

   *this += hulp;
#endif

#ifdef __G2_CON
   hulp.G2(X.gG2());

   *this += hulp;
#endif

   this->proj();

}

/**
 * Implementation of a linear conjugate gradient algoritm for the solution of the overlapmatrix-system
 * S(*this) =  b\n\n 
 * in which S represents the hessian map.
 * @param b righthandside of the equation
 * @return return number of iterations needed to converge to the desired accuracy
 */
int dDPM::solve(dDPM &b){

   *this = 0;

   //de r initialiseren op b
   dDPM r(b);

   double rr = r.ddot(r);
   double rr_old,ward;

   dDPM Sb;

   int cg_iter = 0;

   while(rr > 1.0e-10){

      ++cg_iter;

      Sb.S(b);

      ward = rr/b.ddot(Sb);

      //delta += ward*b
      this->daxpy(ward,b);

      //r -= ward*Hb
      r.daxpy(-ward,Sb);

      //nieuwe r_norm maken
      rr_old = rr;
      rr = r.ddot(r);

      //eerst herschalen van b
      b.dscal(rr/rr_old);

      //dan r er bijtellen
      b += r;

   }

   return cg_iter;

}

/**
 * map a dDPM object onto a dDPM object using the overlapmatrix-map
 * @param ddpm input dDPM object
 */
void dDPM::S(const dDPM &ddpm){

   *this = ddpm;

   dDPM down;
   dPPHM up_pph;
   dPHHM up_phh;

#ifdef __Q2_CON
   dDPM up;
   up.Q('U',ddpm);

   down.Q('D',up);

   *this += down;
#endif

#ifdef __I2_CON
   up_pph.I(ddpm);
   down.I(up_pph);

   *this += down;
#endif

#ifdef __Q1_CON
   up_pph.Q(ddpm);
   down.Q(up_pph);

   *this += down;
#endif

#ifdef __G1_CON
   up_phh.G1(ddpm);
   down.G1(up_phh);

   *this += down;
#endif

#ifdef __G2_CON
   up_phh.G2(ddpm);
   down.G2(up_phh);

   *this += down;
#endif

   this->proj();

}

/**
 * Seperate matrix into two matrices, a positive and negative semidefinite part.
 * @param p positive (plus) output part
 * @param m negative (minus) output part
 */
void dDPM::sep_pm(dDPM &p,dDPM &m){

   for(int l = 0;l < M;++l)
      ddpm[l]->sep_pm(p[l],m[l]);

}
