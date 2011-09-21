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

   int phase_i = get_inco(l,S,S_ab,a,b);

   if(phase_i == 0)
      return 0.0;

   int phase_j = get_inco(l,S,S_cd,c,d);

   if(phase_j == 0)
      return 0.0;

   return phase_i * phase_j * (*this)[l](S,rTPM::gs2t(l,S,S_ab,a,b),rTPM::gs2t(l,S,S_cd,c,d));

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

               i = rTPM::gs2t(l,0,S_ab,l,a);
               j = rTPM::gs2t(l,0,S_cd,l,a);

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

            i = rTPM::gs2t(a,0,0,l,l);

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

            i = rTPM::gs2t(l,0,0,a,a);
            j = rTPM::gs2t(l,0,S_cd,a,l);

            if(a > l)
               phase_j = 1 - 2*S_cd;
            else
               phase_j = 1;

            (*this)[l](0,i,j) = 0.25 * phase_j * ward;
            (*this)[l](0,j,i) = (*this)[l](0,i,j);

         }

         //then the "a" block
         for(int S_lb = 0;S_lb < 2;++S_lb){

            i = rTPM::gs2t(a,0,S_lb,l,a);
            j = rTPM::gs2t(a,0,0,l,l);

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

                  i = rTPM::gs2t(l,0,S_ab,a,b);
                  j = rTPM::gs2t(l,0,S_cd,a,b);

                  (*this)[l](0,i,j) = ward/3.0;
                  (*this)[l](0,j,i) = (*this)[l](0,i,j);

               }

            //then make the rest symmetric
            for(int S_lb = 0;S_lb < 2;++S_lb)
               for(int S_ld = S_lb;S_ld < 2;++S_ld){

                  i = rTPM::gs2t(a,0,S_lb,l,b);
                  j = rTPM::gs2t(a,0,S_ld,l,b);

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

                  i = rTPM::gs2t(b,0,S_al,a,l);
                  j = rTPM::gs2t(b,0,S_cl,a,l);

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

                  i = rTPM::gs2t(l,0,S_ab,a,l);
                  j = rTPM::gs2t(l,0,S_cd,c,a);

                  if(a > l)
                     phase_i = 1 - 2*S_ab;
                  else
                     phase_i = 1;

                  (*this)[l](0,i,j) = phase_i * ward / 3.0;
                  (*this)[l](0,j,i) = (*this)[l](0,i,j);

               }

            //then the "a" block
            for(int S_cl = 0;S_cl < 2;++S_cl){

               i = rTPM::gs2t(a,0,0,l,l);
               j = rTPM::gs2t(a,0,S_cl,c,l);

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

                  i = rTPM::gs2t(l,0,S_ab,a,l);
                  j = rTPM::gs2t(l,0,S_cd,a,d);

                  if(a > l)
                     phase_i = 1 - 2*S_ab;
                  else
                     phase_i = 1;

                  (*this)[l](0,i,j) = phase_i * ward / 3.0;
                  (*this)[l](0,j,i) = (*this)[l](0,i,j);

               }

            //then the "a" block
            for(int S_ld = 0;S_ld < 2;++S_ld){

               i = rTPM::gs2t(a,0,0,l,l);
               j = rTPM::gs2t(a,0,S_ld,l,d);

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

/*
      //next with two equalities: a = b = c: W^l_{aa;ad} <--> W^a_{la;ld}
      for(int a = 0;a < M;++a){

         if(a == l)
            ++a;

         if(a == M)
            break;

         //a = b = c: W^l_{aa;ad} <--> W^a_{la;ld}
         for(int c = 0;c < a;++c){

            if(c == l)
               ++c;

            if(c == a)
               break;

            //first take the average
            for(int S_cd = 0;S_cd < 2;++S_cd){

               double ward = (*this)(l,0,0,a,a,S_cd,c,a);

               for(int S_al = 0;S_al < 2;++S_al)
                  for(int S_cl = 0;S_cl < 2;++S_cl){

                     ward += std::sqrt( (2.0*S_cd + 1.0) * (2.0*S_al + 1.0) * (2.0*S_cl + 1.0) ) * (1 - 2*S_al) * (1 - 2*S_cl) * (1 - 2*S_cd)

                        * _6j[0][S_al] * _6j[S_cd][S_cl] * (*this)(a,0,S_al,a,l,S_cl,c,l);

                  }

               i = rTPM::gs2t(l,0,0,a,a);
               j = rTPM::gs2t(l,0,S_cd,c,a);

               (*this)[l](0,i,j) = 0.5*ward;
               (*this)[l](0,j,i) = (*this)[l](0,i,j);

            }

            //symmetrize the rest
            for(int S_al = 0;S_al < 2;++S_al)
               for(int S_cl = 0;S_cl < 2;++S_cl){

                  i = rTPM::gs2t(a,0,S_al,a,l);
                  j = rTPM::gs2t(a,0,S_cl,c,l);

                  if(a > l)
                     phase_i = 1 - 2*S_al;
                  else
                     phase_i = 1;

                  if(c > l)
                     phase_j = 1 - 2*S_cl;
                  else
                     phase_j = 1;

                  (*this)[a](0,i,j) = 0.0;

                  for(int S_cd = 0;S_cd < 2;++S_cd){

                     (*this)[a](0,i,j) += phase_i*phase_j* std::sqrt( (2.0*S_al + 1.0) * (2.0*S_cl + 1.0) * (2.0*S_cd + 1.0) ) * (1 - 2*S_cd) * (1 - 2*S_al)

                        * (1 - 2*S_cl) * _6j[0][S_al] * _6j[S_cd][S_cl] * (*this)(l,0,0,a,a,S_cd,c,a);

                  }

                  (*this)[a](0,j,i) = (*this)[a](0,i,j);

               }

         }

         for(int d = a + 1;d < M;++d){

            if(d == l)
               d++;
            
            if(d == M)
               break;

            //first take the average
            for(int S_cd = 0;S_cd < 2;++S_cd){

               double ward = (*this)(l,0,0,a,a,S_cd,a,d);

               for(int S_lb = 0;S_lb < 2;++S_lb)
                  for(int S_ld = 0;S_ld < 2;++S_ld){

                     ward += std::sqrt( (2.0*S_cd + 1.0) * (2.0*S_lb + 1.0) * (2.0*S_ld + 1.0) ) * 

                        _6j[0][S_lb] * _6j[S_cd][S_ld] * (*this)(a,0,S_lb,l,a,S_ld,l,d);

                  }

               i = rTPM::gs2t(l,0,0,a,a);
               j = rTPM::gs2t(l,0,S_cd,a,d);

               (*this)[l](0,i,j) = 0.5*ward;
               (*this)[l](0,j,i) = (*this)[l](0,i,j);

            }

            //symmetrize the rest
            for(int S_lb = 0;S_lb < 2;++S_lb)
               for(int S_ld = 0;S_ld < 2;++S_ld){

                  i = rTPM::gs2t(a,0,S_lb,l,a);
                  j = rTPM::gs2t(a,0,S_ld,l,d);

                  if(l > a)
                     phase_i = 1 - 2*S_lb;
                  else
                     phase_i = 1;

                  if(l > d)
                     phase_j = 1 - 2*S_ld;
                  else
                     phase_j = 1;

                  (*this)[a](0,i,j) = 0.0;

                  for(int S_cd = 0;S_cd < 2;++S_cd){

                     (*this)[a](0,i,j) += phase_i*phase_j* std::sqrt( (2.0*S_lb + 1.0) * (2.0*S_ld + 1.0) * (2.0*S_cd + 1.0) ) 

                        * _6j[0][S_lb] * _6j[S_cd][S_ld] * (*this)(l,0,0,a,a,S_cd,a,d);

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

                  i = rTPM::gs2t(l,0,S_ab,l,b);
                  j = rTPM::gs2t(l,0,S_cd,l,d);

                  double ward = mat[S_ab][S_cd];

                  for(int S_lb = 0;S_lb < 2;++S_lb)
                     for(int S_ld = 0;S_ld < 2;++S_ld){

                        ward += std::sqrt( (2.0*S_ab + 1.0) * (2.0*S_cd + 1.0) * (2.0*S_lb + 1.0) * (2.0*S_ld + 1.0) )

                           * _6j[S_ab][S_lb] * _6j[S_cd][S_ld] * mat[S_lb][S_ld];

                     }

                  (*this)[l](0,i,j) = 0.5 * ward;
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

                  i = rTPM::gs2t(l,0,S_ab,a,l);
                  j = rTPM::gs2t(l,0,S_cd,l,d);

                  double ward = mat[S_ab][S_cd];

                  for(int S_al = 0;S_al < 2;++S_al)
                     for(int S_ld = 0;S_ld < 2;++S_ld){

                        ward += std::sqrt( (2.0*S_ab + 1.0) * (2.0*S_cd + 1.0) * (2.0*S_al + 1.0) * (2.0*S_ld + 1.0) ) * (1 - 2*S_al) * (1 - 2*S_ab)

                           * _6j[S_ab][S_al] * _6j[S_cd][S_ld] * mat[S_al][S_ld];

                     }

                  (*this)[l](0,i,j) = 0.5*ward;
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

                  i = rTPM::gs2t(l,0,S_ab,a,l);
                  j = rTPM::gs2t(l,0,S_cd,c,l);

                  double ward = mat[S_ab][S_cd];
                  
                  for(int S_al = 0;S_al < 2;++S_al)
                     for(int S_cl = 0;S_cl < 2;++S_cl){

                        ward += std::sqrt( (2.0*S_ab + 1.0) * (2.0*S_cd + 1.0) * (2.0*S_al + 1.0) * (2.0*S_cl + 1.0) ) * (1 - 2*S_al) * (1 - 2*S_ab)

                           * (1 - 2*S_cd) * (1 - 2*S_cl) * _6j[S_ab][S_al] * _6j[S_cd][S_cl] * mat[S_al][S_cl];

                     }

                  (*this)[l](0,i,j) = 0.5*ward;
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

                     i = rTPM::gs2t(l,0,S_ab,a,b);
                     j = rTPM::gs2t(l,0,S_cd,a,d);

                     (*this)[l](0,i,j) = 0.5*ward;
                     (*this)[l](0,j,i) = (*this)[l](0,i,j);

                  }
                  
                  //then symmetrize the other term
                  for(int S_lb = 0;S_lb < 2;++S_lb)
                     for(int S_ld = 0;S_ld < 2;++S_ld){

                        i = rTPM::gs2t(a,0,S_lb,l,b);
                        j = rTPM::gs2t(a,0,S_ld,l,d);

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

                     i = rTPM::gs2t(l,0,S_ab,a,b);
                     j = rTPM::gs2t(l,0,S_cd,b,d);
                     
                     (*this)[l](0,i,j) = 0.5 * ward;
                     (*this)[l](0,j,i) = (*this)[l](0,i,j);

                  }

               //then symmetrize rest
               for(int S_al = 0;S_al < 2;++S_al)
                  for(int S_ld = 0;S_ld < 2;++S_ld){

                     i = rTPM::gs2t(b,0,S_al,a,l);
                     j = rTPM::gs2t(b,0,S_ld,l,d);

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

                     i = rTPM::gs2t(l,0,S_ab,a,b);
                     j = rTPM::gs2t(l,0,S_cd,c,b);

                     (*this)[l](0,i,j) = 0.5*ward;
                     (*this)[l](0,j,i) = (*this)[l](0,i,j);

                  }

               //then make the rest symmetric
               for(int S_al = 0;S_al < 2;++S_al)
                  for(int S_cl = 0;S_cl < 2;++S_cl){

                     i = rTPM::gs2t(b,0,S_al,a,l);
                     j = rTPM::gs2t(b,0,S_cl,c,l);

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

      //a = l: W^l_{lb;cd}
      for(int b = l + 1;b < M;++b)
         for(int c = l + 1;c < M;++c){

            if(c == b)
               b++;

            if(c == M)
               break;

            for(int d = c + 1;d < M;++d){

               if(d == b)
                  d++;

               if(d == M)
                  break;

               for(int S_ab = 0;S_ab < 2;++S_ab)
                  for(int S_cd = 0;S_cd < 2;++S_cd)
                     mat[S_ab][S_cd] = (*this)(l,0,S_ab,l,b,S_cd,c,d);

               for(int S_ab = 0;S_ab < 2;++S_ab)
                  for(int S_cd = 0;S_cd < 2;++S_cd){

                     double ward = mat[S_ab][S_cd];

                     for(int S_lb = 0;S_lb < 2;++S_lb)
                        ward += std::sqrt( (2.0*S_ab + 1.0) * (2.0*S_lb + 1.0) ) * _6j[S_ab][S_lb] * mat[S_lb][S_cd];

                     i = rTPM::gs2t(l,0,S_ab,l,b);
                     j = rTPM::gs2t(l,0,S_cd,c,d);

                     (*this)[l](0,i,j) = 0.5 * ward;
                     (*this)[l](0,j,i) = (*this)[l](0,i,j);

                  }

            }
         }

      //b = l: W^l_{al;cd}
      for(int a = 0;a < l;++a)
         for(int c = a + 1;c < M;++c){

            if(c == l)
               c++;

            if(c == M)
               break;

            for(int d = c + 1;d < M;++d){

               if(d == l)
                  d++;

               if(d == M)
                  break;

               for(int S_ab = 0;S_ab < 2;++S_ab)
                  for(int S_cd = 0;S_cd < 2;++S_cd)
                     mat[S_ab][S_cd] = (*this)(l,0,S_ab,a,l,S_cd,c,d);

               for(int S_ab = 0;S_ab < 2;++S_ab)
                  for(int S_cd = 0;S_cd < 2;++S_cd){

                     double ward = mat[S_ab][S_cd];

                     for(int S_al = 0;S_al < 2;++S_al)
                        ward += std::sqrt( (2.0*S_ab + 1.0) * (2.0*S_al + 1.0) ) * (1 - 2*S_ab) * (1 - 2*S_al) * _6j[S_ab][S_al] * mat[S_al][S_cd];

                     i = rTPM::gs2t(l,0,S_ab,a,l);
                     j = rTPM::gs2t(l,0,S_cd,c,d);

                     (*this)[l](0,i,j) = 0.5 * ward;
                     (*this)[l](0,j,i) = (*this)[l](0,i,j);

                  }

            }
         }

      //c = l: W^l_{ab;ld}
      for(int a = 0;a < l;++a)
         for(int b = a + 1;b < M;++b){

            if(b == l)
               ++b;

            if(b == M)
               break;

            for(int d = l + 1;d < M;++d){

               if(d == b)
                  ++d;

               if(d == M)
                  break;

               for(int S_ab = 0;S_ab < 2;++S_ab)
                  for(int S_cd = 0;S_cd < 2;++S_cd)
                     mat[S_ab][S_cd] = (*this)(l,0,S_ab,a,b,S_cd,l,d);

               for(int S_ab = 0;S_ab < 2;++S_ab)
                  for(int S_cd = 0;S_cd < 2;++S_cd){

                     double ward = mat[S_ab][S_cd];

                     for(int S_ld = 0;S_ld < 2;++S_ld)
                        ward += std::sqrt( (2.0*S_ab + 1.0) * (2.0*S_ld + 1.0) ) * _6j[S_cd][S_ld] * mat[S_ab][S_ld];

                     i = rTPM::gs2t(l,0,S_ab,a,b);
                     j = rTPM::gs2t(l,0,S_cd,l,d);

                     (*this)[l](0,i,j) = 0.5 * ward;
                     (*this)[l](0,j,i) = (*this)[l](0,i,j);

                  }

            }
         }

      //d = l: W^l_{ab;cl}
      for(int a = 0;a < M;++a)
         for(int b = a + 1;b < M;++b){

            if(b == l)
               ++b;

            if(b == M)
               break;

            for(int c = a + 1;c < l;++c){

               if(c == b)
                  ++c;

               if(c == l)
                  break;

               for(int S_ab = 0;S_ab < 2;++S_ab)
                  for(int S_cd = 0;S_cd < 2;++S_cd)
                     mat[S_ab][S_cd] = (*this)(l,0,S_ab,a,b,S_cd,c,l);

               for(int S_ab = 0;S_ab < 2;++S_ab)
                  for(int S_cd = 0;S_cd < 2;++S_cd){

                     double ward = mat[S_ab][S_cd];

                     for(int S_cl = 0;S_cl < 2;++S_cl)
                        ward += std::sqrt( (2.0*S_ab + 1.0) * (2.0*S_cl + 1.0) ) * (1 - 2*S_cd) * (1 - 2*S_cl) * _6j[S_cd][S_cl] * mat[S_ab][S_cl];

                     i = rTPM::gs2t(l,0,S_ab,a,b);
                     j = rTPM::gs2t(l,0,S_cd,c,l);

                     (*this)[l](0,i,j) = 0.5 * ward;
                     (*this)[l](0,j,i) = (*this)[l](0,i,j);

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

            i = rTPM::gs2t(l,1,1,a,b);

            (*this)[l](1,i,i) = 1.0/3.0 * ( (*this)[l](1,i,i) + (*this)(a,1,1,l,b,1,l,b) + (*this)(b,1,1,a,l,1,a,l) );

            //rest is symmetric
            i = rTPM::gs2t(a,1,1,l,b);

            (*this)[a](1,i,i) = (*this)(l,1,1,a,b,1,a,b);

            i = rTPM::gs2t(b,1,1,a,l);

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

               i = rTPM::gs2t(l,1,1,a,b);
               j = rTPM::gs2t(l,1,1,c,b);

               (*this)[l](1,i,j) = 0.5 * ( (*this)[l](1,i,j) + (*this)(b,1,1,a,l,1,c,l) );
               (*this)[l](1,j,i) = (*this)[l](1,i,j);

               i = rTPM::gs2t(b,1,1,a,l);
               j = rTPM::gs2t(b,1,1,c,l);

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

               i = rTPM::gs2t(l,1,1,a,b);
               j = rTPM::gs2t(l,1,1,b,c);

               (*this)[l](1,i,j) = 0.5 * ( (*this)[l](1,i,j) + (*this)(b,1,1,a,l,1,l,c) );
               (*this)[l](1,j,i) = (*this)[l](1,i,j);

               i = rTPM::gs2t(b,1,1,a,l);
               j = rTPM::gs2t(b,1,1,l,c);

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

               i = rTPM::gs2t(l,1,1,a,b);
               j = rTPM::gs2t(l,1,1,a,c);

               (*this)[l](1,i,j) = 0.5 * ( (*this)[l](1,i,j) + (*this)(a,1,1,l,b,1,l,c) );
               (*this)[l](1,j,i) = (*this)[l](1,i,j);

               i = rTPM::gs2t(a,1,1,l,b);
               j = rTPM::gs2t(a,1,1,l,c);

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
*/
   }

}

/**
 * project onto traceless space, watch out, special trace
 */
void dDPM::proj_fTr(){

   dDPM funit;
   funit.set_funit();

   double ward = this->ddot(funit)/(funit.ddot(funit));

   this->daxpy(-ward,funit);

}

/**
 * total projection, both on traceless space and good third index symmetry
 */
void dDPM::proj(){

   this->proj_W();
   this->proj_fTr();

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
         for(int b = 0;b < M;++b){
            
            for(int S_ab = 0;S_ab < 2;++S_ab)
               for(int S_cd = 0;S_cd < 2;++S_cd){

                     //1)
                     cout << a << "\t" << b << "\t(" << S_ab << ")\t(" << S_cd << ")\t" << (*this)(l,0,S_ab,a,l,S_cd,a,b) << "\t";

                     double ward = 0;

                     for(int S_al = 0;S_al < 2;++S_al)
                        ward += std::sqrt( (2.0*S_ab + 1.0) * (2.0*S_al + 1.0) ) * (1 - 2*S_al) * (1 - 2*S_ab) *_6j[S_ab][S_al] * (*this)(l,0,S_al,a,l,S_cd,a,b);

                     cout << ward << "\t";

                     ward = 0.0;

                     for(int S_ld = 0;S_ld < 2;++S_ld){

                        double hard = std::sqrt(2.0)*std::sqrt( (2.0*S_ab + 1.0) * (2.0*S_cd + 1.0) * (2.0*S_ld + 1.0) ) * _6j[S_ab][0]

                           *_6j[S_cd][S_ld] * (*this)(a,0,0,l,l,S_ld,l,b);

                        if(a == b)
                           hard /= std::sqrt(2.0);

                        if(l == b)
                           hard *= std::sqrt(2.0);

                        ward += hard;

                     }

                     cout << ward << endl;

               }

         }
/*
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
                        for(int S_cl = 0;S_cl < 2;++S_cl)
                           ward += std::sqrt( (2.0*S_al + 1.0) * (2.0*S_cl + 1.0) * (2.0*S_ab + 1.0) * (2.0*S_cd + 1.0) ) * (1 - 2*S_al) * (1 - 2*S_cl)

                              * (1 - 2*S_ab) * (1 - 2*S_cd) * _6j[S_ab][S_al] * _6j[S_cd][S_cl] * (*this)(b,0,S_al,a,l,S_cl,c,l);

                     cout << ward << endl;

                  }

               //2) b = c
               for(int S_ab = 0;S_ab < 2;++S_ab)
                  for(int S_cd = 0;S_cd < 2;++S_cd){

                     cout << "ab;bc\t" << a << "\t" << b << "\t" << c << "\t(" << S_ab << ")\t(" << S_cd << ")\t" << (*this)(l,0,S_ab,a,b,S_cd,b,c) << "\t";

                     double ward = 0.0;

                     for(int S_al = 0;S_al < 2;++S_al)
                        for(int S_ld = 0;S_ld < 2;++S_ld)
                           ward += std::sqrt( (2.0*S_al + 1.0) * (2.0*S_ld + 1.0) * (2.0*S_ab + 1.0) * (2.0*S_cd + 1.0) ) * (1 - 2*S_al) * 

                              (1 - 2*S_ab) * _6j[S_ab][S_al] * _6j[S_cd][S_ld] * (*this)(b,0,S_al,a,l,S_ld,l,c);

                     cout << ward << endl;

                  }

               //3) a = d
               for(int S_ab = 0;S_ab < 2;++S_ab)
                  for(int S_cd = 0;S_cd < 2;++S_cd){

                     cout << "ab;ca\t" << a << "\t" << b << "\t" << c << "\t(" << S_ab << ")\t(" << S_cd << ")\t" << (*this)(l,0,S_ab,a,b,S_cd,c,a) << "\t";

                     double ward = 0.0;

                     for(int S_lb = 0;S_lb < 2;++S_lb)
                        for(int S_cl = 0;S_cl < 2;++S_cl)
                           ward += std::sqrt( (2.0*S_lb + 1.0) * (2.0*S_cl + 1.0) * (2.0*S_ab + 1.0) * (2.0*S_cd + 1.0) ) * (1 - 2*S_cl)

                              * (1 - 2*S_cd) * _6j[S_ab][S_lb] * _6j[S_cd][S_cl] * (*this)(a,0,S_lb,l,b,S_cl,c,l);

                     cout << ward << endl;

                  }

               //4) a = c
               for(int S_ab = 0;S_ab < 2;++S_ab)
                  for(int S_cd = 0;S_cd < 2;++S_cd){

                     cout << "ab;ac\t" << a << "\t" << b << "\t" << c << "\t(" << S_ab << ")\t(" << S_cd << ")\t" << (*this)(l,0,S_ab,a,b,S_cd,a,c) << "\t";

                     double ward = 0.0;

                     for(int S_lb = 0;S_lb < 2;++S_lb)
                        for(int S_ld = 0;S_ld < 2;++S_ld)
                           ward += std::sqrt( (2.0*S_lb + 1.0) * (2.0*S_ld + 1.0) * (2.0*S_ab + 1.0) * (2.0*S_cd + 1.0) ) 

                              * _6j[S_ab][S_lb] * _6j[S_cd][S_ld] * (*this)(a,0,S_lb,l,b,S_ld,l,c);

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
*/
      cout << endl;

      /*
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
       */
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
                        for(int S_cl = 0;S_cl < 2;++S_cl)
                           ward += std::sqrt( (2.0*S_al + 1.0) * (2.0*S_cl + 1.0) * (2.0*S_ab + 1.0) * (2.0*S_cd + 1.0) ) * (1 - 2*S_al) * (1 - 2*S_cl)

                              * (1 - 2*S_ab) * (1 - 2*S_cd) * _6j[S_ab][S_al] * _6j[S_cd][S_cl] * (*this)(b,0,S_al,a,l,S_cl,c,l);

                     if( fabs(ward -  (*this)(l,0,S_ab,a,b,S_cd,c,b)) > 1.0e-12)
                        cout << "ab;cb\t" << a << "\t" << b << "\t" << c << "\t(" << S_ab << ")\t(" << S_cd << ")\t" << (*this)(l,0,S_ab,a,b,S_cd,c,b) << "\t" << ward << endl;

                  }

               //2) b = c
               for(int S_ab = 0;S_ab < 2;++S_ab)
                  for(int S_cd = 0;S_cd < 2;++S_cd){

                     double ward = 0.0;

                     for(int S_al = 0;S_al < 2;++S_al)
                        for(int S_ld = 0;S_ld < 2;++S_ld)
                           ward += std::sqrt( (2.0*S_al + 1.0) * (2.0*S_ld + 1.0) * (2.0*S_ab + 1.0) * (2.0*S_cd + 1.0) ) * (1 - 2*S_al) * 

                              (1 - 2*S_ab) * _6j[S_ab][S_al] * _6j[S_cd][S_ld] * (*this)(b,0,S_al,a,l,S_ld,l,c);

                     if(fabs(ward - (*this)(l,0,S_ab,a,b,S_cd,b,c)) > 1.0e-12)
                        cout << "ab;bc\t" << a << "\t" << b << "\t" << c << "\t(" << S_ab << ")\t(" << S_cd << ")\t" << (*this)(l,0,S_ab,a,b,S_cd,b,c) << "\t" << ward << endl;

                  }

               //3) a = d
               for(int S_ab = 0;S_ab < 2;++S_ab)
                  for(int S_cd = 0;S_cd < 2;++S_cd){

                     double ward = 0.0;

                     for(int S_lb = 0;S_lb < 2;++S_lb)
                        for(int S_cl = 0;S_cl < 2;++S_cl)
                           ward += std::sqrt( (2.0*S_lb + 1.0) * (2.0*S_cl + 1.0) * (2.0*S_ab + 1.0) * (2.0*S_cd + 1.0) ) * (1 - 2*S_cl)

                              * (1 - 2*S_cd) * _6j[S_ab][S_lb] * _6j[S_cd][S_cl] * (*this)(a,0,S_lb,l,b,S_cl,c,l);

                     if(fabs(ward - (*this)(l,0,S_ab,a,b,S_cd,c,a)) > 1.0e-12)
                        cout << "ab;ca\t" << a << "\t" << b << "\t" << c << "\t(" << S_ab << ")\t(" << S_cd << ")\t" << (*this)(l,0,S_ab,a,b,S_cd,c,a) << "\t" << ward << endl;

                  }

               //4) a = c
               for(int S_ab = 0;S_ab < 2;++S_ab)
                  for(int S_cd = 0;S_cd < 2;++S_cd){

                     double ward = 0.0;

                     for(int S_lb = 0;S_lb < 2;++S_lb)
                        for(int S_ld = 0;S_ld < 2;++S_ld)
                           ward += std::sqrt( (2.0*S_lb + 1.0) * (2.0*S_ld + 1.0) * (2.0*S_ab + 1.0) * (2.0*S_cd + 1.0) ) 

                              * _6j[S_ab][S_lb] * _6j[S_cd][S_ld] * (*this)(a,0,S_lb,l,b,S_ld,l,c);

                     if(fabs(ward - (*this)(l,0,S_ab,a,b,S_cd,a,c) ) > 1.0e-12)
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

                     if(fabs(ward -  (*this)(l,0,S_ab,a,l,S_cd,b,c)))
                        cout << "al;cd\t" << a << "\t" << b << "\t" << c << "\t(" << S_ab << ")\t(" << S_cd << ")\t" << (*this)(l,0,S_ab,a,l,S_cd,b,c) << "\t" << ward << endl;

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
      /*
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
       */
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

      //S = 1/2
      for(int i = 0;i < ddpm[l]->gdim(0);++i){

         S_ab = rTPM::gt2s(l,0,i,0);

         a = rTPM::gt2s(l,0,i,1);
         b = rTPM::gt2s(l,0,i,2);

         for(int j = 0;j < ddpm[l]->gdim(0);++j){

            S_cd = rTPM::gt2s(l,0,j,0);

            c = rTPM::gt2s(l,0,j,1);
            d = rTPM::gt2s(l,0,j,2);

            if(S_ab == 0 && S_cd == 0){

               (*this)[l](0,i,j) = tpm(0,a,b,c,d);

               if(a == b){

                  if(c == d){//a == b && c == d

                     if(a == c)
                        (*this)[l](0,i,j) += 0.5 * ( tpm(0,a,l,a,l) + 3.0 * tpm(1,a,l,a,l) );

                  }
                  else{//a == b && c != d

                     if(a == c)
                        (*this)[l](0,i,j) += 0.5*std::sqrt(0.5) * tpm(0,l,b,l,d);

                     if(a == d)
                        (*this)[l](0,i,j) += 0.5*std::sqrt(0.5) * tpm(0,l,b,l,c);

                     if(b == c)
                        (*this)[l](0,i,j) += 1.5*std::sqrt(0.5) * tpm(1,l,a,l,d);

                     if(b == d)
                        (*this)[l](0,i,j) += 1.5*std::sqrt(0.5) * tpm(1,l,a,l,c);

                  }

               }
               else{//a != b

                  if(c == d){//a != b && c == d

                     if(b == c)
                        (*this)[l](0,i,j) += 0.5*std::sqrt(0.5) * tpm(0,a,l,d,l);

                     if(a == c)
                        (*this)[l](0,i,j) += 0.5*std::sqrt(0.5) * tpm(0,b,l,d,l);

                     if(b == d)
                        (*this)[l](0,i,j) += 1.5*std::sqrt(0.5) * tpm(1,l,a,l,c);

                     if(a == d)
                        (*this)[l](0,i,j) += 1.5*std::sqrt(0.5) * tpm(1,l,b,l,c);

                  }

               }

            }
            else if(S_ab == 0 && S_cd == 1){

               (*this)[l](0,i,j) = 0.0;

               if(a == b){

                  if(b == c)
                     (*this)[l](0,i,j) += 0.5*std::sqrt(1.5) * tpm(1,l,a,l,d);

                  if(b == d)
                     (*this)[l](0,i,j) -= 0.5*std::sqrt(1.5) * tpm(1,l,a,l,c);

                  if(a == c)
                     (*this)[l](0,i,j) -= 0.5*std::sqrt(1.5) * tpm(0,l,b,l,d);

                  if(a == d)
                     (*this)[l](0,i,j) += 0.5*std::sqrt(1.5) * tpm(0,l,b,l,c);

               }

            }
            else if(S_ab == 1 && S_cd == 0){

               (*this)[l](0,i,j) = 0.0;

               if(c == d){

                  if(b == c)
                     (*this)[l](0,i,j) += 0.5*std::sqrt(1.5) * tpm(0,a,l,d,l);

                  if(a == c)
                     (*this)[l](0,i,j) -= 0.5*std::sqrt(1.5) * tpm(0,b,l,d,l);

                  if(b == d)
                     (*this)[l](0,i,j) -= 0.5*std::sqrt(1.5) * tpm(1,l,a,l,c);

                  if(a == d)
                     (*this)[l](0,i,j) += 0.5*std::sqrt(1.5) * tpm(1,l,b,l,c);

               }

            }
            else if(S_ab == 1 && S_cd == 1)
               (*this)[l](0,i,j) = tpm(1,a,b,c,d);

            (*this)[l](0,i,j) /= (N - 2.0);

         }
      }

      //S = 3/2
      for(int i = 0;i < ddpm[l]->gdim(1);++i){

         a = rTPM::gt2s(l,1,i,1);
         b = rTPM::gt2s(l,1,i,2);

         for(int j = 0;j < ddpm[l]->gdim(1);++j){

            c = rTPM::gt2s(l,1,j,1);
            d = rTPM::gt2s(l,1,j,2);

            (*this)[l](1,i,j) = tpm(1,a,b,c,d)/(N - 2.0);

         }
      }

   }

   //this->symmetrize();

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

   double ward = N*(N - 1.0)*(N - 2.0)/(2.0*M*(2.0*M - 1)*(2.0*M - 2.0));

   for(int l = 0;l < M;++l){

      for(int S = 0;S < 2;++S){

         for(int i = 0;i < ddpm[l]->gdim(S);++i){

            (*this)[l](S,i,i) = ward;

            for(int j = i + 1;j < ddpm[l]->gdim(S);++j)
               (*this)[l](S,i,j) = (*this)[l](S,j,i) = 0.0;

         }

      }

   }

}

/**
 * different trace that gets the right number "full" trace
 */
double dDPM::ftrace() const{

   double ward = 0.0;

   for(int l = 0;l < M;++l){

      //S = 1/2

      //1) S_ab = 0
      for(int a = 0;a < M;++a){

         if(a == l)
            ++a;

         if(a == M)
            break;

         for(int b = a;b < M;++b){

            if(b == l)
               ++b;

            if(b == M)
               break;

            ward += 2 * (*this)(l,0,0,a,b,0,a,b);

         }
      }

      //2) S_ab = 1;
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

            ward += 2 * (*this)(l,0,1,a,b,1,a,b);

         }
      }


      //S = 3/2
      for(int a = 0;a < M;++a)
         for(int b = a + 1;b < M;++b)
            ward += 4 * (*this)(l,1,1,a,b,1,a,b);

      //extra! extra!
      for(int a = 0;a < M;++a)
         ward += 4.0 * (*this)(a,0,0,l,l,0,l,l);

   }

   return ward;

}

/**
 * Construct the right hand side of the Newton equation for the determination of the search direction, 
 * the negative gradient of the potential:
 * @param t scaling factor of the potential
 * @param ham Hamiltonian of the current problem
 * @param P SUP matrix containing the inverse of the constraint matrices (carrier space matrices).
 */
void dDPM::constr_grad(double t,const dDPM &ham,const SUP &P){

   *this = P.gI1();

   this->dscal(t);

   *this -= ham;

   this->proj();

}

/** 
 * set (*this) equal to the full trace operator
 */
void dDPM::set_funit(){

   int a,b;

   for(int l = 0;l < M;++l)
      for(int S = 0;S < 2;++S)
         for(int i = 0;i < ddpm[l]->gdim(S);++i){

            a = rTPM::gt2s(l,S,i,1);
            b = rTPM::gt2s(l,S,i,2);

            if(a == b)
               (*this)[l](S,i,i) = 3.0;
            else
               (*this)[l](S,i,i) = 1.0;

            for(int j = i + 1;j < ddpm[l]->gdim(S);++j)
               (*this)[l](S,i,j) = (*this)[l](S,j,i) = 0.0;

         }


}

/**
 * solve the Newton equations for the determination of the search direction,
 * @param t scaling factor of the potential
 * @param P SUP matrix containing the inverse of the constraint matrices (carrier space matrices).
 * @param b right hand side (the gradient constructed int dDPM::constr_grad)
 * @return nr of iterations needed to converge to the desired accuracy
 */
int dDPM::solve(double t,const SUP &P,dDPM &b){

   int iter = 0;

   //delta = 0
   *this = 0;

   //residu:
   dDPM r(b);

   //norm van het residu
   double rr = r.ddot(r);

   //enkele variabelen
   double rr_old,ward;

   dDPM Hb;

   while(rr > 1.0e-7){ 

      Hb.H(t,b,P);

      ward = rr/b.ddot(Hb);

      //delta += ward*b
      this->daxpy(ward,b);

      //r -= ward*Hb
      r.daxpy(-ward,Hb);

      //nieuwe variabelen berekenen en oude overdragen
      rr_old = rr;
      rr = r.ddot(r);

      //nieuwe b nog:
      b.dscal(rr/rr_old);

      b += r;

      ++iter;

   }

   return iter;

}

/**
 * The hessian-map of the Newton system:
 * @param t potential scaling factor
 * @param b the dDPM on which the hamiltonian will work, the image will be put in (*this)
 * @param P the SUP matrix containing the constraints, (can be seen as the metric).
 */
void dDPM::H(double t,const dDPM &b,const SUP &P){

   this->L_map(P.gI1(),b);

   this->dscal(t);

   this->proj();

}

/**
 * perform a line search to determine what step size in along the Newton direction is ideal.
 * @param t potential scaling factor
 * @param P SUP matrix containing the inverse of the constraints (carrier space matrices)
 * @param ham Hamiltonian of the problem
 * @return the steplength
 */
double dDPM::line_search(double t,SUP &P,const dDPM &ham){

   double tolerance = 1.0e-5*t;

   if(tolerance < 1.0e-12)
      tolerance = 1.0e-12;

   //neem de wortel uit P
   P.sqrt(1);

   //maak eerst een SUP van delta
   SUP S_delta;

   S_delta.fill(*this);

   //hulpje om dingskes in te steken:
   SUP hulp;
   hulp.L_map(P,S_delta);

   EIG eigen(hulp);

   double a = 0;

   double b = -1.0/eigen.min();

   double c(0);

   double ham_delta = ham.ddot(*this);

   while(b - a > tolerance){

      c = (b + a)/2.0;

      if( (ham_delta - t*eigen.lsfunc(c)) < 0.0)
         a = c;
      else
         b = c;

   }

   return c;

}

/**
 * perform a line search what step size in along a certain direction (*this) minimizes the potential, this one is used for extrapolation.
 * @param t potential scaling factor
 * @param W dDPM containing the current approximation of the W object
 * @param ham Hamiltonian of the problem
 * @return the steplength
 */
double dDPM::line_search(double t,const dDPM &W,const dDPM &ham){

   SUP P;

   P.fill(W);

   P.invert();

   return this->line_search(t,P,ham);

}

/**
 * a function that returns the inproduct of a dDPM object with a TPM object. Test for ddot function.
 * Needed because there are problems when one of the sp indices equals l
 * @param input TPM object
 * @return the inproduct
 */
double dDPM::inprod(const TPM &tpm){

   double ward = 0.0;

   for(int l = 0;l < M;++l){

      for(int a = 0;a < M;++a)
         for(int b = a;b < M;++b)
            for(int c = 0;c < M;++c)
               for(int d = c;d < M;++d)
                  ward += 2.0 * (*this)(l,0,0,a,b,0,c,d) * tpm(0,a,b,c,d);

      for(int a = 0;a < M;++a)
         for(int b = a + 1;b < M;++b)
            for(int c = 0;c < M;++c)
               for(int d = c + 1;d < M;++d)
                  ward += 2.0 * (*this)(l,0,1,a,b,1,c,d) * tpm(1,a,b,c,d);

      for(int a = 0;a < M;++a)
         for(int b = a + 1;b < M;++b)
            for(int c = 0;c < M;++c)
               for(int d = c + 1;d < M;++d)
                  ward += 4.0 * (*this)(l,1,1,a,b,1,c,d) * tpm(1,a,b,c,d);

   }

   return ward/(N - 2.0);

}
