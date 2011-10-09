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

int ssdTPM::M;
int ssdTPM::N;

/**
 * initialize the static variables
 * @param M_in dimension of sp space
 * @param N_in nr of particles
 */
void ssdTPM::init(int M_in,int N_in){

   M = M_in;
   N = N_in;

}

/**
 * standard constructor: constructs M SPM object with parameter l = 0 -> M-1
 */
ssdTPM::ssdTPM() {

   ssdtpm = new SPM * [M];

   for(int l = 0;l < M;++l)
      ssdtpm[l] = new SPM();

}

/**
 * copy constructor: constructs an array and copies the content of ssdtpm_c in it.
 * @param ssdtpm_c object that will be copied into this.
 */
ssdTPM::ssdTPM(const ssdTPM &ssdtpm_c) { 

   ssdtpm = new SPM * [M];

   for(int l = 0;l < M;++l)
      ssdtpm[l] = new SPM(ssdtpm_c[l]);
   
}

/**
 * destructor
 */
ssdTPM::~ssdTPM(){

   for(int l = 0;l < M;++l)
      delete ssdtpm[l];

   delete [] ssdtpm;
   
}

/**
 * @return nr of particles
 */
int ssdTPM::gN() const {

   return N;

}

/**
 * @return dimension of sp space
 */
int ssdTPM::gM() const {

   return M;

}

/**
 * acces to the individual SPM objects
 * @param l the specific SPM object you want
 * @return the SPM object with parameter l
 */
SPM &ssdTPM::operator[](int l){

   return *ssdtpm[l];

}

/**
 * acces to the individual SPM objects: the const version
 * @param l the specific SPM object you want
 * @return the SPM object with parameter l
 */
const SPM &ssdTPM::operator[](int l) const{

   return *ssdtpm[l];

}

ostream &operator<<(ostream &output,const ssdTPM &ssdtpm_p){

   for(int l = 0;l < ssdtpm_p.gM();++l){

      output << std::endl;
      output << l << std::endl;
      output << std::endl;

      output << ssdtpm_p[l] << endl;

   }

   return output;

}

/**
 * @return the trace of the ssdTPM object
 */
double ssdTPM::trace() const{

   double ward = 0.0;

   for(int l = 0;l < M;++l)
      ward += ssdtpm[l]->trace();

   return ward;

}

/**
 * Scale the ssdTPM with parameter alpha
 * @param alpha scalefactor
 */
void ssdTPM::dscal(double alpha){

   for(int l = 0;l < M;++l)
      ssdtpm[l]->dscal(alpha);

}

/**
 * overload the equality operator
 * @param ssdtpm_c The ssdTPM you want to be copied into this
 */
ssdTPM &ssdTPM::operator=(const ssdTPM &ssdtpm_c){

   for(int l = 0;l < M;++l)
      *ssdtpm[l] = ssdtpm_c[l];

   return *this;

}

/**
 * Make all the number in your ssdTPM equal to the number a, e.g. usefull for initialization (ssdTPM ssdtpm = 0)
 * @param a the number
 */
ssdTPM &ssdTPM::operator=(double a){

   for(int l = 0;l < M;++l)
      *ssdtpm[l] = a;

   return *this;

}

/**
 * overload the += operator for ssdTPM's
 * @param ssdtpm_pl The ssdTPM you want to add to this
 */
ssdTPM &ssdTPM::operator+=(const ssdTPM &ssdtpm_pl){

   for(int l = 0;l < M;++l)
      *ssdtpm[l] += ssdtpm_pl[l];

   return *this;

}

/**
 * overload the -= operator for ssdTPM's
 * @param ssdtpm_m The ssdTPM you want to deduct from this
 */
ssdTPM &ssdTPM::operator-=(const ssdTPM &ssdtpm_m){

   for(int l = 0;l < M;++l)
      *ssdtpm[l] -= ssdtpm_m[l];

   return *this;

}

/**
 * add the ssdtpm ssdtpm_pl times the constant alpha to *this
 * @param alpha the constant to multiply the ssdtpm_pl with
 * @param ssdtpm_pl the ssdTPM to be multiplied by alpha and added to this
 */
ssdTPM &ssdTPM::daxpy(double alpha,const ssdTPM &ssdtpm_pl){

   for(int l = 0;l < M;++l)
      ssdtpm[l]->daxpy(alpha,ssdtpm_pl[l]);

   return *this;

}

/**
 * Multiply symmetric ssdtpm object left en right with symmetric ssdtpm map to 
 * form another symmetric ssdtpm and put it in (*this): this = map*object*map
 * @param map ssdtpm that will be multiplied to the left en to the right of ssdtpm object
 * @param object central ssdtpm
 */
void ssdTPM::L_map(const ssdTPM &map,const ssdTPM &object){

   for(int l = 0;l < M;++l)
      ssdtpm[l]->L_map(map[l],object[l]);
   
}

/**
 * ssdTPM product of two general matrices A en B, put result in this
 * @param A left ssdtpm
 * @param B right ssdtpm
 */
ssdTPM &ssdTPM::mprod(const ssdTPM &A,const ssdTPM &B){

   for(int l = 0;l < M;++l)
      ssdtpm[l]->mprod(A[l],B[l]);

   return *this;

}

/**
 * Take the square root out of the positive semidefinite ssdtpm, destroys original ssdtpm, square root will be put in (*this)
 * @param option = 1, positive square root, = -1, negative square root.
 */
void ssdTPM::sqrt(int option){

   for(int l = 0;l < M;++l)
      ssdtpm[l]->sqrt(option);

}

/**
 * @return inproduct of (*this) ssdtpm with ssdtpm_i, defined as Tr (A B)
 * @param ssdtpm_i input ssdtpm
 */
double ssdTPM::ddot(const ssdTPM &ssdtpm_i) const{

   double ward = 0.0;

   for(int l = 0;l < M;++l)
      ward += ssdtpm[l]->ddot(ssdtpm_i[l]);

   return ward;

}

/**
 * Invert positive semidefinite symmetric ssdtpm is stored in (*this), original ssdtpm (*this) is destroyed
 */
void ssdTPM::invert(){

   for(int l = 0;l < M;++l)
      ssdtpm[l]->invert();

}

/**
 * copy upper in lower part of ssdTPM object
 */
void ssdTPM::symmetrize(){

   for(int l = 0;l < M;++l)
      ssdtpm[l]->symmetrize();

}

/**
 * Fill the matrix with random numbers.
 */
void ssdTPM::fill_Random(){

   for(int l = 0;l < M;++l)
      ssdtpm[l]->fill_Random();

}

/**
 * map a dDPM matrix on a ssdTPM by using a special "bar" function, the different blocks of a dDPM matrix.
 * special because the intermediate spin is not required to be diagonal, instead the usual [][]{...} is used.
 * @param scale the ssdTPM with this number
 * @param ddpm input dDPM object
 */
void ssdTPM::bar(double scale,const dDPM &ddpm){

   double ward,hard;

   for(int l = 0;l < M;++l){

      for(int a = 0;a < M;++a)
         for(int c = a;c < M;++c){

            (*this)[l](a,c) = 0.0;

            //first S = 1/2
            for(int S_ab = 0;S_ab < 2;++S_ab)
               for(int S_cd = 0;S_cd < 2;++S_cd){

                  ward = 0.0;

                  for(int b = 0;b < M;++b){

                     hard = ddpm(l,0,S_ab,a,b,S_cd,c,b);

                     if(a == b)
                        hard *= std::sqrt(2.0);

                     if(c == b)
                        hard *= std::sqrt(2.0);

                     ward += hard;

                  }

                  (*this)[l](a,c) += std::sqrt( (2.0*S_ab + 1.0) * (2.0*S_cd + 1.0) ) * Tools::g6j(1,1,2*S_ab,1,1,2*S_cd) * ward;

               }

            //then S = 3/2
            for(int b = 0;b < M;++b)
               (*this)[l](a,c) -= 2.0 * ddpm(l,1,1,a,b,1,c,b);


            //finally scale
            (*this)[l](a,c) *= scale;

         }

   }

   this->symmetrize();

}

/**
 * map a dPPHM matrix on a ssdTPM by using a special "bar" function, the different blocks of a dPPHM matrix.
 * special because the intermediate spin is not required to be diagonal, instead the usual [][] is used.
 * @param scale the ssdTPM with this number
 * @param dpphm input dPPHM object
 */
void ssdTPM::bar(double scale,const dPPHM &dpphm){

   double ward,hard;

   for(int l = 0;l < M;++l){

      for(int a = 0;a < M;++a)
         for(int c = a;c < M;++c){

            (*this)[l](a,c) = 0.0;

            for(int S_ab = 0;S_ab < 2;++S_ab)
               for(int S_cd = 0;S_cd < 2;++S_cd){

                  ward = 0.0;

                  for(int b = 0;b < M;++b){

                     hard = dpphm(l,0,S_ab,a,b,S_cd,c,b);

                     if(a == b)
                        hard *= std::sqrt(2.0);

                     if(c == b)
                        hard *= std::sqrt(2.0);

                     ward += hard;

                  }

                  (*this)[l](a,c) += std::sqrt( (2.0*S_ab + 1.0) * (2.0*S_cd + 1.0) ) * ward;

               }

            //finally scale
            (*this)[l](a,c) *= scale;

         }

   }

   this->symmetrize();

}

/**
 * map a dPHHM matrix on a ssdTPM by using a spinsummed "skew-bar" function on the different blocks of a dPHHM matrix.
 * spinsummed because the intermediate spin is not required to be diagonal, instead the usual [][] is used.
 * BEWARE: matrix not symmetric in general!
 * @param scale the ssdTPM with this number
 * @param dphhm input dPHHM object
 */
void ssdTPM::skew_bar(double scale,const dPHHM &dphhm){

   double ward;

   for(int l = 0;l < M;++l){

      for(int a = 0;a < M;++a)
         for(int c = 0;c < M;++c){

            (*this)[l](a,c) = 0.0;

            for(int S_bl = 0;S_bl < 2;++S_bl)
               for(int S_dl = 0;S_dl < 2;++S_dl){

                  ward = 0.0;

                  for(int b = 0;b < M;++b)
                     ward += dphhm(l,0,S_bl,b,b,S_dl,a,c);

                  (*this)[l](a,c) += (1 - 2*S_bl) * std::sqrt( (2.0*S_bl + 1.0) * (2.0*S_dl + 1.0) ) * ward;

               }

            //finally scale
            (*this)[l](a,c) *= scale;

         }

   }

}
