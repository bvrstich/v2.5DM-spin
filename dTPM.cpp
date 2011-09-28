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

int dTPM::M;
int dTPM::N;

/**
 * initialize the static variables
 * @param M_in dimension of sp space
 * @param N_in nr of particles
 */
void dTPM::init(int M_in,int N_in){

   M = M_in;
   N = N_in;

}

/**
 * standard constructor: constructs M SPM object with parameter l = 0 -> M-1
 */
dTPM::dTPM() {

   dtpm = new SPM * [M];

   for(int l = 0;l < M;++l)
      dtpm[l] = new SPM();

}

/**
 * copy constructor: constructs an array and copies the content of dtpm_c in it.
 * @param dtpm_c object that will be copied into this.
 */
dTPM::dTPM(const dTPM &dtpm_c) { 

   dtpm = new SPM * [M];

   for(int l = 0;l < M;++l)
      dtpm[l] = new SPM(dtpm_c[l]);
   
}

/**
 * destructor
 */
dTPM::~dTPM(){

   for(int l = 0;l < M;++l)
      delete dtpm[l];

   delete [] dtpm;
   
}

/**
 * @return nr of particles
 */
int dTPM::gN() const {

   return N;

}

/**
 * @return dimension of sp space
 */
int dTPM::gM() const {

   return M;

}

/**
 * acces to the individual SPM objects
 * @param l the specific SPM object you want
 * @return the SPM object with parameter l
 */
SPM &dTPM::operator[](int l){

   return *dtpm[l];

}

/**
 * acces to the individual SPM objects: the const version
 * @param l the specific SPM object you want
 * @return the SPM object with parameter l
 */
const SPM &dTPM::operator[](int l) const{

   return *dtpm[l];

}

ostream &operator<<(ostream &output,const dTPM &dtpm_p){

   for(int l = 0;l < dtpm_p.gM();++l){

      output << std::endl;
      output << l << std::endl;
      output << std::endl;

      output << dtpm_p[l] << endl;

   }

   return output;

}

/**
 * @return the trace of the dTPM object
 */
double dTPM::trace() const{

   double ward = 0.0;

   for(int l = 0;l < M;++l)
      ward += dtpm[l]->trace();

   return ward;

}

/**
 * Scale the dTPM with parameter alpha
 * @param alpha scalefactor
 */
void dTPM::dscal(double alpha){

   for(int l = 0;l < M;++l)
      dtpm[l]->dscal(alpha);

}

/**
 * overload the equality operator
 * @param dtpm_c The dTPM you want to be copied into this
 */
dTPM &dTPM::operator=(const dTPM &dtpm_c){

   for(int l = 0;l < M;++l)
      *dtpm[l] = dtpm_c[l];

   return *this;

}

/**
 * Make all the number in your dTPM equal to the number a, e.g. usefull for initialization (dTPM dtpm = 0)
 * @param a the number
 */
dTPM &dTPM::operator=(double a){

   for(int l = 0;l < M;++l)
      *dtpm[l] = a;

   return *this;

}

/**
 * overload the += operator for dTPM's
 * @param dtpm_pl The dTPM you want to add to this
 */
dTPM &dTPM::operator+=(const dTPM &dtpm_pl){

   for(int l = 0;l < M;++l)
      *dtpm[l] += dtpm_pl[l];

   return *this;

}

/**
 * overload the -= operator for dTPM's
 * @param dtpm_m The dTPM you want to deduct from this
 */
dTPM &dTPM::operator-=(const dTPM &dtpm_m){

   for(int l = 0;l < M;++l)
      *dtpm[l] -= dtpm_m[l];

   return *this;

}

/**
 * add the dtpm dtpm_pl times the constant alpha to *this
 * @param alpha the constant to multiply the dtpm_pl with
 * @param dtpm_pl the dTPM to be multiplied by alpha and added to this
 */
dTPM &dTPM::daxpy(double alpha,const dTPM &dtpm_pl){

   for(int l = 0;l < M;++l)
      dtpm[l]->daxpy(alpha,dtpm_pl[l]);

   return *this;

}

/**
 * Multiply symmetric dtpm object left en right with symmetric dtpm map to 
 * form another symmetric dtpm and put it in (*this): this = map*object*map
 * @param map dtpm that will be multiplied to the left en to the right of dtpm object
 * @param object central dtpm
 */
void dTPM::L_map(const dTPM &map,const dTPM &object){

   for(int l = 0;l < M;++l)
      dtpm[l]->L_map(map[l],object[l]);
   
}

/**
 * dTPM product of two general matrices A en B, put result in this
 * @param A left dtpm
 * @param B right dtpm
 */
dTPM &dTPM::mprod(const dTPM &A,const dTPM &B){

   for(int l = 0;l < M;++l)
      dtpm[l]->mprod(A[l],B[l]);

   return *this;

}

/**
 * Take the square root out of the positive semidefinite dtpm, destroys original dtpm, square root will be put in (*this)
 * @param option = 1, positive square root, = -1, negative square root.
 */
void dTPM::sqrt(int option){

   for(int l = 0;l < M;++l)
      dtpm[l]->sqrt(option);

}

/**
 * @return inproduct of (*this) dtpm with dtpm_i, defined as Tr (A B)
 * @param dtpm_i input dtpm
 */
double dTPM::ddot(const dTPM &dtpm_i) const{

   double ward = 0.0;

   for(int l = 0;l < M;++l)
      ward += dtpm[l]->ddot(dtpm_i[l]);

   return ward;

}

/**
 * Invert positive semidefinite symmetric dtpm is stored in (*this), original dtpm (*this) is destroyed
 */
void dTPM::invert(){

   for(int l = 0;l < M;++l)
      dtpm[l]->invert();

}

/**
 * copy upper in lower part of dTPM object
 */
void dTPM::symmetrize(){

   for(int l = 0;l < M;++l)
      dtpm[l]->symmetrize();

}

/**
 * Fill the matrix with random numbers.
 */
void dTPM::fill_Random(){

   for(int l = 0;l < M;++l)
      dtpm[l]->fill_Random();

}

/**
 * map a dDPM matrix on a dTPM by free-"barring" the different blocks of a dDPM matrix, it is called freebar
 * because the intermediate spin is not required to be diagonal, instead the usual [][]{...} is used.
 * @param scale the dTPM with this number
 * @param ddpm input dDPM object
 */
void dTPM::freebar(double scale,const dDPM &ddpm){

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

                  (*this)[l](a,c) += std::sqrt( (2.0*S_ab + 1.0) * (2.0*S_cd + 1.0) ) * Tools::g6j(S_ab,S_cd) * ward;

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
