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
 * standard constructor: constructs M xSPM object with parameter l = 0 -> M-1
 */
dTPM::dTPM() {

   dtpm = new xSPM * [M];

   for(int l = 0;l < M;++l)
      dtpm[l] = new xSPM();

}

/**
 * copy constructor: constructs an array and copies the content of dtpm_c in it.
 * @param dtpm_c object that will be copied into this.
 */
dTPM::dTPM(const dTPM &dtpm_c) { 

   dtpm = new xSPM * [M];

   for(int l = 0;l < M;++l)
      dtpm[l] = new xSPM(dtpm_c[l]);
   
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
 * acces to the individual xSPM objects
 * @param l the specific xSPM object you want
 * @return the xSPM object with parameter l
 */
xSPM &dTPM::operator[](int l){

   return *dtpm[l];

}

/**
 * acces to the individual xSPM objects: the const version
 * @param l the specific xSPM object you want
 * @return the xSPM object with parameter l
 */
const xSPM &dTPM::operator[](int l) const{

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
 * special "bar" function that maps a dDPM on a dTPM object, see notes for info.
 * @param scale the factor you scale the dTPM with
 * @param ddpm input dDPM object
 */
void dTPM::bar(double scale,const dDPM &ddpm){

   double ward,hard;

   for(int l = 0;l < M;++l){

      for(int b = 0;b < M;++b)
         for(int d = b;d < M;++d){

            //first Z = 0
            (*this)[l](0,b,d) = 0.0;

            //only S = 1/2 contribution
            for(int S_ab = 0;S_ab < 2;++S_ab)
               for(int S_cd = 0;S_cd < 2;++S_cd){

                  ward = 0.0;

                  for(int a = 0;a < M;++a){

                     hard = ddpm(l,0,S_ab,a,b,S_cd,a,d);

                     if(a == b)
                        hard *= std::sqrt(2.0);

                     if(a == d)
                        hard *= std::sqrt(2.0);

                     ward += hard;

                  }

                  (*this)[l](0,b,d) += 2.0 * std::sqrt( (2.0*S_ab + 1.0) * (2.0*S_cd + 1.0) ) * Tools::g6j(0,0,S_ab,0) * Tools::g6j(0,0,S_cd,0) * ward;

               }

            (*this)[l](0,b,d) *= scale;

            //then Z = 1
            (*this)[l](1,b,d) = 0.0;

            //first S = 1/2 contribution
            for(int S_ab = 0;S_ab < 2;++S_ab)
               for(int S_cd = 0;S_cd < 2;++S_cd){

                  ward = 0.0;

                  for(int a = 0;a < M;++a){

                     hard = ddpm(l,0,S_ab,a,b,S_cd,a,d);

                     if(a == b)
                        hard *= std::sqrt(2.0);

                     if(a == d)
                        hard *= std::sqrt(2.0);

                     ward += hard;

                  }

                  (*this)[l](1,b,d) += 2.0 * std::sqrt( (2.0*S_ab + 1.0) * (2.0*S_cd + 1.0) ) * Tools::g6j(0,0,S_ab,1) * Tools::g6j(0,0,S_cd,1) * ward;

               }

            //finally S = 3/2 contribution
            for(int a = 0;a < M;++a)
               (*this)[l](1,b,d) += 4.0/3.0 * ddpm(l,1,1,a,b,1,a,d);

            (*this)[l](1,b,d) *= scale;

         }

   }

   this->symmetrize();

}

/**
 * special "bar" function that maps a dPPHM on a dTPM object, see notes for info.
 * @param scale the factor you scale the dTPM with
 * @param dpphm input dPPHM object
 */
void dTPM::bar(double scale,const dPPHM &dpphm){

   double ward,hard;

   for(int l = 0;l < M;++l){

      for(int Z = 0;Z < 2;++Z){

         for(int b = 0;b < M;++b)
            for(int d = b;d < M;++d){

               (*this)[l](Z,b,d) = 0.0;

               for(int S = 0;S < 2;++S){

                  for(int S_ab = 0;S_ab < 2;++S_ab)
                     for(int S_cd = 0;S_cd < 2;++S_cd){

                        ward = 0.0;

                        for(int a = 0;a < M;++a){

                           hard = dpphm(l,S,S_ab,a,b,S_cd,a,d);

                           if(a == b)
                              hard *= std::sqrt(2.0);

                           if(a == d)
                              hard *= std::sqrt(2.0);

                           ward += hard;

                        }

                        (*this)[l](Z,b,d) += (2*(S + 0.5) + 1.0) * std::sqrt( (2.0*S_ab + 1.0) * (2.0*S_cd + 1.0) ) * 
                        
                           Tools::g9j(S,Z,S_ab,S_cd) * ward;

                     }

               }

               (*this)[l](Z,b,d) *= scale;

            }

      }
   }

   this->symmetrize();

}

/**
 * special "bar" function that maps a dPHHM on a dTPM object, see notes for info.
 * @param scale the factor you scale the dTPM with
 * @param dphhm input dPHHM object
 */
void dTPM::bar(double scale,const dPHHM &dphhm){

   for(int l = 0;l < M;++l){

      for(int Z = 0;Z < 2;++Z){

         for(int b = 0;b < M;++b)
            for(int d = b;d < M;++d){

               (*this)[l](Z,b,d) = 0.0;

               for(int S = 0;S < 2;++S)
                  for(int a = 0;a < M;++a)
                     (*this)[l](Z,b,d) += (2*(S + 0.5) + 1.0)/(2*Z + 1.0) * dphhm(l,S,Z,a,b,Z,a,d);
               

               (*this)[l](Z,b,d) *= scale;

            }

      }
   }

   this->symmetrize();

}


/**
 * special "skew-bar" function that maps a dPHHM on a dTPM object, see notes for info.
 * BEWARE: not symmetrical
 * @param scale the factor you scale the dTPM with
 * @param dphhm input dPHHM object
 */
void dTPM::skew_bar(double scale,const dPHHM &dphhm){

   double ward;

   for(int l = 0;l < M;++l){

      for(int S_dl = 0;S_dl < 2;++S_dl){

         for(int b = 0;b < M;++b)
            for(int d = 0;d < M;++d){

               (*this)[l](S_dl,b,d) = 0.0;

               for(int S_bl = 0;S_bl < 2;++S_bl){

                  ward = 0.0;

                  for(int a = 0;a < M;++a)
                     ward += dphhm(l,0,S_bl,a,a,S_dl,b,d);

                  (*this)[l](S_dl,b,d) += std::sqrt( (2*S_bl + 1.0) / (2*S_dl + 1.0) ) * (1 - 2*S_dl) * (1 - 2*S_bl) * ward;

               }

               (*this)[l](S_dl,b,d) *= scale;

            }

      }
   }

}

/**
 * special spinsummed "bar" function that maps a dPHHM on a dTPM object, see notes for info.
 * spinsummed because the intermediate spin is not required to be diagonal, and to have a different name then the other bar...
 * @param scale the factor you scale the dTPM with
 * @param dphhm input dPPHM object
 */
void dTPM::ssbar(double scale,const dPHHM &dphhm){

   double ward;

   for(int l = 0;l < M;++l){

      for(int Z = 0;Z < 2;++Z){

         for(int a = 0;a < M;++a)
            for(int c = a;c < M;++c){

               (*this)[l](Z,a,c) = 0.0;

               for(int S = 0;S < 2;++S)
                  for(int S_bl = 0;S_bl < 2;++S_bl)
                     for(int S_dl = 0;S_dl < 2;++S_dl){

                        ward = 0.0;

                        for(int b = 0;b < M;++b)
                           ward += dphhm(l,S,S_bl,a,b,S_dl,c,b);

                        (*this)[l](Z,a,c) += (2.0*(S + 0.5) + 1.0) * std::sqrt( (2.0*S_bl + 1.0) * (2.0*S_dl + 1.0) ) * (1 - 2*S_bl) * (1 - 2*S_dl)

                           * Tools::g9j(S,Z,S_dl,S_bl) * ward;

                     }


               (*this)[l](Z,a,c) *= scale;

            }

      }
   }

   this->symmetrize();

}
