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

int dSPM::M;
int dSPM::N;

/**
 * initialize the static variables
 * @param M_in dimension of sp space
 * @param N_in nr of particles
 */
void dSPM::init(int M_in,int N_in){

   M = M_in;
   N = N_in;

}

/**
 * standard constructor: constructs M rSPM object with parameter l = 0 -> M-1
 */
dSPM::dSPM() {

   dspm = new double [M];

}

/**
 * copy constructor: constructs an array and copies the content of dspm_c in it.
 * @param dspm_c object that will be copied into this.
 */
dSPM::dSPM(const dSPM &dspm_c) { 

   dspm = new double [M];

#pragma omp parallel for
   for(int l = 0;l < M;++l)
      dspm[l] = dspm_c[l];
   
}

/**
 * destructor
 */
dSPM::~dSPM(){

   delete [] dspm;
   
}

/**
 * @return nr of particles
 */
int dSPM::gN() const {

   return N;

}

/**
 * @return dimension of spatial sp space
 */
int dSPM::gM() const {

   return M;

}

/**
 * acces to the individual numbers
 * @param l the specific number object you want
 * @return the number
 */
double &dSPM::operator[](int l){

   return dspm[l];

}

/**
 * acces to the individual numbers
 * @param l the specific number object you want
 * @return the number
 */
double dSPM::operator[](int l) const {

   return dspm[l];

}

ostream &operator<<(ostream &output,const dSPM &dspm_p){

   for(int l = 0;l < dspm_p.gM();++l)
      output << l << "\t" << dspm_p[l] << endl;

   return output;

}

/**
 * construct the dSPM object from a dDPM by tracing the different blocks
 * @param scale all the traces with this value
 * @param ddpm the input dDPM
 */
void dSPM::trace(double scale,const dDPM &ddpm){

#pragma omp parallel for
   for(int l = 0;l < M;++l)
      dspm[l] = scale * ddpm[l].trace();

}

/**
 * construct the dSPM object from a dPPHM by tracing the different blocks
 * @param scale all the traces with this value
 * @param dpphm the input dPPHM
 */
void dSPM::trace(double scale,const dPPHM &dpphm){

#pragma omp parallel for
   for(int l = 0;l < M;++l)
      dspm[l] = scale * dpphm[l].trace();

}

/**
 * construct the dSPM object from a dPHHM by skew-tracing the different blocks
 * @param scale all the traces with this value
 * @param dphhm the input dPHHM
 */
void dSPM::skew_trace(double scale,const dPHHM &dphhm){

#pragma omp parallel for
   for(int l = 0;l < M;++l){

      dspm[l] = 0.0;

      for(int S_bl = 0;S_bl < 2;++S_bl)
         for(int S_dl = 0;S_dl < 2;++S_dl){

            double ward = 0.0;

            for(int a = 0;a < M;++a)
               for(int c = 0;c < M;++c)
                  ward += dphhm(l,0,S_bl,a,a,S_dl,c,c);

            dspm[l] += std::sqrt( (2.0*S_bl + 1.0) * (2.0*S_dl + 1.0) ) * (1 - 2*S_bl) * (1 - 2*S_dl) * ward;

         }

      dspm[l] *= scale;

   }

}
