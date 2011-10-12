#include <iostream>
#include <fstream>
#include <cmath>
#include <gsl/gsl_sf_coupling.h>

#include "include.h"

using std::cout;
using std::endl;

dDPM *Tools::unit;

double *Tools::x6j;
double *Tools::x9j;

double *Tools::C;

int Tools::M;
int Tools::N;


/**
 * initialize and allocate the static variables
 * @param M_in nr of sites
 * @param N_in nr of particles
 */
void Tools::init(int M_in,int N_in){

   M = M_in;
   N = N_in;

   x6j = new double [16];

   x9j = new double [16];

   C = new double [16];

   //fill the x6j list
   for(int S_ = 0;S_ < 2;++S_)
      for(int S = 0;S < 2;++S)
         for(int S_ab = 0;S_ab < 2;++S_ab)
            for(int S_cd = 0;S_cd < 2;++S_cd)
               x6j[8*S_ + 4*S + 2*S_ab + S_cd] = gsl_sf_coupling_6j(2*S_ + 1,1,2*S_ab,2*S + 1,1,2*S_cd);

   //fill the x9j list
   for(int S = 0;S < 2;++S)
      for(int Z = 0;Z < 2;++Z)
         for(int S_ab = 0;S_ab < 2;++S_ab)
            for(int S_cd = 0;S_cd < 2;++S_cd)
               x9j[8*S + 4*Z + 2*S_ab + S_cd] = gsl_sf_coupling_9j(1,1,2*Z,2*S + 1,2*S_ab,1,2*S_cd,1,1);

   //fill the C list
   for(int S = 0;S < 2;++S)
      for(int S_ = 0;S_ < 2;++S_)
         for(int S_ab = 0;S_ab < 2;++S_ab)
            for(int S_cd = 0;S_cd < 2;++S_cd){

               C[8*S + 4*S_ + 2*S_ab + S_cd] = 0.0;

               for(int j = 0;j < 3;++j){

                  C[8*S + 4*S_ + 2*S_ab + S_cd] += (2*j + 1.0) * (1 - 2*(j%2)) * gsl_sf_coupling_6j(2*S_ab,2*S_cd,2*j,2*S_ + 1,2*S + 1,1) 
                  
                     * gsl_sf_coupling_6j(2*S_ab,2*S_cd,2*j,2*S + 1,2*S_ + 1,1);

               }

            }

   unit = new dDPM();

   *unit = 0.0;

   for(int l = 0;l < M;++l)
      for(int S = 0;S < 2;++S)
         for(int i = 0;i < (*unit)[l].gdim(S);++i)
            (*unit)[l](S,i,i) = 1.0;

   unit->proj_W();

}

/**
 * deallocate the static lists and objects
 */
void Tools::clear(){

   delete unit;

   delete [] x6j;

   delete [] x9j;

   delete [] C;

}

/**
 * @return the projected unit matrix
 */
dDPM &Tools::gunit(){

   return *unit;

}

/**
 * function that returns stuff needed by the I2 map
 * @param S first dp spin index (0 or 1 means 1/2 or 3/2)
 * @param S_ second dp spin index (0 or 1 means 1/2 or 3/2)
 * @param S_ab intermediate spinindex (0 or 1 means 0 or 1)
 * @param S_cd intermediate spinindex (0 or 1 means 0 or 1)
 * @return the entry with indices S,S_,S_ab,S_cd
 */
double Tools::gC(int S,int S_,int S_ab,int S_cd){

   return C[8*S + 4*S_ + 2*S_ab + S_cd];

}

/**
 * store only the six j symbols needed
 * @param S_ 0 or 1 means 1/2 or 3/2
 * @param S 0 or 1 means 1/2 or 3/2
 * @param S_ab 0 or 1 means 0 or 1
 * @param S_cd 0 or 1 means 0 or 1
 * @return the 6j symbol
 */
double Tools::g6j(int S_,int S,int S_ab,int S_cd){

   return x6j[8*S_ + 4*S + 2*S_ab + S_cd];

}

/**
 * store only the nine j symbols needed
 * @param S 0 or 1 means 1/2 or 3/2
 * @param Z 0 or 1 means 0 or 1
 * @param S_ab 0 or 1 means 0 or 1
 * @param S_cd 0 or 1 means 0 or 1
 * @return the 6j symbol
 */
double Tools::g9j(int S,int Z,int S_ab,int S_cd){

   return x9j[8*S + 4*Z + 2*S_ab + S_cd];

}
