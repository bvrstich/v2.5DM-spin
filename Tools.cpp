#include <iostream>
#include <fstream>
#include <cmath>

#include "include.h"

dDPM *Tools::unit;
double **Tools::_6j;

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

   unit = new dDPM();

   *unit = 0.0;

   for(int l = 0;l < M;++l)
      for(int S = 0;S < 2;++S)
         for(int i = 0;i < (*unit)[l].gdim(S);++i)
            (*unit)[l](S,i,i) = 1.0;

   unit->proj_W();

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
 * deallocate the static lists and objects
 */
void Tools::clear(){

   delete unit;

   for(int S = 0;S < 2;++S)
      delete [] _6j[S];

   delete [] _6j;

}

/**
 * @param S_ab first parameter
 * @param S_cd second parameter
 * @return the 6j symbol with entries S_ab and S_cd
 */
double Tools::g6j(int S_ab,int S_cd){

   return _6j[S_ab][S_cd];

}

/**
 * @return the projected unit matrix
 */
dDPM &Tools::gunit(){

   return *unit;

}
