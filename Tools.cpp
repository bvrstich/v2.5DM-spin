#include <iostream>
#include <fstream>
#include <cmath>

#include "include.h"

using std::cout;
using std::endl;

dDPM *Tools::unit;

double ******Tools::_6j;
double ******Tools::_3j;

double ****Tools::C;

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
   _6j = new double ***** [5];

   for(int j1 = 0;j1 < 5;++j1){

      _6j[j1] = new double **** [5];

      for(int j2 = 0;j2 < 5;++j2){

         _6j[j1][j2] = new double *** [5];

         for(int j3 = 0;j3 < 5;++j3){

            _6j[j1][j2][j3] = new double ** [5];

            for(int j4 = 0;j4 < 5;++j4){

               _6j[j1][j2][j3][j4] = new double * [5];

               for(int j5 = 0;j5 < 5;++j5)
                  _6j[j1][j2][j3][j4][j5] = new double [5];

            }
         }
      }
   }

   //allocate
   _3j = new double ***** [4];

   for(int j1 = 0;j1 < 4;++j1){

      _3j[j1] = new double **** [4];

      for(int j2 = 0;j2 < 4;++j2){

         _3j[j1][j2] = new double *** [4];

         for(int j3 = 0;j3 < 4;++j3){

            _3j[j1][j2][j3] = new double ** [7];

            for(int m1 = 0;m1 < 7;++m1){

               _3j[j1][j2][j3][m1] = new double * [7];

               for(int m2 = 0;m2 < 7;++m2)
                  _3j[j1][j2][j3][m1][m2] = new double [7];

            }
         }
      }
   }

   C = new double *** [2];

   for(int S = 0;S < 2;++S){

      C[S] = new double ** [2];

      for(int S_ = 0;S_ < 2;++S_){

         C[S][S_] = new double * [2];

         for(int S_ab = 0;S_ab < 2;++S_ab)
            C[S][S_][S_ab] = new double [2];

      }
   }

   //fill the 6j list
   for(long int j1 = 0;j1 < 5;++j1)
      for(long int j2 = 0;j2 < 5;++j2)
         for(long int j3 = 0;j3 < 5;++j3)
            for(long int j4 = 0;j4 < 5;++j4)
               for(long int j5 = 0;j5 < 5;++j5)
                  for(long int j6 = 0;j6 < 5;++j6)
                     _6j[j1][j2][j3][j4][j5][j6] = x6j_(&j1,&j2,&j3,&j4,&j5,&j6);

   //fill the 3j list
   for(long int j1 = 0;j1 < 4;++j1)
      for(long int j2 = 0;j2 < 4;++j2)
         for(long int j3 = 0;j3 < 4;++j3)
            for(long int m1 = 0;m1 < 7;++m1)
               for(long int m2 = 0;m2 < 7;++m2)
                  for(long int m3 = 0;m3 < 7;++m3){

                     long int M1 = m1 - 3;
                     long int M2 = m2 - 3;
                     long int M3 = m3 - 3;

                     _3j[j1][j2][j3][m1][m2][m3] = x3j_(&j1,&j2,&j3,&M1,&M2,&M3);

                  }
   
   for(int S = 0;S < 2;++S)
      for(int S_ = 0;S_ < 2;++S_)
         for(int S_ab = 0;S_ab < 2;++S_ab)
            for(int S_cd = 0;S_cd < 2;++S_cd){

               C[S][S_][S_ab][S_cd] = 0.0;

               for(int j = 0;j < 3;++j)
                  C[S][S_][S_ab][S_cd] += (2*j + 1.0) * (1 - 2*(j%2)) * g6j(2*S_ab,2*S_cd,2*j,2*S_ + 1,2*S + 1,1) * g6j(2*S_ab,2*S_cd,2*j,2*S + 1,2*S_ + 1,1);

            }

}

/**
 * deallocate the static lists and objects
 */
void Tools::clear(){

   delete unit;

   //delete the 6j list
   for(long int j1 = 0;j1 < 4;++j1){

      for(long int j2 = 0;j2 < 4;++j2){

         for(long int j3 = 0;j3 < 4;++j3){

            for(long int j4 = 0;j4 < 4;++j4){

               for(long int j5 = 0;j5 < 4;++j5)
                  delete [] _6j[j1][j2][j3][j4][j5];

               delete [] _6j[j1][j2][j3][j4];

            }

            delete [] _6j[j1][j2][j3];

         }

         delete [] _6j[j1][j2];

      }

      delete [] _6j[j1];

   }

   delete [] _6j;

   //delete the 3j list
   for(long int j1 = 0;j1 < 4;++j1){

      for(long int j2 = 0;j2 < 4;++j2){

         for(long int j3 = 0;j3 < 4;++j3){

            for(long int m1 = 0;m1 < 7;++m1){

               for(long int m2 = 0;m2 < 7;++m2)
                  delete [] _3j[j1][j2][j3][m1][m2];

               delete [] _3j[j1][j2][j3][m1];

            }

            delete [] _3j[j1][j2][j3];

         }

         delete [] _3j[j1][j2];

      }

      delete [] _3j[j1];

   }

   delete [] _3j;

   for(int S = 0;S < 2;++S){

      for(int S_ = 0;S_ < 2;++S_){

         for(int S_ab = 0;S_ab < 2;++S_ab)
            delete [] C[S][S_][S_ab];

         delete [] C[S][S_];

      }

      delete [] C[S];

   }

   delete [] C;

}

/**
 * watch out, the j's you have to enter are 2 times the physical angular momentum! (i.e. S = 1/2 means you enter j = 1)
 * @param j1 first j
 * @param j2 second j
 * @param j3 third j
 * @param j4 fourth j
 * @param j5 fifth j
 * @param j6 sixt j
 * @return the 6j symbol
 */
double Tools::g6j(int j1,int j2,int j3,int j4,int j5,int j6){

   return _6j[j1][j2][j3][j4][j5][j6];

}

/**
 * watch out, the j's and m's you have to enter are 2 times the physical angular momentum! (i.e. S = 1/2 means you enter j = 1)
 * m's can be negative!
 * @param j1 first j
 * @param j2 second j
 * @param j3 third j
 * @param m1 first m
 * @param m2 second m
 * @param m3 third m
 * @return the 3j symbol
 */
double Tools::g3j(int j1,int j2,int j3,int m1,int m2,int m3){

   return _3j[j1][j2][j3][m1 + 3][m2 + 3][m3 + 3];

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

   return C[S][S_][S_ab][S_cd];

}
