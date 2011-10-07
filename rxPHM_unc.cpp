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

int rxPHM_unc::M;
int rxPHM_unc::N;

vector< vector<int> > *rxPHM_unc::ph2s;
int ****rxPHM_unc::s2ph;

/**
 * initialize the static variables and allocate and initialize the static lists
 * @param M_in nr of sp orbs
 * @param N_in nr of particles
 */
void rxPHM_unc::init(int M_in,int N_in){

   M = M_in;
   N = N_in;

   //allocate
   ph2s = new vector< vector<int> > [M];

   s2ph = new int *** [M];

   for(int k = 0;k < M;++k){

      s2ph[k] = new int ** [2];

      for(int s = 0;s < 2;++s){

         s2ph[k][s] = new int * [2*M];

         for(int a = 0;a < 2*M;++a)
            s2ph[k][s][a] = new int [2*M];

      }

   }

   //initialize
   int i;

   vector<int> v(3);

   for(int k = 0;k < M;++k){

      i = 0;

      //0 is up, 1 is down
      for(int s = 0;s < 2;++s){

         for(int a = 0;a < 2*M;++a)
            for(int b = 0;b < 2*M;++b){

               if(b == 2*k + s)
                  b++;

               if(b == 2*M)
                  break;

               v[0] = s;
               v[1] = a;
               v[2] = b;

               ph2s[k].push_back(v);

               s2ph[k][s][a][b] = i;

               ++i;

            }

      }

   }

}

/**
 * deallocate the static lists
 */
void rxPHM_unc::clear(){

   for(int k = 0;k < M;++k){

      for(int s = 0;s < 2;++s){

         for(int a = 0;a < 2*M;++a)
            delete [] s2ph[k][s][a];

         delete [] s2ph[k][s];

      }

      delete [] s2ph[k];

   }

   delete [] s2ph;

   delete [] ph2s;

}

/**
 * standard constructor
 * @param l the parameter l that determines which sp index is blocked out
 */
rxPHM_unc::rxPHM_unc(int l) : Matrix(ph2s[l].size()) {

   this->l = l;

}

/**
 * copy constructor: constructs Matrix object and fills it with the content of matrix rxphm_c
 * @param rxphm_c object that will be copied into this.
 */
rxPHM_unc::rxPHM_unc(const rxPHM_unc &rxphm_c) : Matrix(rxphm_c){ 

   this->l = rxphm_c.gl();

}

/**
 * destructor
 * 
 */
rxPHM_unc::~rxPHM_unc(){ }

/**
 * access the elements of the matrix in sp mode, antisymmetry is automatically accounted for:
 * @param s_ab spinindex for the row
 * @param a first sp index that forms the tp row index i together with b
 * @param b second sp index that forms the tp row index i together with a
 * @param c first sp index that forms the tp column index j together with d
 * @param s_cd spinindex for the column
 * @param d second sp index that forms the tp column index j together with c
 * @return the number on place rxPHM_unc(i,j) with the right phase.
 */
double rxPHM_unc::operator()(int s_ab,int a,int b,int s_cd,int c,int d) const{

   if(b == 2*l + s_ab)
      return 0.0;

   if(d == 2*l + s_cd)
      return 0.0;

   return (*this)(s2ph[l][s_ab][a][b],s2ph[l][s_cd][c][d]);

}

ostream &operator<<(ostream &output,const rxPHM_unc &rxphm_p){

   for(int i = 0;i < rxphm_p.gn();++i)
      for(int j = 0;j < rxphm_p.gn();++j){

         output << i << "\t" << j << "\t|\t" << rxphm_p.ph2s[rxphm_p.l][i][0] << "\t" << rxphm_p.ph2s[rxphm_p.l][i][1] << "\t" << rxphm_p.ph2s[rxphm_p.l][i][2]

            << "\t" << rxphm_p.ph2s[rxphm_p.l][j][0] << "\t" << rxphm_p.ph2s[rxphm_p.l][j][1] << "\t" << rxphm_p.ph2s[rxphm_p.l][j][2]

            << "\t" << rxphm_p(i,j) << endl;

      }

   return output;

}

/**
 * @return number of particles
 */
int rxPHM_unc::gN() const{

   return N;

}

/**
 * @return number of sp orbitals
 */
int rxPHM_unc::gM() const{

   return M;

}

/**
 * @return the parameter l
 */
int rxPHM_unc::gl() const{

   return l;

}

/**
 * @param k which type is the rxPHM_unc
 * @param i the tp index
 * @param option ==0 return s_ab, == 1 return a, == 2 return b
 * @return the sp indices corresponding to the tp index i with parameter k
 */
int rxPHM_unc::gph2s(int k,int i,int option){

   return ph2s[k][i][option];

}

/**
 * @param k which type is the rxPHM_unc
 * @param s_ab the spinindex
 * @param a the first sp index
 * @param b the second sp index
 * @return the tp index corresponding to the sp indices a and b with parameter k
 */
int rxPHM_unc::gs2ph(int k,int s_ab,int a,int b){

   return s2ph[k][s_ab][a][b];

}

/**
 * test function that prints out the basis list
 */
void rxPHM_unc::print_basis(){

   for(int k = 0;k < M;++k){

      cout << endl;
      cout << "l = " << k << endl;
      cout << endl;

      for(unsigned int i = 0;i < ph2s[k].size();++i)
         cout << i << "\t|\t" << ph2s[k][i][0] << "\t(" << ph2s[k][i][1]/2 << "," << ph2s[k][i][1]%2 << ")\t(" << ph2s[k][i][2]/2 << "," << ph2s[k][i][2]%2 << ")" << endl;

   }

}
