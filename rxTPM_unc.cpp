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

int rxTPM_unc::M;
int rxTPM_unc::N;

vector< vector<int> > *rxTPM_unc::t2s;
int ****rxTPM_unc::s2t;

/**
 * initialize the static variables and allocate and initialize the static lists
 * @param M_in nr of sp orbs
 * @param N_in nr of particles
 */
void rxTPM_unc::init(int M_in,int N_in){

   M = M_in;
   N = N_in;

   //allocate
   t2s = new vector< vector<int> > [M];

   s2t = new int *** [M];

   for(int k = 0;k < M;++k){

      s2t[k] = new int ** [2];

      for(int s = 0;s < 2;++s){

         s2t[k][s] = new int * [2*M];

         for(int a = 0;a < 2*M;++a)
            s2t[k][s][a] = new int [2*M];

      }

   }

   //initialize
   int i;

   vector<int> v(3);

   for(int k = 0;k < M;++k){

      i = 0;

      //0 is up, 1 is down
      for(int s = 0;s < 2;++s){

         for(int a = 0;a < 2*M;++a){

            if(a == 2*k + s)
               a++;

            if(a == 2*M)
               break;

            for(int b = a + 1;b < 2*M;++b){

               if(b == 2*k + s)
                  b++;

               if(b == 2*M)
                  break;

               v[0] = s;
               v[1] = a;
               v[2] = b;

               t2s[k].push_back(v);

               s2t[k][s][a][b] = i;
               s2t[k][s][b][a] = i;

               ++i;

            }

         }

      }

   }

}

/**
 * deallocate the static lists
 */
void rxTPM_unc::clear(){

   for(int k = 0;k < M;++k){

      for(int s = 0;s < 2;++s){

         for(int a = 0;a < 2*M;++a)
            delete [] s2t[k][s][a];

         delete [] s2t[k][s];

      }

      delete [] s2t[k];

   }

   delete [] s2t;

   delete [] t2s;

}

/**
 * standard constructor
 * @param l the parameter l that determines which sp index is blocked out
 */
rxTPM_unc::rxTPM_unc(int l) : Matrix(t2s[l].size()) {

   this->l = l;

}

/**
 * copy constructor: constructs Matrix object and fills it with the content of matrix rxtpm_c
 * @param rxtpm_c object that will be copied into this.
 */
rxTPM_unc::rxTPM_unc(const rxTPM_unc &rxtpm_c) : Matrix(rxtpm_c){ 

   this->l = rxtpm_c.gl();

}

/**
 * destructor
 * 
 */
rxTPM_unc::~rxTPM_unc(){ }

/**
 * access the elements of the matrix in sp mode, antisymmetry is automatically accounted for:
 * @param s_ab spinindex for the row
 * @param a first sp index that forms the tp row index i together with b
 * @param b second sp index that forms the tp row index i together with a
 * @param c first sp index that forms the tp column index j together with d
 * @param s_cd spinindex for the column
 * @param d second sp index that forms the tp column index j together with c
 * @return the number on place rxTPM_unc(i,j) with the right phase.
 */
double rxTPM_unc::operator()(int s_ab,int a,int b,int s_cd,int c,int d) const{

   int phase_i = get_inco(l,s_ab,a,b);

   if(phase_i == 0)
      return 0.0;

   int phase_j = get_inco(l,s_cd,c,d);

   if(phase_j == 0)
      return 0.0;

   return phase_i * phase_j * (*this)(s2t[l][s_ab][a][b],s2t[l][s_cd][c][d]);

}

/**
 * order the dDP basis, switch the indices to get the symmetry-part I keep in my basis and add the correct phase.
 * @param k parameter that indicates which spatial indices are blocked out
 * @param s third index spin
 * @param a first index
 * @param b second index
 */
int rxTPM_unc::get_inco(int k,int &s,int &a,int &b){

   if(a == b)
      return 0;

   if(2*k + s == a || 2*k + s == b)
      return 0;

   if(a > b)
      return -1;
   else
      return 1;

}

ostream &operator<<(ostream &output,const rxTPM_unc &rxtpm_p){

   for(int i = 0;i < rxtpm_p.gn();++i)
      for(int j = 0;j < rxtpm_p.gn();++j){

         output << i << "\t" << j << "\t|\t" << rxtpm_p.t2s[rxtpm_p.l][i][0] << "\t" << rxtpm_p.t2s[rxtpm_p.l][i][1] << "\t" << rxtpm_p.t2s[rxtpm_p.l][i][2]

            << "\t" << rxtpm_p.t2s[rxtpm_p.l][j][0] << "\t" << rxtpm_p.t2s[rxtpm_p.l][j][1] << "\t" << rxtpm_p.t2s[rxtpm_p.l][j][2]

            << "\t" << rxtpm_p(i,j) << endl;

      }

   return output;

}

/**
 * @return number of particles
 */
int rxTPM_unc::gN() const{

   return N;

}

/**
 * @return number of sp orbitals
 */
int rxTPM_unc::gM() const{

   return M;

}

/**
 * @return the parameter l
 */
int rxTPM_unc::gl() const{

   return l;

}

/**
 * @param k which type is the rxTPM_unc
 * @param i the tp index
 * @param option ==0 return s_ab, == 1 return a, == 2 return b
 * @return the sp indices corresponding to the tp index i with parameter k
 */
int rxTPM_unc::gt2s(int k,int i,int option){

   return t2s[k][i][option];

}

/**
 * @param k which type is the rxTPM_unc
 * @param s_ab the spinindex
 * @param a the first sp index
 * @param b the second sp index
 * @return the tp index corresponding to the sp indices a and b with parameter k
 */
int rxTPM_unc::gs2t(int k,int s_ab,int a,int b){

   return s2t[k][s_ab][a][b];

}

/**
 * test function that prints out the basis list
 */
void rxTPM_unc::print_basis(){

   for(int k = 0;k < M;++k){

      cout << endl;
      cout << "l = " << k << endl;
      cout << endl;

      for(unsigned int i = 0;i < t2s[k].size();++i)
         cout << i << "\t|\t" << t2s[k][i][0] << "\t(" << t2s[k][i][1]/2 << "," << t2s[k][i][1]%2 << ")\t(" << t2s[k][i][2]/2 << "," << t2s[k][i][2]%2 << ")" << endl;

   }

}
