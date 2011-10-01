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

int xTPM_unc::M;
int xTPM_unc::N;

vector< vector<int> > xTPM_unc::t2s;
int ***xTPM_unc::s2t;

/**
 * initialize the static variables and allocate and initialize the static lists
 * @param M_in nr of sp orbs
 * @param N_in nr of particles
 */
void xTPM_unc::init(int M_in,int N_in){

   M = M_in;
   N = N_in;


   s2t = new int ** [2];

   for(int s = 0;s < 2;++s){

      s2t[s] = new int * [2*M];

      for(int a = 0;a < 2*M;++a)
         s2t[s][a] = new int [2*M];

   }

   //initialize
   int i;

   vector<int> v(3);

   i = 0;

   //0 is up, 1 is down
   for(int s = 0;s < 2;++s){

      for(int a = 0;a < 2*M;++a)
         for(int b = a + 1;b < 2*M;++b){

            v[0] = s;
            v[1] = a;
            v[2] = b;

            t2s.push_back(v);

            s2t[s][a][b] = i;
            s2t[s][b][a] = i;

            ++i;

         }

   }

}

/**
 * deallocate the static lists
 */
void xTPM_unc::clear(){

   for(int s = 0;s < 2;++s){

      for(int a = 0;a < 2*M;++a)
         delete [] s2t[s][a];

      delete [] s2t[s];

   }

   delete [] s2t;

}

/**
 * standard constructor
 * @param l the parameter l that determines which sp index is blocked out
 */
xTPM_unc::xTPM_unc() : Matrix(t2s.size()) { }

/**
 * copy constructor: constructs Matrix object and fills it with the content of matrix rxtpm_c
 * @param rxtpm_c object that will be copied into this.
 */
xTPM_unc::xTPM_unc(const xTPM_unc &rxtpm_c) : Matrix(rxtpm_c){ 

}

/**
 * destructor
 * 
 */
xTPM_unc::~xTPM_unc(){ }

/**
 * access the elements of the matrix in sp mode, antisymmetry is automatically accounted for:
 * @param s_ab spinindex for the row
 * @param a first sp index that forms the tp row index i together with b
 * @param b second sp index that forms the tp row index i together with a
 * @param c first sp index that forms the tp column index j together with d
 * @param s_cd spinindex for the column
 * @param d second sp index that forms the tp column index j together with c
 * @return the number on place xTPM_unc(i,j) with the right phase.
 */
double xTPM_unc::operator()(int s_ab,int a,int b,int s_cd,int c,int d) const{

   int phase_i = get_inco(s_ab,a,b);

   if(phase_i == 0)
      return 0.0;

   int phase_j = get_inco(s_cd,c,d);

   if(phase_j == 0)
      return 0.0;

   return phase_i * phase_j * (*this)(s2t[s_ab][a][b],s2t[s_cd][c][d]);

}

/**
 * order the dDP basis, switch the indices to get the symmetry-part I keep in my basis and add the correct phase.
 * @param s third index spin
 * @param a first index
 * @param b second index
 */
int xTPM_unc::get_inco(int s,int a,int b){

   if(a == b)
      return 0;

   if(a > b)
      return -1;
   else
      return 1;

}

ostream &operator<<(ostream &output,const xTPM_unc &rxtpm_p){

   for(int i = 0;i < rxtpm_p.gn();++i)
      for(int j = 0;j < rxtpm_p.gn();++j){

         output << i << "\t" << j << "\t|\t" << rxtpm_p.t2s[i][0] << "\t" << rxtpm_p.t2s[i][1] << "\t" << rxtpm_p.t2s[i][2]

            << "\t" << rxtpm_p.t2s[j][0] << "\t" << rxtpm_p.t2s[j][1] << "\t" << rxtpm_p.t2s[j][2] << "\t" << rxtpm_p(i,j) << endl;

      }

   return output;

}

/**
 * @return number of particles
 */
int xTPM_unc::gN() const{

   return N;

}

/**
 * @return number of sp orbitals
 */
int xTPM_unc::gM() const{

   return M;

}

/**
 * @param i the tp index
 * @param option ==0 return s_ab, == 1 return a, == 2 return b
 * @return the sp indices corresponding to the tp index i with parameter k
 */
int xTPM_unc::gt2s(int i,int option){

   return t2s[i][option];

}

/**
 * @param s_ab the spinindex
 * @param a the first sp index
 * @param b the second sp index
 * @return the tp index corresponding to the sp indices a and b with parameter k
 */
int xTPM_unc::gs2t(int s_ab,int a,int b){

   return s2t[s_ab][a][b];

}

/**
 * test function that prints out the basis list
 */
void xTPM_unc::print_basis(){

   for(unsigned int i = 0;i < t2s.size();++i)
      cout << i << "\t|\t" << t2s[i][0] << "\t(" << t2s[i][1]/2 << "," << t2s[i][1]%2 << ")\t(" << t2s[i][2]/2 << "," << t2s[i][2]%2 << ")" << endl;

}
