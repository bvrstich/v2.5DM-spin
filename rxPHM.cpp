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

int rxPHM::M;
int rxPHM::N;

vector< vector<int> > **rxPHM::ph2s;
int *****rxPHM::s2ph;

/**
 * initialize the static variables and allocate and initialize the static lists
 * @param M_in nr of sp orbs
 * @param N_in nr of particles
 */
void rxPHM::init(int M_in,int N_in){

   M = M_in;
   N = N_in;

   //allocate
   ph2s = new vector< vector<int> > * [M];

   for(int k = 0;k < M;++k)
      ph2s[k] = new vector< vector<int> > [2];

   s2ph = new int **** [M];

   for(int k = 0;k < M;++k){

      s2ph[k] = new int *** [2];

      for(int S = 0;S < 2;++S){

         s2ph[k][S] = new int ** [2];

         for(int S_bl = 0;S_bl < 2;++S_bl){

            s2ph[k][S][S_bl] = new int * [M];

            for(int a = 0;a < M;++a)
               s2ph[k][S][S_bl][a] = new int [M];

         }

      }

   }

   vector<int> v(3);

   int ph;

   for(int k = 0;k < M;++k){

      ph = 0;

      //S == 1/2

      //S_bl = 0;
      for(int a = 0;a < M;++a)
         for(int b = 0;b < M;++b){

            v[0] = 0;//S_bl
            v[1] = a;
            v[2] = b;

            ph2s[k][0].push_back(v);

            s2ph[k][0][0][a][b] = ph;

            ++ph;

         }

      //S_bl = 1
      for(int a = 0;a < M;++a)
         for(int b = 0;b < M;++b){

            if(b == k)
               b++;

            if(b == M)
               break;

            v[0] = 1;//S_bl
            v[1] = a;
            v[2] = b;

            ph2s[k][0].push_back(v);

            s2ph[k][0][1][a][b] = ph;

            ++ph;

         }

      //S == 3/2: only S_bl = 1
      ph = 0;

      for(int a = 0;a < M;++a)
         for(int b = 0;b < M;++b){

            if(b == k)
               b++;

            if(b == M)
               break;

            v[0] = 1;//S_bl
            v[1] = a;
            v[2] = b;

            ph2s[k][1].push_back(v);

            s2ph[k][1][1][a][b] = ph;

            ++ph;

         }

   }

}

/**
 * deallocate the static lists
 */
void rxPHM::clear(){

   for(int k = 0;k < M;++k)
      delete [] ph2s[k];

   delete [] ph2s;

   for(int k = 0;k < M;++k){

      for(int S = 0;S < 2;++S){

         for(int S_bl = 0;S_bl < 2;++S_bl){

            for(int a = 0;a < M;++a)
               delete [] s2ph[k][S][S_bl][a];

            delete [] s2ph[k][S][S_bl];

         }

         delete [] s2ph[k][S];

      }

      delete [] s2ph[k];

   }

   delete [] s2ph;

}

/**
 * standard constructor, constructs two blocks, the S = 1/2 and S = 3/2
 * @param l the parameter l that determines which spatial sp indices are blocked out
 */
rxPHM::rxPHM(int l) : BlockMatrix(2) {

   this->l = l;

   this->setMatrixDim(0,ph2s[l][0].size(),2);
   this->setMatrixDim(1,ph2s[l][1].size(),4);

}

/**
 * copy constructor: constructs Matrix object and fills it with the content of blockmatrix rxphm_c
 * @param rxphm_c object that will be copied into this.
 */
rxPHM::rxPHM(const rxPHM &rxphm_c) : BlockMatrix(rxphm_c){ 

   this->l = rxphm_c.gl();

}
/**
 * destructor
 */
rxPHM::~rxPHM(){ }

ostream &operator<<(ostream &output,const rxPHM &rxphm_p){

   for(int S = 0;S < 2;++S){

      output << endl;
      output << "S = " << 2*S + 1 << "/" << 2 << endl;
      output << endl;

      for(int i = 0;i < rxphm_p.gdim(S);++i)
         for(int j = 0;j < rxphm_p.gdim(S);++j){

            output << i << "\t" << j << "\t|\t(" << rxphm_p.ph2s[rxphm_p.gl()][S][i][0] << ")\t" << 
            
               rxphm_p.ph2s[rxphm_p.gl()][S][i][1] << "\t" << rxphm_p.ph2s[rxphm_p.gl()][S][i][2]

               << "\t;\t(" << rxphm_p.ph2s[rxphm_p.gl()][S][j][0] << ")\t" << rxphm_p.ph2s[rxphm_p.gl()][S][j][1]
               
               << "\t" << rxphm_p.ph2s[rxphm_p.gl()][S][j][2] << "\t" 

               << rxphm_p[S](i,j) << endl;

         }

   }

   return output;

}

/**
 * @return number of particles
 */
int rxPHM::gN() const{

   return N;

}

/**
 * @return number of sp orbitals
 */
int rxPHM::gM() const{

   return M;

}

/**
 * @return the parameter l
 */
int rxPHM::gl() const{

   return l;

}

/**
 * @param k which type is the rxPHM
 * @param S the dp spin
 * @param i the tp index
 * @param option == 0 return a, == 1 return b
 * @return the sp indices corresponding to the tp index i with parameter l
 */
int rxPHM::gph2s(int k,int S,int i,int option){

   return ph2s[k][S][i][option];

}

/**
 * @param k which type is the rxPHM
 * @param S the dp spin
 * @param S_ab intermediate spin
 * @param a the first sp index
 * @param b the second sp index
 * @return the tp index corresponding to the sp indices a and b with parameter l
 */
int rxPHM::gs2ph(int k,int S,int S_ab,int a,int b){

   return s2ph[k][S][S_ab][a][b];

}

/**
 * test if basis is correctly constructed
 */
void rxPHM::print_basis(){

   for(int k = 0;k < M;++k){

      cout << endl;
      cout << "l =\t" << k << endl;
      cout << endl;

      for(int S = 0;S < 2;++S){

         cout << endl;
         cout << "S =\t" << 2*S + 1 << "/" << 2 << endl;
         cout << endl;

         for(unsigned int i = 0;i < ph2s[k][S].size();++i)
            cout << i << "\t|\t(" << ph2s[k][S][i][0] << ")\t" << ph2s[k][S][i][1] << "\t" << ph2s[k][S][i][2] << endl;

      }

   }

}
