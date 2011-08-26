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

int rTPM::M;
int rTPM::N;

vector< vector<int> > **rTPM::t2s;
int *****rTPM::s2t;

/**
 * initialize the static variables and allocate and initialize the static lists
 * @param M_in nr of sp orbs
 * @param N_in nr of particles
 */
void rTPM::init(int M_in,int N_in){

   M = M_in;
   N = N_in;

   //allocate
   t2s = new vector< vector<int> > * [M];

   for(int k = 0;k < M;++k)
      t2s[k] = new vector< vector<int> > [2];

   s2t = new int **** [M];

   for(int k = 0;k < M;++k){

      s2t[k] = new int *** [2];

      for(int S = 0;S < 2;++S){

         s2t[k][S] = new int ** [2];

         for(int S_ab = 0;S_ab < 2;++S_ab){

            s2t[k][S][S_ab] = new int * [M];

            for(int a = 0;a < M;++a)
               s2t[k][S][S_ab][a] = new int [M];

         }

      }

   }

   vector<int> v(3);

   int tp;

   for(int k = 0;k < M;++k){

      tp = 0;

      //S == 1/2
      for(int S_ab = 0;S_ab < 2;++S_ab){

         for(int a = 0;a < M;++a){

            if(a == k)
               ++a;

            if(a == M)
               break;

            for(int b = a + S_ab;b < M;++b){

               if(b == k)
                  ++b;

               if(b == M)
                  break;

               v[0] = S_ab;//S_ab
               v[1] = a;
               v[2] = b;

               t2s[k][0].push_back(v);

               s2t[k][0][S_ab][a][b] = tp;
               s2t[k][0][S_ab][b][a] = tp;

               ++tp;

            }

         }

      }

      tp = 0;

      //S == 3/2
      for(int a = 0;a < M;++a){

         if(a == k)
            ++a;

         if(a == M)
            break;

         for(int b = a + 1;b < M;++b){

            if(b == k)
               ++b;

            if(b == M)
               break;

            v[0] = 1;//S_ab
            v[1] = a;
            v[2] = b;

            t2s[k][1].push_back(v);

            s2t[k][1][1][a][b] = tp;
            s2t[k][1][1][b][a] = tp;

            ++tp;

         }
      }

   }

}

/**
 * deallocate the static lists
 */
void rTPM::clear(){

   for(int k = 0;k < M;++k)
      delete [] t2s[k];

   delete [] t2s;

   for(int k = 0;k < M;++k){

      for(int S = 0;S < 2;++S){

         for(int S_ab = 0;S_ab < 2;++S_ab){

            for(int a = 0;a < M;++a)
               delete [] s2t[k][S][S_ab][a];

            delete [] s2t[k][S][S_ab];

         }

         delete [] s2t[k][S];

      }

      delete [] s2t[k];

   }

   delete [] s2t;

}

/**
 * standard constructor, constructs two blocks, the S = 1/2 and S = 3/2
 * @param l the parameter l that determines which spatial sp indices are blocked out
 */
rTPM::rTPM(int l) : BlockMatrix(2) {

   this->l = l;

   this->setMatrixDim(0,t2s[l][0].size(),2);
   this->setMatrixDim(1,t2s[l][1].size(),4);

}

/**
 * copy constructor: constructs Matrix object and fills it with the content of blockmatrix rtpm_c
 * @param rtpm_c object that will be copied into this.
 */
rTPM::rTPM(const rTPM &rtpm_c) : BlockMatrix(rtpm_c){ 

   this->l = rtpm_c.gl();

}
/**
 * destructor
 */
rTPM::~rTPM(){ }

ostream &operator<<(ostream &output,const rTPM &rtpm_p){

   for(int S = 0;S < 2;++S){

      output << endl;
      output << "S = " << 2*S + 1 << "/" << 2 << endl;
      output << endl;

      for(int i = 0;i < rtpm_p.gdim(S);++i)
         for(int j = 0;j < rtpm_p.gdim(S);++j){

            output << i << "\t" << j << "\t|\t(" << rtpm_p.t2s[rtpm_p.gl()][S][i][0] << ")\t" << rtpm_p.t2s[rtpm_p.gl()][S][i][1] << "\t" << rtpm_p.t2s[rtpm_p.gl()][S][i][2]

               << "\t;\t(" << rtpm_p.t2s[rtpm_p.gl()][S][j][0] << ")\t" << rtpm_p.t2s[rtpm_p.gl()][S][j][1] << "\t" << rtpm_p.t2s[rtpm_p.gl()][S][j][2] << "\t" 
               
               << rtpm_p[S](i,j) << endl;

         }

   }

   return output;

}

/**
 * @return number of particles
 */
int rTPM::gN() const{

   return N;

}

/**
 * @return number of sp orbitals
 */
int rTPM::gM() const{

   return M;

}

/**
 * @return the parameter l
 */
int rTPM::gl() const{

   return l;

}

/**
 * @param k which type is the rTPM
 * @pararm S the dp spin
 * @param i the tp index
 * @param option == 0 return a, == 1 return b
 * @return the sp indices corresponding to the tp index i with parameter l
 */
int rTPM::gt2s(int k,int S,int i,int option){

   return t2s[k][S][i][option];

}
 
/**
 * @param k which type is the rTPM
 * @param S the dp spin
 * @pararm S_ab intermediate spin
 * @param a the first sp index
 * @param b the second sp index
 * @return the tp index corresponding to the sp indices a and b with parameter l
 */
int rTPM::gs2t(int k,int S,int S_ab,int a,int b){

   return s2t[k][S][S_ab][a][b];

}
