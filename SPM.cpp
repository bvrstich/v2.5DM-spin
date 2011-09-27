#include <iostream>
#include <fstream>
#include <cmath>

using std::ostream;
using std::endl;

#include "include.h"

int SPM::M;
int SPM::N;

/**
 * initalize the static variables
 */
void SPM::init(int M_in,int N_in){

   M = M_in;
   N = N_in;

}

/**
 * constructor, makes matrix of dimension M, there will be two degenerate blocks +1/2,-1/2. So
 * only one matrix is needed of dimension M is needed to represent the SPM.
 * @param M dimension of single particle space
 * @param N Nr of particles
 */
SPM::SPM() : Matrix(M) { }

/**
 * copy constructor
 * @param spm_copy content of this matrix will be copied into the constructed matrix
 */
SPM::SPM(const SPM &spm_copy) : Matrix(spm_copy) { }

/**
 * TPM constructor: Creates a SPM initialized on the "bar" of the TPM.
 * @param scale the factor u want the SPM to be scaled with (1/N-1 for normal sp density matrix)
 * @param tpm the TPM out of which the SPM will be initiated.
 */
SPM::SPM(double scale,const TPM &tpm) : Matrix(tpm.gM()) {

   this->bar(scale,tpm);

}

/**
 * destructor
 */
SPM::~SPM(){

}

/**
 * @return nr of particles
 */
int SPM::gN() const{

   return N;

}

/**
 * @return dimension of sp space
 */
int SPM::gM() const{

   return M;

}

ostream &operator<<(ostream &output,const SPM &spm_p){

   for(int i = 0;i < spm_p.gn();++i)
      for(int j = 0;j < spm_p.gn();++j)
         output << i << "\t" << j << "\t" << spm_p(i,j) << endl;

   return output;

}

/**
 * Trace out a set of indices to create the "bar" matrix of a TPM
 * @param scale the factor u want the SPM to be scaled with (1/N-1 for normal sp density matrix)
 * @param tpm the TPM out of which the SPM will be filled
 */
void SPM::bar(double scale,const TPM &tpm){

   //hulpvariabele
   double ward;

   for(int a = 0;a < M;++a)
      for(int c = a;c < M;++c){

         (*this)(a,c) = 0.0;

         for(int b = 0;b < M;++b){

            //S = 0 stuk
            ward = tpm(0,a,b,c,b);

            if(a == b)
               ward *= std::sqrt(2.0);

            if(c == b)
               ward *= std::sqrt(2.0);

            (*this)(a,c) += ward;

            //S = 1 stuk: hier kan nooit a = b en c = d wegens antisymmetrie
            (*this)(a,c) += 3.0*tpm(1,a,b,c,b);

         }

         //nog schalen
         (*this)(a,c) *= 0.5*scale;

      }

   this->symmetrize();

}

/**
 * map a dPPHM matrix on a SPM by tracing the second indices together with the diagonal third (see symmetry.pdf)
 * @param scale the SPM with this number
 * @param dpphm input dPPHM object
 */
void SPM::breve(double scale,const dDPM &ddpm){

   double ward,hard;

   for(int b = 0;b < M;++b)
      for(int d = b;d < M;++d){

         (*this)(b,d) = 0.0;

         //first S = 1/2
         for(int S_ab = 0;S_ab < 2;++S_ab)
            for(int S_cd = 0;S_cd < 2;++S_cd){

               ward = 0.0;

               for(int l = 0;l < M;++l){

                  hard = ddpm(l,0,S_ab,l,b,S_cd,l,d);

                  if(b == l)
                     hard *= std::sqrt(2.0);

                  if(d == l)
                     hard *= std::sqrt(2.0);

                  ward += hard;

               }

               (*this)(b,d) += std::sqrt( (2.0*S_ab + 1.0) * (2.0*S_cd + 1.0) ) * Tools::g6j(S_ab,S_cd) * ward;

            }

         //then S = 3/2
         for(int l = 0;l < M;++l)
            (*this)(b,d) -= 2.0 * ddpm(l,1,1,l,b,1,l,d);

         (*this)(b,d) *= scale;

      }

   this->symmetrize();

}
