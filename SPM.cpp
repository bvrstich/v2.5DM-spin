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
 * @param M_in input nr of sites
 * @param N_in input nr of particles
 */
void SPM::init(int M_in,int N_in){

   M = M_in;
   N = N_in;

}

/**
 * constructor, makes matrix of dimension M, there will be two degenerate blocks +1/2,-1/2. So
 * only one matrix is needed of dimension M is needed to represent the SPM.
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
 * map a dDPM matrix on a SPM by tracing the first indices together with the diagonal third (see symmetry.pdf)
 * @param scale the SPM with this number
 * @param ddpm input dDPM object
 */
void SPM::breve(double scale,const dDPM &ddpm){

   double ward,hard;

   for(int b = 0;b < M;++b)
      for(int d = b;d < M;++d){

         (*this)(b,d) = 0.0;

         //only S = 1/2 contribution
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

               (*this)(b,d) += std::sqrt( (2.0*S_ab + 1.0) * (2.0*S_cd + 1.0) ) * Tools::g6j(0,0,S_ab,S_cd) * ward;

            }

         (*this)(b,d) *= scale;

      }

   this->symmetrize();

}

/**
 * map a dPPHM matrix on a SPM by tracing the first indices together with the diagonal third (see symmetry.pdf)
 * @param scale the SPM with this number
 * @param dpphm input dPPHM object
 */
void SPM::breve(double scale,const dPPHM &dpphm){

   double ward,hard;

   for(int b = 0;b < M;++b)
      for(int d = b;d < M;++d){

         (*this)(b,d) = 0.0;

         //only S = 1/2 contribution
         for(int S_ab = 0;S_ab < 2;++S_ab)
            for(int S_cd = 0;S_cd < 2;++S_cd){

               ward = 0.0;

               for(int l = 0;l < M;++l){

                  hard = dpphm(l,0,S_ab,l,b,S_cd,l,d);

                  if(b == l)
                     hard *= std::sqrt(2.0);

                  if(d == l)
                     hard *= std::sqrt(2.0);

                  ward += hard;

               }

               (*this)(b,d) += std::sqrt( (2.0*S_ab + 1.0) * (2.0*S_cd + 1.0) ) * ward;

            }

         (*this)(b,d) *= scale;

      }

   this->symmetrize();

}

/**
 * map a dPHHM matrix on a SPM by tracing the first indices together with the diagonal third (see symmetry.pdf)
 * @param scale the SPM with this number
 * @param dphhm input dPHHM object
 */
void SPM::breve(double scale,const dPHHM &dphhm){

   double ward;

   for(int b = 0;b < M;++b)
      for(int d = b;d < M;++d){

         (*this)(b,d) = 0.0;

         //only S = 1/2 contribution
         for(int S_ab = 0;S_ab < 2;++S_ab)
            for(int S_cd = 0;S_cd < 2;++S_cd){

               ward = 0.0;

               for(int l = 0;l < M;++l)
                  ward += dphhm(l,0,S_ab,l,b,S_cd,l,d);

               (*this)(b,d) += std::sqrt( (2.0*S_ab + 1.0) * (2.0*S_cd + 1.0) ) * ward;

            }

         (*this)(b,d) *= scale;

      }

   this->symmetrize();

}

/**
 * map a dPHHM matrix on a SPM by tracing the second indices together with the diagonal third (see symmetry.pdf)
 * this is for the G2 down, slightly different breve then for the G1 down, the trace is over the second index, so breve_si...
 * @param scale the SPM with this number
 * @param dphhm input dPHHM object
 */
void SPM::breve_si(double scale,const dPHHM &dphhm){

   double ward;

   for(int a = 0;a < M;++a)
      for(int c = a;c < M;++c){

         (*this)(a,c) = 0.0;

         for(int S = 0;S < 2;++S)
            for(int S_bl = 0;S_bl < 2;++S_bl){

               ward = 0.0;

               for(int l = 0;l < M;++l)
                  ward += dphhm(l,S,S_bl,a,l,S_bl,c,l);

               (*this)(a,c) += (1 - 2*S_bl) * (2*(S + 0.5) + 1.0) * ward;

            }

         (*this)(a,c) *= scale;

      }

   this->symmetrize();

}

/**
 * map a dPHHM matrix on a SPM by doing a normal double bar
 * @param scale the SPM with this number
 * @param dphhm input dPHHM object
 */
void SPM::barbar(double scale,const dPHHM &dphhm){

   for(int a = 0;a < M;++a)
      for(int c = a;c < M;++c){

         (*this)(a,c) = 0.0;

         for(int l = 0;l < M;++l)
            for(int S = 0;S < 2;++S)
               for(int S_bl = 0;S_bl < 2;++S_bl)
                  for(int b = 0;b < M;++b)
                     (*this)(a,c) += (2*(S + 0.5) + 1.0) * dphhm(l,S,S_bl,a,b,S_bl,c,b);
         
         (*this)(a,c) *= scale;

      }

   this->symmetrize();

}
