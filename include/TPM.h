#ifndef TPM_H
#define TPM_H

#include <iostream>
#include <fstream>

using std::ostream;

#include "BlockMatrix.h"

class dDPM;
class dPPHM;

/**
 * @author Brecht Verstichel
 * @date 19-04-2010\n\n
 * This class TPM is a class written for two particle matrices with spinsymmetry included, it inherits alle the function from its mother 
 * BlockMatrix, some special member functions and two lists that give the relationship between the sp and the tp 
 * basis.
 */
class TPM : public BlockMatrix {

   /**
    * Output stream operator overloaded, the usage is simple, if you want to print to a file, make an
    * ifstream object and type:\n\n
    * object << tpm_p << endl;\n\n
    * For output onto the screen type: \n\n
    * cout << tpm_p << endl;\n\n
    * @param output The stream to which you are writing (e.g. cout)
    * @param tpm_p the TPM you want to print
    */
   friend ostream &operator<<(ostream &output,const TPM &tpm_p);

   public:
      
      TPM();

      //copy constructor
      TPM(const TPM &);

      //destructor
      virtual ~TPM();

      void constr_lists();

      using BlockMatrix::operator=;

      using BlockMatrix::operator();

      //easy to access the numbers, in sp mode and with spin quantumnumer
      double operator()(int S,int a,int b,int c,int d) const;

      int gN() const;

      int gM() const;

      void hubbard(double U);

      void Q(int option,const TPM &);

      void Q(int option,double A,double B,double C,const TPM &);

      void set_unit();

      void set_S_2();

      double spin() const;

      void bar(double scale,const dDPM &);

      void bar(double scale,const dPPHM &);

      static void test_basis();

      static void init(int,int);

      static void clear();

   private:

      //!static list of dimension [2][dim[i]][2] that takes in a tp index i and a spinquantumnumber S, and returns two sp indices: a = t2s[S][i][0] and b = t2s[S][i][1]
      static vector< vector<int> > *t2s;

      //!static list of dimension [2][M/2][M/2] that takes two sp indices a,b and a spinquantumnumber S, and returns a tp index i: i = s2t[S][a][b]
      static int ***s2t;

      //!list of 6j symbols needed.
      static double **_6j;

      //!nr of particles
      static int N;

      //!dimension of sp hilbert space
      static int M;

};

#endif
