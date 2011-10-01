#ifndef xTPM_unc_H
#define xTPM_unc_H

#include <iostream>
#include <fstream>
#include <vector>

using std::ostream;
using std::vector;

#include "Matrix.h"

/**
 * @author Brecht Verstichel
 * @date 27-05-2011\n\n
 * This class xTPM_unc is a class written for two particle matrices with extra parameter l, which the sp indices cannot be equal to.
 */
class xTPM_unc : public Matrix {

   /**
    * Output stream operator overloaded
    * @param output The stream to which you are writing (e.g. cout)
    * @param rtpm_p the xTPM_unc you want to print
    */
   friend ostream &operator<<(ostream &output,const xTPM_unc &rtpm_p);

   public:
      
      //constructor
      xTPM_unc();

      //copy constructor
      xTPM_unc(const xTPM_unc &);

      //destructor
      virtual ~xTPM_unc();

      using Matrix::operator=;

      using Matrix::operator();

      //easy to access the numbers, in sp mode
      double operator()(int s_ab,int a,int b,int s_cd,int c,int d) const;

      static int get_inco(int,int,int);

      //geef N terug
      int gN() const;

      //geef M terug
      int gM() const;

      static void init(int,int);

      static void clear();

      static void print_basis();

      static int gt2s(int,int);

      static int gs2t(int,int,int);

   private:

      //!static list that takes in a parameter a tp index i and returns two sp indices a,b and a spinindex s_l
      static vector< vector<int> > t2s;

      //!static list that takes in a spinindex s_l and two sp indices a,b and returns a tp index i
      static int ***s2t;

      //!nr of particles
      static int N;

      //!nr of spatial sp orbs
      static int M;

};

#endif
