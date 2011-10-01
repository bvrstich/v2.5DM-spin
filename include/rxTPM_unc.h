#ifndef rxTPM_unc_H
#define rxTPM_unc_H

#include <iostream>
#include <fstream>
#include <vector>

using std::ostream;
using std::vector;

#include "Matrix.h"

/**
 * @author Brecht Verstichel
 * @date 27-05-2011\n\n
 * This class rxTPM_unc is a class written for two particle matrices with extra parameter l, which the sp indices cannot be equal to.
 */
class rxTPM_unc : public Matrix {

   /**
    * Output stream operator overloaded
    * @param output The stream to which you are writing (e.g. cout)
    * @param rtpm_p the rxTPM_unc you want to print
    */
   friend ostream &operator<<(ostream &output,const rxTPM_unc &rtpm_p);

   public:
      
      //constructor
      rxTPM_unc(int);

      //copy constructor
      rxTPM_unc(const rxTPM_unc &);

      //destructor
      virtual ~rxTPM_unc();

      using Matrix::operator=;

      using Matrix::operator();

      //easy to access the numbers, in sp mode
      double operator()(int s_ab,int a,int b,int s_cd,int c,int d) const;

      static int get_inco(int,int &,int &,int &);

      //geef N terug
      int gN() const;

      //geef M terug
      int gM() const;

      //geef l terug
      int gl() const;

      static void init(int,int);

      static void clear();

      static void print_basis();

      static int gt2s(int,int,int);

      static int gs2t(int,int,int,int);

   private:

      //!static list that takes in a parameter l and a tp index i and returns two sp indices a,b and a spinindex s_l
      static vector< vector<int> > *t2s;

      //!static list that takes in a parameter l, a spinindex s_l and two sp indices a,b and returns a tp index i
      static int ****s2t;

      //!parameter l, spatial sp index that is blocked out.
      int l;

      //!nr of particles
      static int N;

      //!nr of spatial sp orbs
      static int M;

};

#endif
