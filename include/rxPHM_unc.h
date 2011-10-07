#ifndef rxPHM_unc_H
#define rxPHM_unc_H

#include <iostream>
#include <fstream>
#include <vector>

using std::ostream;
using std::vector;

#include "Matrix.h"

/**
 * @author Brecht Verstichel
 * @date 27-05-2011\n\n
 * This class rxPHM_unc is a class written for two particle matrices with extra parameter l, which the sp indices cannot be equal to.
 */
class rxPHM_unc : public Matrix {

   /**
    * Output stream operator overloaded
    * @param output The stream to which you are writing (e.g. cout)
    * @param rphm_p the rxPHM_unc you want to print
    */
   friend ostream &operator<<(ostream &output,const rxPHM_unc &rxphm_p);

   public:
      
      //constructor
      rxPHM_unc(int);

      //copy constructor
      rxPHM_unc(const rxPHM_unc &);

      //destructor
      virtual ~rxPHM_unc();

      using Matrix::operator=;

      using Matrix::operator();

      //easy to access the numbers, in sp mode
      double operator()(int s_ab,int a,int b,int s_cd,int c,int d) const;

      //geef N terug
      int gN() const;

      //geef M terug
      int gM() const;

      //geef l terug
      int gl() const;

      static void init(int,int);

      static void clear();

      static void print_basis();

      static int gph2s(int,int,int);

      static int gs2ph(int,int,int,int);

   private:

      //!static list that takes in a parameter l and a ph index i and returns two sp indices a,b and a spinindex s_l
      static vector< vector<int> > *ph2s;

      //!static list that takes in a parameter l, a spinindex s_l and two sp indices a,b and returns a ph index i
      static int ****s2ph;

      //!parameter l, spatial sp index that is blocked out.
      int l;

      //!nr of particles
      static int N;

      //!nr of spatial sp orbs
      static int M;

};

#endif
