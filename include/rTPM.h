#ifndef rTPM_H
#define rTPM_H

#include <iostream>
#include <fstream>
#include <vector>

using std::ostream;
using std::vector;

#include "Matrix.h"

class dDPM;

/**
 * @author Brecht Verstichel
 * @date 04-08-2011\n\n
 * This class rTPM is a class written for spinsymmetrical two particle matrices with extra parameter l, which the spatial sp indices cannot be equal to.
 */
class rTPM : public BlockMatrix {

   /**
    * Output stream operator overloaded
    * @param output The stream to which you are writing (e.g. cout)
    * @param rtpm_p the rTPM you want to print
    */
   friend ostream &operator<<(ostream &output,const rTPM &rtpm_p);

   public:
      
      //constructor
      rTPM(int);

      //copy constructor
      rTPM(const rTPM &);

      //destructor
      virtual ~rTPM();

      using BlockMatrix::operator=;

      using BlockMatrix::operator();

      //easy to access the numbers, in sp mode
      double operator()(int S,int S_ab,int a,int b,int S_cd,int c,int d) const;

      //geef N terug
      int gN() const;

      //geef M terug
      int gM() const;

      //geef l terug
      int gl() const;

      static void init(int,int);

      static void clear();

      static int gt2s(int,int,int,int);

      static int gs2t(int,int,int,int,int);

   private:

      //!static list that takes in a paramater l, a dp-spinindex S and a tp index i and returns two sp indices a and b and intermediate spin S_ab
      static vector< vector<int> > **t2s;

      //!static list that takes in a parameter l, a dp-spinindex S, intermediate spin S_ab and two sp indices a,b and returns a tp index i
      static int *****s2t;

      //!parameter l, sp index that is blocked out.
      int l;

      //!nr of particles
      static int N;

      //!dimension of sp hilbert space
      static int M;

};

#endif
