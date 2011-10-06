#ifndef rxPHM_H
#define rxPHM_H

#include <iostream>
#include <fstream>
#include <vector>

using std::ostream;
using std::vector;

#include "Matrix.h"

class dPHHM;

/**
 * @author Brecht Verstichel
 * @date 29-09-2011\n\n
 * This class rxPHM is a class written for the blocks of the dPHHM matrices, it is called rxPHM because it is a reduced expansion of a PHM matrix.
 * expanded for the same reason as xTPM and reduced because there is an extra parameter l which the second sp index cannot be equal to if intermediate spin S_bl = 1.
 */
class rxPHM : public BlockMatrix {

   /**
    * Output stream operator overloaded
    * @param output The stream to which you are writing (e.g. cout)
    * @param rxphm_p the rxPHM you want to print
    */
   friend ostream &operator<<(ostream &output,const rxPHM &rxphm_p);

   public:
      
      //constructor
      rxPHM(int);

      //copy constructor
      rxPHM(const rxPHM &);

      //destructor
      virtual ~rxPHM();

      using BlockMatrix::operator=;

      using BlockMatrix::operator();

      //geef N terug
      int gN() const;

      //geef M terug
      int gM() const;

      //geef l terug
      int gl() const;

      static void init(int,int);

      static void clear();

      static int gph2s(int,int,int,int);

      static int gs2ph(int,int,int,int,int);

      static void print_basis();

   private:

      //!static list that takes in a paramater l, a phh-spinindex S and a ph index i and returns two sp indices a and b and intermediate spin S_bl
      static vector< vector<int> > **ph2s;

      //!static list that takes in a parameter l, a phh-spinindex S, intermediate spin S_bl and two sp indices a,b and returns a ph index i
      static int *****s2ph;

      //!parameter l, sp index that is blocked out
      int l;

      //!nr of particles
      static int N;

      //!dimension of sp hilbert space
      static int M;

};

#endif
