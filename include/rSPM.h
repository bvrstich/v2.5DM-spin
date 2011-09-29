#ifndef rSPM_H
#define rSPM_H

#include <iostream>
#include <fstream>

using std::ostream;

#include "BlockMatrix.h"

/**
 * @author Brecht Verstichel
 * @date 19-04-2010\n\n
 * This class rSPM is a class written for the blocks of the dTPM object, which is a TPM object with the second index diagonal
 * It consists of two blocks with equal dimensions and inherits from BlockMatrix.
 */
class rSPM : public BlockMatrix {

   public:
      
      rSPM();

      //copy constructor
      rSPM(const rSPM &);

      //destructor
      virtual ~rSPM();

      void constr_lists();

      using BlockMatrix::operator=;

      using BlockMatrix::operator();

      int gN() const;

      int gM() const;

      static void init(int,int);

   private:

      //!nr of particles
      static int N;

      //!dimension of sp hilbert space
      static int M;

};

#endif
