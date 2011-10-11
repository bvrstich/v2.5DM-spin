#ifndef EIG_H
#define EIG_H

#include <iostream>
#include <fstream>

using std::ostream;

#include "SUP.h"

/**
 * @author Brecht Verstichel
 * @date 30-05-2011\n\n
 * This class, EIG is a "block"-vector over the carrierspace's of the active condtions. It contains room
 * to store the eigenvalues and special member function that work with these eigenvalues.
 * This class should only be used when a SUP matrix has been diagonalized, some functions could give strange results when the EIG object is filled
 * with random numbers.\n\n
 */
class EIG{

   /**
    * Output stream operator overloaded,
    * @param output The stream to which you are writing (e.g. cout)
    * @param eig_p the EIG you want to print
    */
   friend ostream &operator<<(ostream &output,const EIG &eig_p);

   public:

      //constructor met initialisatie op 
      EIG(SUP &);
      
      //copy constructor
      EIG(const EIG &);

      //destructor
      ~EIG();

      int gN() const;

      int gM() const;

      int gdim() const;

      //overload equality operator
      EIG &operator=(const EIG &);

      dDPV &gv_I1();

      const dDPV &gv_I1() const;

#ifdef __Q2_CON
      dDPV &gv_Q2();

      const dDPV &gv_Q2() const;
#endif

#ifdef __I2_CON
      dPPHV &gv_I2();

      const dPPHV &gv_I2() const;
#endif

#ifdef __Q1_CON
      dPPHV &gv_Q1();

      const dPPHV &gv_Q1() const;
#endif

#ifdef __G1_CON
      dPHHV &gv_G1();

      const dPHHV &gv_G1() const;
#endif

#ifdef __G2_CON
      dPHHV &gv_G2();

      const dPHHV &gv_G2() const;
#endif

      double min() const;

      double max() const;

      double lsfunc(double) const;

      static void init(int,int);

   private:

      //!pointer to a dDPV object that contains the eigenvalues of the I1 condition of a SUP matrix 
      dDPV *v_I1;

#ifdef __Q2_CON
      //!pointer to a dDPV object that contains the eigenvalues of the Q2 condition of a SUP matrix 
      dDPV *v_Q2;
#endif

#ifdef __I2_CON
      //!pointer to a dPPHV object that contains the eigenvalues of the I2 condition of a SUP matrix 
      dPPHV *v_I2;
#endif

#ifdef __Q1_CON
      //!pointer to a dPPHV object that contains the eigenvalues of the Q1 condition of a SUP matrix 
      dPPHV *v_Q1;
#endif

#ifdef __G1_CON
      //!pointer to a dPPHV object that contains the eigenvalues of the G1 condition of a SUP matrix 
      dPHHV *v_G1;
#endif

#ifdef __G2_CON
      //!pointer to a dPPHV object that contains the eigenvalues of the G2 condition of a SUP matrix 
      dPHHV *v_G2;
#endif

      //!number of particles
      static int N;

      //!dimension of sp space
      static int M;

      //!total dimension of the EIG object
      static int dim;

};

#endif
