#ifndef dSPM_H
#define dSPM_H

#include <iostream>
#include <fstream>

using std::ostream;

class dDPM;

/**
 * @author Brecht Verstichel
 * @date 23-06-2011\n\n
 * This class dSPM is an array of doubles, it will contain the traces of the different blocks of the dDPM and consorts.
 */
class dSPM {

   /**
    * Output stream operator overloaded
    * @param output The stream to which you are writing (e.g. cout)
    * @param dspm_p the dSPM you want to print
    */
   friend ostream &operator<<(ostream &output,const dSPM &dspm_p);

   public:

      //constructor
      dSPM();

      //copy constructor
      dSPM(const dSPM &);

      //destructor
      virtual ~dSPM();

      double &operator[](int);

      double operator[](int) const;

      int gN() const;

      int gM() const;

      void trace(double,const dDPM &);

      void trace(double,const dPPHM &);

      static void init(int,int);

   private:

      //!nr of particles
      static int N;

      //!dimension of spatial sp space (nr of spatial orbs)
      static int M;

      //!array of doubles
      double *dspm;

};

#endif
