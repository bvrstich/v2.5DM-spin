#ifndef dPPHV_H
#define dPPHV_H

#include <iostream>
#include <fstream>

using std::ostream;

#include "Vector.h"

/**
 * @author Brecht Verstichel
 * @date 30-05-2011\n\n
 * This class dPPHV is a class written the eigenvalues of the dPPHM object.
 */
class dPPHV {

   /**
    * Output stream operator overloaded
    * @param output The stream to which you are writing (e.g. cout)
    * @param dpphv_p the dPPHV you want to print
    */
   friend ostream &operator<<(ostream &output,const dPPHV &dpphv_p);

   public:

      //standard constructor
      dPPHV(dPPHM &);
      
      //copy constructor
      dPPHV(const dPPHV &);

      //destructor
      virtual ~dPPHV();

      //geef N terug
      int gN() const;

      //geef M terug
      int gM() const;

      BlockVector<xTPM> &operator[](int);

      const BlockVector<xTPM> &operator[](int) const;

      dPPHV &operator=(const dPPHV &);

      dPPHV &operator=(double );

      dPPHV &operator+=(const dPPHV &);

      dPPHV &operator-=(const dPPHV &);

      dPPHV &daxpy(double alpha,const dPPHV &);

      double sum() const;

      double log_product() const;

      double ddot(const dPPHV &) const;

      void dscal(double alpha);

      double min() const;
      
      double max() const;

      double lsfunc(double) const;

      static void init(int,int);

   private:
      
      //!double pointer to vector<xTPM> that contain all the eigenvalues of a dPPHM object
      BlockVector<xTPM> **dpphv;

      //!nr of particles
      static int N;

      //!dimension of sp hilbert space
      static int M;

};

#endif
