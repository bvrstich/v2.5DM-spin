#ifndef dPHHV_H
#define dPHHV_H

#include <iostream>
#include <fstream>

using std::ostream;

#include "Vector.h"

/**
 * @author Brecht Verstichel
 * @date 30-05-2011\n\n
 * This class dPHHV is a class written the eigenvalues of the dPHHM object.
 */
class dPHHV {

   /**
    * Output stream operator overloaded
    * @param output The stream to which you are writing (e.g. cout)
    * @param ddpv_p the dPHHV you want to print
    */
   friend ostream &operator<<(ostream &output,const dPHHV &ddpv_p);

   public:

      //standard constructor
      dPHHV(dPHHM &);
      
      //copy constructor
      dPHHV(const dPHHV &);

      //destructor
      virtual ~dPHHV();

      //geef N terug
      int gN() const;

      //geef M terug
      int gM() const;

      BlockVector<rxPHM> &operator[](int);

      const BlockVector<rxPHM> &operator[](int) const;

      dPHHV &operator=(const dPHHV &);

      dPHHV &operator=(double );

      dPHHV &operator+=(const dPHHV &);

      dPHHV &operator-=(const dPHHV &);

      dPHHV &daxpy(double alpha,const dPHHV &);

      double sum() const;

      double log_product() const;

      double ddot(const dPHHV &) const;

      void dscal(double alpha);

      double min() const;
      
      double max() const;

      double lsfunc(double) const;

      static void init(int,int);

   private:
      
      //!double pointer to vector<rxPHM> that contain all the eigenvalues of a dPHHM object
      BlockVector<rxPHM> **dphhv;

      //!nr of particles
      static int N;

      //!dimension of sp hilbert space
      static int M;

};

#endif
