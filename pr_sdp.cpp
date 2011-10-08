/**
 * @mainpage 
 * This is an implementation of the dual only, potential reduction interior point method
 * for optimizing the 2.5DM using the I1, I2, Q1, Q2, G1 and G2 N-representability conditions.
 * Compiling can be done with the options I1 , I, I1Q2, IQ, IQG1, IQG2 and IQG with logical consequences for the program.
 * @author Brecht Verstichel, Ward Poelmans
 * @date 10-06-2011
 */

#include <iostream>
#include <fstream>
#include <cmath>

using std::cout;
using std::endl;
using std::ifstream;
using std::ofstream;

//includes all important headers and defines which conditions are
//going to be used:
#include "include.h"

/**
 * In the main the actual program is run.\n 
 * We start from the unity 2.5DM normed on the particle number and minimize the ojective function:\n\n
 * Tr (Gamma H) - t * ln(det P(Gamma)) \n\n
 * Once the minimum is found the parameter t is reduced and a new search is initiated,
 * this goes on until convergence is reached.\n
 * The potential is minimized using the Newton-Raphson method and the resulting linear system
 * is solved via the linear conjugate gradient method.
 */
int main(void){

   srand(time(NULL));

   cout.precision(10);

   const int M = 4;//nr of spatial orbitals
   const int N = 4;//nr of particles

   rxTPM::init(M,N);
   TPM::init(M,N);
   SPM::init(M,N);
   xSPM::init(M,N);
   PHM::init(M,N);
   xTPM::init(M,N);
   rxPHM::init(M,N);
   dDPM::init(M,N);
   dSPM::init(M,N);
   dTPM::init(M,N);
   ssdTPM::init(M,N);
   dPPHM::init(M,N);
   dPHHM::init(M,N);
   dDPV::init(M,N);
   dPPHV::init(M,N);
   dPHHV::init(M,N);
   EIG::init(M,N);
   SUP::init(M,N);
   Tools::init(M,N);

   dDPM ddpm;
   ddpm.fill_Random();
   
   ddpm.proj_W();

   dPHHM dphhm;
   dphhm.fill_Random();

   dDPM G1_down;
   G1_down.G1(dphhm);

   dPHHM G1_up;
   G1_up.G1(ddpm);

   cout << G1_up.ddot(dphhm) << "\t" << G1_down.ddot(ddpm) << endl;
/*
   //hamiltoniaan
   dDPM ham;
   ham.hubbard(1.0);

   dDPM W;
   W.unit();

   dDPM backup(W);

   double t = 1.0;
   double tolerance = 1.0e-5;

   //outer iteration: scaling of the potential barrier
   while(t > 1.0e-10){

      cout << t << "\t" << W.trace() << "\t" << W.ddot(ham) << "\t";

      int nr_cg_iter = 0;
      int nr_newton_iter = 0;

      double convergence = 1.0;

      //inner iteration: 
      //Newton's method for finding the minimum of the current potential
      while(convergence > tolerance){

         ++nr_newton_iter;

         SUP P;

         P.fill(W);

         P.invert();

         //eerst -gradient aanmaken:
         dDPM grad;

         grad.constr_grad(t,ham,P);

         //dit wordt de stap:
         dDPM delta;

         //los het hessiaan stelsel op:
         nr_cg_iter += delta.solve(t,P,grad);

         //line search
         double a = delta.line_search(t,P,ham);

         //W += a*delta;
         W.daxpy(a,delta);

         convergence = a*a*delta.ddot(delta);

      }

      cout << nr_newton_iter << "\t" << nr_cg_iter << endl;

      t /= 1.5;

      //what is the tolerance for the newton method?
      tolerance = 1.0e-5*t;

      if(tolerance < 1.0e-12)
         tolerance = 1.0e-12;

      //extrapolatie:
      dDPM extrapol(W);

      extrapol -= backup;

      //overzetten voor volgende stap
      backup = W;

      double b = extrapol.line_search(t,W,ham);

      W.daxpy(b,extrapol);

   }

   cout << endl;

   cout << "Final Energy:\t" << ham.ddot(W) << endl;
*/
   Tools::clear();
   dDPM::clear();
   xTPM::clear();
   PHM::clear();
   TPM::clear();
   rxTPM::clear();
   rxPHM::clear();

   return 0;

}
