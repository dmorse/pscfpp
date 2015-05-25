
/*
* PSCF++ - Polymer Self-Consistent Field Theory Simulator
*
* Copyright 2013, David Morse (morse012@.umn.edu)
* Distributed under the terms of the GNU General Public License.
*/
   
   /**
   * Interface for a species.
   */
   class Species
   {
   
   public:
   
      enum Ensemble {UNKNOWN, CLOSED, OPEN};
   
      Species(System& system);
      setEnsemble(Ensemble& ensemble);
      setPhi();
      setMu();
   
      /**
      * Solve modified diffusion equation.
      */
      virtual void compute() = 0;
   
      /**
      * Get overall occupied volume fraction.
      */
      double phi();
   
      /**
      * Get chemical potential.
      */
      double mu();
   
      /**
      * Overall partition function.
      */
      double q();
   
   protected:
   
      double phi_;
      double mu_;
      double q_
      bool isComputed_;
   
   };
   
