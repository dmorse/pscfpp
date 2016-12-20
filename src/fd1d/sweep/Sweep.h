#ifndef FD1D_SWEEP_H
#define FD1D_SWEEP_H

/*
* PSCF - Polymer Self-Consistent Field Theory
*
* Copyright 2016, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include <util/param/ParamComposite.h>    // base class
#include <util/global.h>                  

namespace Pscf {
namespace Fd1d
{

   class Mixture;
   class Domain;
   class System;
   class Iterator;

   using namespace Util;

   /**
   * Solve a sequence of problems along a line in parameter space.
   *
   * \ingroup Pscf_Fd1d_Module
   */
   class Sweep : public ParamComposite
   {

   public:

      /**
      * Default Constructor.
      * 
      * Objects instantiated with this constructor must also call
      * the setSystem() function.
      */
      Sweep();

      /**
      * Constructor.
      * 
      * \param system parent System object.
      */
      Sweep(System& system);

      /**
      * Destructor.
      */
      ~Sweep();

      /**
      * Create association with the parent System.
      *
      * Use iff instantiated with default constructor.
      * 
      * \param system parent System object.
      */
      void setSystem(System& system);

      /**
      * Read ns and baseFileName parameters.
      *
      * \param in input stream
      */
      virtual void readParameters(std::istream& in);

      /**
      * Setup operation at beginning sweep.
      */
      virtual void setup(){};

      /**
      * Set system parameters.
      *
      * \param s path length coordinate, in range [0,1]
      */
      virtual void setState(double s) = 0;

      /**
      * Output information after obtaining a converged solution.
      *
      * \param stateFileName base name of output files
      * \param s value of path length parameter s
      */
      virtual void outputSolution(std::string const & stateFileName, double s);

      /**
      * Output data to a running summary.
      *
      * \param outFile  output file, open for writing
      * \param i  integer index
      * \param s  value of path length parameter s
      */
      virtual void outputSummary(std::ostream& outFile, int i, double s);

      /**
      * Iterate to solution.
      */
      virtual void solve();

   protected:

      /// Get parent System by reference.
      System & system();

      /// Get parent System by const reference.
      System const & system() const;
     
      /// Get associated Mixture by reference. 
      Mixture& mixture();

      /// Get associated Mixture by const reference. 
      Mixture const & mixture() const;
      
      /// Get associated Domain by reference. 
      Domain& domain();

      /// Get associated Domain by const reference. 
      Domain const & domain() const;
      
      /// Get associated Iterator by reference. 
      Iterator& iterator();

      /// Get associated Iterator by const reference. 
      Iterator const & iterator() const;
     
   protected:

      /// Number of steps. 
      int ns_;

      /// Mode for comparison to homogeneous system (none -> -1)
      int homogeneousMode_;

      /// Base name for output files
      std::string baseFileName_;

   private:

      /// Pointer to parent System object.
      System* systemPtr_;

      /// Pointer to associated Mixture object.
      Mixture* mixturePtr_;

      /// Pointer to associated Domain object.
      Domain* domainPtr_;

      /// Pointer to associated Iterator object.
      Iterator* iteratorPtr_;

   };

   inline
   System& Sweep::system()
   {
      UTIL_ASSERT(systemPtr_);
      return *systemPtr_;
   }

   inline
   System const & Sweep::system() const
   {
      UTIL_ASSERT(systemPtr_);
      return *systemPtr_;
   }

   inline
   Mixture& Sweep::mixture()
   {
      UTIL_ASSERT(mixturePtr_);
      return *mixturePtr_;
   }

   inline
   Mixture const & Sweep::mixture() const
   {
      UTIL_ASSERT(mixturePtr_);
      return *mixturePtr_;
   }

   inline
   Domain& Sweep::domain()
   {
      UTIL_ASSERT(domainPtr_);
      return *domainPtr_;
   }

   inline
   Domain const & Sweep::domain() const
   {
      UTIL_ASSERT(domainPtr_);
      return *domainPtr_;
   }

   inline
   Iterator& Sweep::iterator()
   {
      UTIL_ASSERT(iteratorPtr_);
      return *iteratorPtr_;
   }

   inline
   Iterator const & Sweep::iterator() const
   {
      UTIL_ASSERT(iteratorPtr_);
      return *iteratorPtr_;
   }

} // namespace Fd1d
} // namespace Pscf
#endif
