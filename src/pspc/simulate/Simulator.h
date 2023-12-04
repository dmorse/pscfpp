#ifndef PSPC_SIMULATOR_H
#define PSPC_SIMULATOR_H

/*
* PSCF - Polymer Self-Consistent Field Theory
*
* Copyright 2016 - 2022, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include <util/param/ParamComposite.h>     // base class

#include <prdc/cpu/RField.h>               // memmber (template arg)
#include <util/random/Random.h>            // member
#include <util/containers/DArray.h>        // member (template)
#include <util/containers/DMatrix.h>       // member (template)

namespace Pscf {
namespace Pspc {

   template <int D> class System;

   using namespace Util;
   using namespace Prdc::Cpu;

   /**
   * Field theoretic simulator (base class).
   *
   * \ingroup Pspc_Simulate_Module
   */
   template <int D>
   class Simulator : public ParamComposite
   {

   public:

      /**
      * Constructor.
      *
      * \param system parent System
      */
      Simulator(System<D>& system);

      /**
      * Destructor.
      */
      ~Simulator();

      /**
      * Allocate required memory.
      */
      void allocate();

      /**
      * Read parameters for a simulation.
      *
      * \param in input parameter stream
      */
      virtual void readParameters(std::istream &in);

      /// \name Primary Actions: Simulation and Analysis
      ///@{

      /**
      * Perform a field theoretic Monte-Carlo simulation.
      *
      * Perform a field theoretic simulation of nSteps using the
      * partial saddle-point approximation.
      *
      * The default implemention is a do-nothing placeholder that throws
      * an error if called, and must be re-implemented by subclasses.
      *
      * \param nStep  number of simulation steps
      */
      virtual void simulate(int nStep);

      /**
      * Read and analyze a trajectory file.
      *
      * This function uses an instance of the TrajectoryReader class
      * specified by the "classname" argument to read a trajectory 
      * file with the specified named "filename". The function must
      * open file, performs the analysis, and closes the file before
      * returning.
      *
      * The default implemention is a do-nothing placeholder that throws
      * an error if called, and must be re-implemented by subclasses.
      *
      * \param min  first frame number
      * \param max  last frame number
      * \param classname  name of TrajectoryReader class
      * \param filename  name of trajectory file
      */
      virtual void analyze(int min, int max,
                           std::string classname,
                           std::string filename);

      /**
      * Output timing results
      */
      virtual void outputTimers(std::ostream& out);

      /**
      * Clear timers
      */
      virtual void clearTimers();

      /**
      * Return the current simulation step index.
      */
      long iStep();

      ///@}
      /// \name Projected Chi Matrix
      ///@{

      /**
      * Perform eigenvalue analysis of projected chi matrix.
      */
      void analyzeChi();

      /**
      * Get an array of the eigenvalues of the projected chi matrix.
      *
      * The projected chi matrix is given by the matrix product P*chi*P,
      * where P is the symmetric projection matrix that projects onto 
      * the subspace perpendicular to the vector e = (1,1,...,1). The 
      * projected chi matrix is singular, and has a zero eigenvalue with
      * associated eigenvector e. By convention, this zero eigenvalue 
      * and its eigenvector e are listed last, with index nMonomer - 1.
      */
      DArray<double> const & chiEvals() const
      {  return chiEvals_; }

      /**
      * Get a single eigenvalue of the projected chi matrix.
      *
      * \param i index of eigenvalue (0, ... , nMonomer - 1)
      */
      double chiEval(int i ) const
      {  return chiEvals_[i]; }

      /**
      * Get the matrix of all eigenvectors of the projected chi matrix.
      *
      * This function returns the nMonomer x nMonomer matrix of the
      * eigenvectors of the projected chi Matrix, in which each row
      * is an eigenvector. The first (row) index of this matrix thus
      * identifies an eigenvector, while the second (column) index 
      * identifiers the monomer type of each component of such a vector. 
      *
      * Each eigenvector is normalized such that the sum of the squares of 
      * all elements is equal to nMonomer, the number of monomer types.  
      * The sign of each vector is chosen so as to make the first (0) 
      * component non-negative.  The last eigenvector is always the null 
      * vector e = (1,1,...,1).
      * 
      * For nMonomer = 2, the two eigenvectors are (1,-1) and (1,1).
      */
      DMatrix<double> const & chiEvecs() const
      {  return chiEvecs_; }

      /**
      * Get components of S in a basis of eigenvectors of chiP.
      *
      * The value of component \f$ S_{a} \f$ may be expressed using 
      * Einstein summation convention as
      * \f[
      *     S_{a} \equiv \frac{1}{M^2}} v_{ai}\chi_{ij}e_{j}
      * \f]
      * for any \f$ a = 0, \ldots, M - 1 \f$, where M = nMonomer (the
      * number of monomer types), \f$ e_{j} =1 \f$ for any j, and 
      * \f$ v_{ai} \f$ is component associated with monomer type i of 
      * eigenvector a of the projected chi matrix,
      */
      DArray<double> const & sc() const;

      /**
      * Get a single component of the S vector.
      *
      * This function retrieves on component of the vector defined in
      * the documentation for function sc().
      * 
      * \param i  eigenvector index (0, ..., nMonomer - 1)
      */
      double sc(int i) const;

      ///@}
      /// \name Field Theoretic Hamiltonian 
      ///@{

      /**
      * Compute the Hamiltonian used in Monte-Carlo simulations.
      */
      void computeHamiltonian();

      /**
      * Get the Hamiltonian used in Monte-Carlo simulations.
      */
      double hamiltonian() const;

      /**
      * Get ideal gas contribution (lnQ) to MC Hamiltonian.
      */
      double idealHamiltonian() const;

      /**
      * Get the quadratci field contribution (HW) to MC Hamiltonian.
      */
      double fieldHamiltonian() const;

      /**
      * Has the MC Hamiltonian been computed for current w and c fields?
      */
      bool hasHamiltonian() const;

      /**
      * Clear field eigen-components and hamiltonian components.
      */
      void clearData();

      ///@}
      /// \name Chemical Potential Field (W Field) Components
      ///@{

      /**
      * Compute eigenvector components of the current w fields.
      *
      * Compute and store the components of the values of the w fields
      * on nodes of a real-space grid (r-grid) in a basis of the
      * eigenvectors of the projected chi matrix. 
      */
      void computeWc();

      /**
      * Get one eigenvector component of the current w fields.
      *
      * Each component is a point-wise projection of the w fields onto a
      * corresponding eigenvector of the projected chi matrix. The last
      * index, i = nMonomer - 1, corresponds to the Lagrange multiplier
      * pressure component, with associated eigenvector [1, 1, ..., 1].
      */
      DArray< RField<D> > const & wc() const;

      /**
      * Get one eigenvector component of the current w fields.
      *
      * See documentation of function wc().
      *
      * \param i eigenvector / eigenvalue index
      */
      RField<D> const & wc(int i) const;

      /**
      * Are eigen-components of current w fields valid ?
      */
      bool hasWc() const;

      ///@}
      /// \name Monomer Concentration Field (C-Field) Components
      ///@{

      /**
      * Compute eigenvector components of the current c fields.
      *
      * Compute and store the components of the values of the c fields
      * on nodes of a real-space grid (r-grid) in a basis of the
      * eigenvectors of the projected chi matrix. 
      */
      void computeCc();

      /**
      * Get all eigenvector component of the current w fields.
      *
      * Each component is a point-wise projection of the c fields onto a
      * corresponding eigenvector of the projected chi matrix. 
      */
      DArray< RField<D> > const & cc() const;

      /**
      * Get one eigenvector component of the current c fields.
      *
      * See documentation of function cc().
      *
      * \param i eigenvector / eigenvalue index
      */
      RField<D> const & cc(int i) const;

      /**
      * Are eigen-components of current c fields valid ?
      */
      bool hasCc() const;

      ///@}
      /// \name Functional Derivatives of H[W]
      ///@{

      /**
      * Compute functional derivatives of the Hamiltonian.
      *
      * Compute and store the functional derivatives of the Hamiltonian
      * with respect to eigenvector components of the w fields.
      */
      void computeDc();

      /**
      * Get all of the current d fields.
      */
      DArray< RField<D> > const & dc() const;

      /**
      * Get one eigenvector component of the current d fields.
      *
      * \param i eigenvector / eigenvalue index
      */
      RField<D> const & dc(int i) const;

      /**
      * Are the current d fields valid ?
      */
      bool hasDc() const;

      ///@}
      /// \name Miscellaneous
      ///@{

      /**
      * Get parent system by reference.
      */
      System<D>& system();

      /**
      * Get random number generator by reference.
      */
      Random& random();

      ///@}

   protected:

      using Util::ParamComposite::setClassName;

      /**
      * Random number generator
      */
      Random random_;

      /**
      * Eigenvector components of w fields on a real space grid.
      *
      * Each field component corresponds to a point-wise projection of w
      * onto an eigenvector of the projected chi matrix.
      */
      DArray< RField<D> > wc_;

      /**
      * Eigenvector components of c fields on a real space grid.
      *
      * Each field component corresponds to a point-wise projection of c
      * onto an eigenvector of the projected chi matrix.
      */
      DArray< RField<D> > cc_;

      /**
      * Components of d fields on a real space grid.
      *
      * Each field component is the functional derivative of H[W]
      * with respect to one eigenvector w-field component.
      */
      DArray< RField<D> > dc_;

      /**
      * Field theoretic Hamiltonian H[W] (extensive value).
      */
      double hamiltonian_;

      /**
      * Ideal gas contribution (lnQ) to Hamiltonian H[W]
      */
      double idealHamiltonian_;

      /**
      * Field contribution (H_W) to Hamiltonian
      */
      double fieldHamiltonian_;

      /**
      * Simulation step counter.
      */
      long iStep_;

      /**
      * Has the Hamiltonian been computed for the current w and c fields?
      */
      bool hasHamiltonian_;

      /**
      * Have eigen-components of the current w fields been computed ?
      */
      bool hasWc_;

      /**
      * Have eigen-components of the current c fields been computed ?
      */
      bool hasCc_;

      /**
      * Have functional derivatives of H[W] been computed ?
      */
      bool hasDc_;

   private:

      /**
      * Projected chi matrix
      *
      * Projected matrix chiP_ = P*chi*P, where P = I - e e^{T} / M is
      * a projection matrix that projects onto the subspace orthogonal
      * to the vector e = [1, ... , 1]^{T}, where M = nMonomer.
      */
      DMatrix<double> chiP_;

      /**
      * Eigenvectors of the projected chi matrix.
      *
      * Each row (identified by first index) is an eigenvector. 
      * The last eigenvector, with index nMonomer - 1, is always the
      * vector e = [1, 1, ...., 1].
      */
      DMatrix<double> chiEvecs_;

      /**
      * Eigenvalues of the projected chi matrix.
      *
      * The last eigenvalue, with index nMonomer - 1, is always zero.
      */
      DArray<double>  chiEvals_;

      /**
      * Components of vector s = chi*e in a basis of eigenvectors.
      *
      * Component sc_[a] is given by sc_[a] = v_{a} chi e / M^2, 
      * where e = [1 1 ... 1]^{T}, v_{a} is a row vector representation
      * of eigenvector a of the projected chi matrix, given by row a of 
      * chiEvecs_, and M = nMonomer.
      */
      DArray<double>  sc_;

      /**
      * Pointer to the parent system.
      */
      System<D>* systemPtr_;

   };

   // Inline functions

   // Get the random number generator.
   template <int D>
   inline Random& Simulator<D>::random()
   {  return random_; }

   // Get the parent System.
   template <int D>
   inline System<D>& Simulator<D>::system()
   {  return *systemPtr_; }

   // Get the precomputed Hamiltonian
   template <int D>
   inline double Simulator<D>::hamiltonian() const
   {
      UTIL_CHECK(hasHamiltonian_);
      return hamiltonian_;
   }

   // Get the ideal gas component of the precomputed Hamiltonian
   template <int D>
   inline double Simulator<D>::idealHamiltonian() const
   {
      UTIL_CHECK(hasHamiltonian_);
      return idealHamiltonian_;
   }

   // Get the W field component of the precomputed Hamiltonian.
   template <int D>
   inline double Simulator<D>::fieldHamiltonian() const
   {
      UTIL_CHECK(hasHamiltonian_);
      return fieldHamiltonian_;
   }

   // Has the Hamiltonian been computed for the current w fields ?
   template <int D>
   inline bool Simulator<D>::hasHamiltonian() const
   {  return hasHamiltonian_; }

   // Return all eigencomponents of the w fields.
   template <int D>
   inline DArray< RField<D> > const & Simulator<D>::wc() const
   {  return wc_; }

   // Return a single eigencomponent of the w fields.
   template <int D>
   inline RField<D> const & Simulator<D>::wc(int i) const
   {  return wc_[i]; }

   // Have eigen-components of current w fields been computed?
   template <int D>
   inline bool Simulator<D>::hasWc() const
   {  return hasWc_; }

   // Return all eigencomponents of the current c fields.
   template <int D>
   inline DArray< RField<D> > const & Simulator<D>::cc() const
   {  return cc_; }

   // Return a single eigencomponent of the current c fields.
   template <int D>
   inline RField<D> const & Simulator<D>::cc(int i) const
   {  return cc_[i]; }

   // Have eigen-components of current c fields been computed?
   template <int D>
   inline bool Simulator<D>::hasCc() const
   {  return hasCc_; }

   // Return all eigencomponents of the current c fields.
   template <int D>
   inline DArray< RField<D> > const & Simulator<D>::dc() const
   {  return dc_; }

   // Return a single eigencomponent of the current c fields.
   template <int D>
   inline RField<D> const & Simulator<D>::dc(int i) const
   {  return dc_[i]; }

   // Have eigen-components of current c fields been computed?
   template <int D>
   inline bool Simulator<D>::hasDc() const
   {  return hasDc_; }

   // Clear all data (eigen-components of w field and Hamiltonian)
   template <int D>
   inline void Simulator<D>::clearData()
   {
      hasHamiltonian_ = false;
      hasWc_ = false;
      hasCc_ = false;
      hasDc_ = false;
   }

   template <int D>
   inline long Simulator<D>::iStep()
   {  return iStep_; }

   #ifndef PSPC_SIMULATOR_TPP
   // Suppress implicit instantiation
   extern template class Simulator<1>;
   extern template class Simulator<2>;
   extern template class Simulator<3>;
   #endif

}
}
#endif
