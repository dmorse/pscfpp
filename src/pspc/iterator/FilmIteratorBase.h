#ifndef PSPC_FILM_ITERATOR_BASE_H
#define PSPC_FILM_ITERATOR_BASE_H

/*
* PSCF - Polymer Self-Consistent Field Theory
*
* Copyright 2016 - 2021, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include "pspc/iterator/Iterator.h"        // base class
#include "util/containers/FSArray.h"       // container
#include <string>
#include <iostream>

namespace Pscf {
namespace Pspc
{

   template <int D>
   class System;
   
   using namespace Util;

   /**
   * Descriptor for a FilmIterator object. This base class template defines
   * all traits of a FilmIterator that do not depend on D, the dimension of
   * the system. The subclasses of FilmIteratorBase are the partial
   * specializations of FilmIterator for 1D, 2D, and 3D.
   * 
   * If the user chooses a FilmIterator as their iterator, then the 
   * system will contain two parallel hard surfaces ("walls"), confining
   * the polymers/solvents to a "thin film" region of the unit cell.
   * This only affects the iterator, not the rest of SCFT, so we isolate
   * the imposition of the thin film constraint to this subclass of 
   * iterator. This is essentially a wrapper for any other type of iterator
   * (e.g., AmIterator), that adds the additional functionality required
   * to impose the thin film constraint properly.
   * 
   * FilmIterator is generalized to be compatible with any iterator within
   * it, as long as the iterator can impose 1) a mask that confines the 
   * polymers/solvents to a certain region of space, and 2) an external 
   * field. A FilmIterator object in the param file is created by appending
   * "Film" to the end of the name of the iterator that is stored inside of
   * FilmIterator (e.g., "AmIteratorFilm{" will create a FilmIterator with
   * an AmIterator object inside of it). 
   *
   * \ingroup Pspc_Iterator_Module
   */
   template <int D, typename IteratorType>
   class FilmIteratorBase : public Iterator<D>
   {

   public:

      /**
      * Constructor.
      */
      FilmIteratorBase(System<D>& system);

      /**
      * Destructor.
      */
      ~FilmIteratorBase();

      /**
      * Read and initialize.
      *
      * \param in  input parameter stream
      */
      void readParameters(std::istream& in);

      /**
      * Allocate required memory, perform necessary checks to ensure user 
      * input is compatible with a film constraint, and create the mask / 
      * external fields that will be used to represent the walls during 
      * iteration.
      */
      void setup();

      /**
      * Iterate to a solution
      */
      int solve();

      /**
      * Return const reference to the real iterator within this FilmIterator
      */
      IteratorType const & iterator() const;

      /** 
      * Return whether unit cell is flexible.
      */ 
      bool isFlexible() const;

      /**
      * Determine the indices of each flexible lattice parameter, based on
      * normalVecId and unitCell definitions in param file. Assumes that
      * isFlexible == true, and gives a warning if none of the parameters
      * are actually able to be varied given the thin film constraints. 
      * Stores resulting array in flexibleParams_ member of this object,
      * as well as the flexibleParams_ member of the iterator within this
      * object.
      * 
      * This function varies depending on D, the dimensionality of the
      * system. Therefore, it is implemented in the partial class 
      * specializations in FilmIterator, rather than in this base class.
      */
      virtual void setFlexibleParams() = 0;

      /**
      * Check that user-defined lattice basis vectors (stored in the
      * Domain<D> object associated with this FilmIterator class)
      * are compatible with the thin film constraint
      * 
      * This function varies depending on D, the dimensionality of the
      * system. Therefore, it is implemented in the partial class 
      * specializations in FilmIterator, rather than in this base class.
      */
      virtual void checkLatticeVectors() const = 0;

      /**
      * Generates the field representation of the walls, based on the values
      * of wallThickness and interfaceThickness that were input by the user.
      * Then, passes this wall field into the iterator_ object to be used
      * as a mask during iteration, and also passes the corresponding
      * external potential fields into iterator_ if isAthermal() == false.
      */
      void generateWallFields();

      /**
      * Checks whether the lattice parameters have been updated since the
      * last call of generateWallFields(), and if the parameters have
      * changed then calls generateWallFields() again to update them. 
      * 
      * Also updates the external fields if the wall/polymer chi parameters
      * have been updated since external fields were last generated.
      */
     void updateWallFields();

      /**
      * Iterates over all polymers and solvents in the Mixture, and 
      * multiplies their volume fractions (phi values) by (1-phi_w), where 
      * phi_w is the volume fraction occupied by the walls. This correction 
      * ensures that the system will still be able to converge to a 
      * solution that satisfies the incompressibility constraint. 
      */
      void adjustPhiVals();

      /**
      * Check that user-defined space group is compatible with the
      * thin film constraint
      */
      void checkSpaceGroup() const;

      /**
      * Check whether or not the walls are chemically identical using
      * the chi array.
      */
      bool isSymmetric() const;

      /**
      * Check whether or not the walls are athermal (this is only true
      * if all of the values in the chi array are zero)
      */
      bool isAthermal() const;

      /**
      * Set the value of chi between species s and wall w
      * 
      * \param w  wall index, 0 or 1
      * \param s  species index, 0 <= id < nVertex
      * \param chi  value of chi(s,w)
      */
      void setChi(int s, int w, double chi);

      /**
      * Write the concentration field of the walls to an output stream, 
      * in basis format.
      * 
      * \param out  output stream
      */
      void writeWallField(std::ostream &out) const;

      /**
      * Write the concentration field of the walls to a file specified 
      * by 'filename', in basis format.
      * 
      * \param filename  name of file to which data will be written
      */
      void writeWallField(std::string filename) const;

      /**
      * Write the concentration field of the walls to an output stream, in 
      * r-grid format.
      * 
      * \param out  output stream
      */
      void writeWallFieldRGrid(std::ostream &out) const;

      /**
      * Write the concentration field of the walls to a file specified by 
      * 'filename', in r-grid format.
      * 
      * \param filename  name of file to which data will be written
      */
      void writeWallFieldRGrid(std::string filename) const;

      /**
      * Returns the field that describes the mask imposed upon the unit 
      * cell, which represents the wall. Field is in symmetry-adapted 
      * basis format.
      */
      DArray<double> const & maskBasis() const;

      /**
      * Returns the array of fields that describes the external potential 
      * field felt by each monomer species, which is a result of chemical 
      * interactions between the wall and the monomer species. Fields are 
      * in symmetry-adapted basis format.
      */
      DArray< DArray<double> > const & externalFields() const;

      /**
      * Returns the field that describes the external potential field felt 
      * by monomer species monomerId, which is a result of chemical 
      * interactions between the wall and the monomer species. Fields are 
      * in symmetry-adapted basis format.
      *
      * \param monomerId integer monomer type index
      */
      DArray<double> const & externalField(int monomerId) const;

      /**
      * Get const value of normalVecId
      */
      int normalVecId() const;

      /**
      * Get const value of interfaceThickness
      */
      double interfaceThickness() const;

      /**
      * Get const value of wallThickness
      */
      double wallThickness() const;

      /**
      * Get const chi matrix by reference
      */
      DMatrix<double> const & chi() const;

      /**
      * Get the const chi parameter between wall w and species s by 
      * reference
      * 
      * \param s  species index, 0 <= id < nVertex
      * \param w  wall index, 0 or 1
      */
      double chi(int s, int w) const;

      using Iterator<D>::maskPhi;

   protected:

      /**
      * Return reference to the real iterator within this FilmIterator
      */
      IteratorType& iterator();

      /**
      * Generate external fields only, and pass them into iterator object 
      * for use during calculation. This is called by generateWallFields().
      */
      void generateExternalFields();

      using Iterator<D>::system;
      using Iterator<D>::setClassName;
      using Iterator<D>::isFlexible_;
      using Iterator<D>::hasMask_;
      using Iterator<D>::hasExternalFields_;
      using ParamComposite::read;
      using ParamComposite::readOptional;
      using ParamComposite::readOptionalDMatrix;
      using ParamComposite::setParent;
      using ParamComposite::addComponent;

   private:

      /// The actual iterator that does all the work
      IteratorType iterator_;

      /// Lattice parameters associated with the current maskBasis
      FSArray<double, 6> parameters_;

      /// Wall chi array associated with the current externalFields
      DMatrix<double> chiFields_;

      /// Lattice basis vector that is normal to the walls
      int normalVecId_;

      /// Interface thickness
      double t_;

      /// Wall thickness
      double T_;

      /// Chi matrix
      DMatrix<double> chi_;

      /// Name of file to which to write wall field in basis format
      std::string wallBasisOutFile_;

      /// Name of file to which to write wall field in r-grid format
      std::string wallRGridOutFile_;

   };

   // Inline member functions

   // Return reference to iterator within this FilmIterator
   template <int D, typename IteratorType>
   inline IteratorType& FilmIteratorBase<D, IteratorType>::iterator()
   {  return iterator_; }

   // Return const reference to iterator within this FilmIterator
   template <int D, typename IteratorType>
   inline 
   IteratorType const & FilmIteratorBase<D, IteratorType>::iterator() const
   {  return iterator_; }

   // Return whether unit cell is flexible.
   template <int D, typename IteratorType>
   inline bool FilmIteratorBase<D, IteratorType>::isFlexible() const
   {  return iterator_.isFlexible(); }

   // Set value of chi between species s and wall w
   template <int D, typename IteratorType>
   inline 
   void FilmIteratorBase<D, IteratorType>::setChi(int s, int w, double chi)
   {  chi_(s,w) = chi; }

   // Get const mask (wall) field by reference
   template <int D, typename IteratorType>
   inline 
   DArray<double> const & FilmIteratorBase<D, IteratorType>::maskBasis() 
   const
   {  return iterator_.maskBasis(); }

   // Get const array of external fields by reference
   template <int D, typename IteratorType>
   inline 
   DArray< DArray<double> > const & 
   FilmIteratorBase<D, IteratorType>::externalFields() const
   {  return iterator_.externalFields(); }

   // Get const external field for species monomerId by reference
   template <int D, typename IteratorType>
   inline 
   DArray<double> const & 
   FilmIteratorBase<D, IteratorType>::externalField(int monomerId)
   const
   {  return iterator_.externalField(monomerId); }

   // Get value of normalVecId
   template <int D, typename IteratorType>
   inline int FilmIteratorBase<D, IteratorType>::normalVecId() const
   {  return normalVecId_; }

   // Get value of interfaceThickness
   template <int D, typename IteratorType>
   inline double FilmIteratorBase<D, IteratorType>::interfaceThickness() 
   const
   {  return t_; }

   // Get value of wallThickness
   template <int D, typename IteratorType>
   inline double FilmIteratorBase<D, IteratorType>::wallThickness() const
   {  return T_; }

   // Get chi matrix by const reference
   template <int D, typename IteratorType>
   inline DMatrix<double> const & FilmIteratorBase<D, IteratorType>::chi() 
   const
   {  return chi_; }

   // Get the chi parameter between wall w and species s by reference
   template <int D, typename IteratorType>
   inline double FilmIteratorBase<D, IteratorType>::chi(int s, int w) 
   const
   {  return chi_(s, w); }

} // namespace Pspc
} // namespace Pscf

#include "FilmIteratorBase.tpp"
#endif