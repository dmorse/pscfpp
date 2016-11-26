#ifndef PSCF_HOMOGENEOUS_MIXTURE_H
#define PSCF_HOMOGENEOUS_MIXTURE_H

/*
* PSCF - Molecule Self-Consistent Field Theory
*
* Copyright 2016, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include <util/param/ParamComposite.h>    // base class

#include <pscf/chem/Monomer.h>            // Member template argument
#include <pscf/homogeneous/Molecule.h>    // Member template argument
#include <util/containers/DArray.h>       // Member template

namespace Pscf {
namespace Homogeneous {

   using namespace Util;

   /**
   * A spatially homogeneous mixture.
   *
   * \ingroup Pscf_Homogeneous_Module
   */
   class Mixture : public ParamComposite
   {
   public:

      /**
      * Constructor.
      */
      Mixture();

      /**
      * Destructor.
      */
      ~Mixture();

      /**
      * Read parameters from file and initialize.
      *
      * \param in input parameter file
      */
      virtual void readParameters(std::istream& in);

      /**
      * Set the number of molecular species and allocate memory.
      */
      void setNMolecule(int nMonomer);
      
      /**
      * Set the number of monomer types.
      */
      void setNMonomer(int nMonomer);

      /// \name Accessors (by non-const reference)
      //@{
 
      /**
      * Get a molecule object.
      *
      * \param id integer molecule species index (0 <= id < nMolecule)
      */
      Molecule& molecule(int id);

      /**
      * Get number of molecule species.
      */
      int nMolecule() const;

      /**
      * Get number of monomer types.
      */
      int nMonomer() const;

      //@}

      /**
      * Validate all data structures.
      *
      * Throw an exception if an error is found.
      */
      void validate() const;

   private:

      /**
      * Array of molecule species solver objects.
      *
      * Array capacity = nMolecule.
      */
      DArray<Molecule> molecules_;

      /**
      * Array of molecular volume fractions.
      */
      DArray<double> phi_;

      /**
      * Array of molecular chemical potentials. 
      */
      DArray<double> mu_;

      /**
      * Array of monomer volume fractions.
      */
      DArray<double> f_;

      /**
      * Number of molecule species.
      */
      int nMolecule_;

      /**
      * Number of monomer types (maximum monomer id + 1).
      */
      int nMonomer_;

   };

   // Inline member functions

   inline Molecule& Mixture::molecule(int id)
   {  return molecules_[id]; }

   inline int Mixture::nMolecule() const
   {  return nMolecule_; }

   inline int Mixture::nMonomer() const
   {  return nMonomer_; }

}
}
#endif
