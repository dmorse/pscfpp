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

      /// \name Accessors (by non-const reference)
      //@{
 
      /**
      * Get a Monomer type descriptor.
      *
      * \param id integer monomer type index (0 <= id < nMonomer)
      */
      Monomer& monomer(int id);

      /**
      * Get a molecule object.
      *
      * \param id integer molecule species index (0 <= id < nMolecule)
      */
      Molecule& molecule(int id);

      //@}
      /// \name Accessors (by value)
      //@{
 
      /**
      * Get number of monomer types.
      */
      int nMonomer() const;

      /**
      * Get number of molecule species.
      */
      int nMolecule() const;

      //@}

   private:

      /**
      * Array of monomer type descriptors.
      */
      DArray<Monomer> monomers_;

      /**
      * Array of molecule species solver objects.
      *
      * Array capacity = nMolecule.
      */
      DArray<Molecule> molecules_;

      /**
      * Number of monomer types.
      */
      int nMonomer_; 

      /**
      * Number of molecule species.
      */
      int nMolecule_;

   };

   // Inline member functions

   inline int Mixture::nMonomer() const
   {  return nMonomer_; }

   inline int Mixture::nMolecule() const
   {  return nMolecule_; }

   inline Monomer& Mixture::monomer(int id)
   {  return monomers_[id]; }

   inline Molecule& Mixture::molecule(int id)
   {  return molecules_[id]; }

}
}
#endif
