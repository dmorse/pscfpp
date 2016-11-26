/*
* PSCF - Molecule Self-Consistent Field Theory
*
* Copyright 2016, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include "Mixture.h"

namespace Pscf {
namespace Homogeneous {

   using namespace Util;

   /*
   * Constructor.
   */
   Mixture::Mixture()
    : ParamComposite(),
      molecules_(),
      nMolecule_(0),
      nMonomer_(0)
   {  setClassName("Mixture"); }

   /*
   * Destructor.
   */
   Mixture::~Mixture()
   {}

   /*
   * Read all parameters and initialize.
   */
   void Mixture::readParameters(std::istream& in)
   {
      // Precondition
      UTIL_ASSERT(molecules_.capacity() == 0);

      read<int>(in, "nMonomer", nMonomer_);
      f_.allocate(nMonomer_);

      read<int>(in, "nMolecule", nMolecule_);
      molecules_.allocate(nMolecule_);
      phi_.allocate(nMolecule_);
      mu_.allocate(nMolecule_);
      for (int i = 0; i < nMolecule_; ++i) {
         readParamComposite(in, molecules_[i]);
      }

      validate();
   }

   void Mixture::setNMolecule(int nMolecule)
   {
      UTIL_ASSERT(molecules_.capacity() == 0);
      nMolecule_ = nMolecule;
      molecules_.allocate(nMolecule_);
      phi_.allocate(nMolecule_);
      mu_.allocate(nMolecule_);
   }

   void Mixture::setNMonomer(int nMonomer)
   {
      UTIL_ASSERT(nMonomer_ == 0);
      nMonomer_ = nMonomer;
      f_.allocate(nMonomer_);
   }

   /*
   * Check validity and complete initialization.
   */
   void Mixture::validate() const
   {
      UTIL_ASSERT(nMolecule_ > 0);
      UTIL_ASSERT(nMonomer_ > 0);
      for (int i = 0; i < nMolecule_; ++i) {
         Molecule const & mol = molecules_[i];
         UTIL_ASSERT(mol.nClump() > 0);
         for (int j = 0; j < mol.nClump(); ++j) {
            UTIL_ASSERT(mol.clump(j).monomerId() < nMonomer_);
         }
      }
   }

} // namespace Homogeneous
} // namespace Pscf
