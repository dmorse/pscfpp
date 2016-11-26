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
      monomers_(),
      molecules_(),
      nMonomer_(0), 
      nMolecule_(0)
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
      // Monomers
      read<int>(in, "nMonomer", nMonomer_);
      monomers_.allocate(nMonomer_);
      readDArray< Monomer >(in, "monomers", monomers_, nMonomer_);

      // Molecules 
      read<int>(in, "nMolecule", nMolecule_);
      molecules_.allocate(nMolecule_);
      for (int i = 0; i < nMolecule_; ++i) {
         readParamComposite(in, molecules_[i]);
      }

   }

} // namespace Homogeneous
} // namespace Pscf
