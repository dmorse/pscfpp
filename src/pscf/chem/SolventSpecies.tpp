/*
* PSCF - Polymer Self-Consistent Field
*
* Copyright 2015 - 2025, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include "SolventSpecies.h"

namespace Pscf {

   /*
   * Constructor
   */
   template <typename WT>
   SolventSpecies<WT>::SolventSpecies()
    : Species<WT>(),
      monomerId_(-1),
      size_(0.0)
   {  ParamComposite::setClassName("SolventSpecies"); }

   /*
   * Read contents of parameter file block
   */
   template <typename WT>
   void SolventSpecies<WT>::readParameters(std::istream& in)
   {
      ParamComposite::read<int>(in, "monomerId", monomerId_);
      ParamComposite::read<double>(in, "size", size_);

      // Read phi or mu (but not both) and set ensemble accordingly
      Species<WT>::readParameters(in);
   }

   /*
   * Set the monomer type id for this solvent species.
   */ 
   template <typename WT>
   void SolventSpecies<WT>::setMonomerId(int monomerId)
   {  monomerId_ = monomerId; }
  
   /*
   * Set the size parameter for this solvent.
   */ 
   template <typename WT>
   void SolventSpecies<WT>::setSize(double size)
   {  size_ = size; }

}
