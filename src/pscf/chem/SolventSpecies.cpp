/*
* PSCF - Polymer Self-Consistent Field Theory
*
* Copyright 2016 - 2022, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include "SolventSpecies.h"

namespace Pscf {

   /*
   * Constructor
   */
   SolventSpecies::SolventSpecies()
    : Species(),
      monomerId_(-1),
      size_(0.0)
   {  setClassName("SolventSpecies"); }

   /*
   * Destructor
   */
   SolventSpecies::~SolventSpecies()
   {}

   /*
   * Read contents of parameter file block
   */
   void SolventSpecies::readParameters(std::istream& in)
   {
      read<int>(in, "monomerId", monomerId_);
      read<double>(in, "size", size_);

      // Read phi or mu (but not both) and set ensemble accordingly
      Species::readParameters(in);
   }

   /*
   * Set the monomer type id for this solvent species.
   */ 
   void SolventSpecies::setMonomerId(int monomerId)
   {  monomerId_ = monomerId; }
  
   /*
   * Set the size parameter for this solvent.
   */ 
   void SolventSpecies::setSize(double size)
   {  size_ = size; }

}
