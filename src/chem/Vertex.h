#ifndef PSCF_CHEM_VERTEX_H
#define PSCF_CHEM_VERTEX_H

/*
* PFTS - Polymer Field Theory Simulator
*
* Copyright 2013, David Morse (morse012@.umn.edu)
* Distributed under the terms of the GNU General Public License.
*/

namespace Pfts{ 
namespace Chem{

   /**
   * A junction or chain end in a block copolymer.
   */
   class Vertex
   {
   public:
   
      void addBlock(Block& block, int end);
      int degree();
   
   private:
   
      std::vector<int> blockIds_;
   
   };

} 
} 
#endif 
