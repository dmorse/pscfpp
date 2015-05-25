#ifndef PSCF_CHEM_VERTEX_H
#define PSCF_CHEM_VERTEX_H

/*
* PFTS - Polymer Field Theory Simulator
*
* Copyright 2013, David Morse (morse012@.umn.edu)
* Distributed under the terms of the GNU General Public License.
*/

#include <vector>

namespace Pfts{ 
namespace Chem{

   class Block;

   /**
   * A junction or chain end in a block polymer.
   */
   class Vertex
   {
   public:

      Vertex();
   
      void addBlock(int blockId, int end);
      int degree();
   
   private:
   
      std::vector<int> blockIds_;
   
   };

} 
} 
#endif 
