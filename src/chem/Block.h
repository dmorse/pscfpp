#ifndef PSCF_CHEM_BLOCK_H
#define PSCF_CHEM_BLOCK_H

/*
* PFTS - Polymer Field Theory Simulator
*
* Copyright 2013, David Morse (morse012@.umn.edu)
* Distributed under the terms of the GNU General Public License.
*/

#include <util/containers/Pair.h>

namespace Pfts{ 
namespace Chem{

   using namespace Util;

   /**
   * A linear homopolymer block within a block copolymer.
   */
   class Block
   {
   
      Block();
   
      void setId(unsigned int id);
      void setVertexIds(int VertexAId, int VertexBId);
      void setMonomerId(int monomerId);
      void setLength(double length);
   
      int id() const;
      int monomerId() const;
      const Pair<int>& vertexIds() const;
      int vertexId() const;
      double length() const;
    
   private:
   
      int id_;
      int monomerId_;
      Pair<int> vertexIds_;
      double length_;

   };

} 
} 
#endif 
