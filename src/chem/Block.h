#ifndef PSCF_CHEM_BLOCK_H
#define PSCF_CHEM_BLOCK_H

/*
* PSCF++ - Polymer Self-Consistent Field Theory Simulator
*
* Copyright 2013, David Morse (morse012@.umn.edu)
* Distributed under the terms of the GNU General Public License.
*/

namespace Pscf{ 
namespace Core{

   /**
   * A linear homopolymer block within a block copolymer.
   */
   class Block
   {
   
      Block();
   
      setId(unsigned int id);
      setVertexIds(int VertexAId, int VertexBId);
      setMonomerId(int monomerId);
      setLength(double length);
   
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
