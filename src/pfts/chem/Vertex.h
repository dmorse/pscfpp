#ifndef PFTS_VERTEX_H
#define PFTS_VERTEX_H

/*
* PFTS - Polymer Field Theory Simulator
*
* Copyright 2013, David Morse (morse012@.umn.edu)
* Distributed under the terms of the GNU General Public License.
*/

#include <util/containers/GArray.h>
#include <util/containers/Pair.h>

namespace Pfts{ 

   class Block;
   using namespace Util;

   /**
   * A junction or chain end in a block polymer.
   */
   class Vertex
   {
   public:

      Vertex();
      ~Vertex();
   
      void setId(int id);
      void addBlock(const Block& block);

      int id() const;
      int size() const;
      const Pair<int>& inSolverId(int i) const;
      const Pair<int>& outSolverId(int i) const;
   
   private:
   
      GArray< Pair<int> > inSolverIds_;
      GArray< Pair<int> > outSolverIds_;
      int id_;
   
   };

   inline int Vertex::id() const
   {  return id_; }

   inline int Vertex::size() const
   {  return outSolverIds_.size(); }

   inline 
   const Pair<int>& Vertex::inSolverId(int i) const
   {  return inSolverIds_[i]; }

   inline 
   const Pair<int>& Vertex::outSolverId(int i) const
   {  return outSolverIds_[i]; }

} 
#endif 
