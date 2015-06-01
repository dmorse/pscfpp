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
   
      void setId(int id);
      void addBlock(const Block& block);

      int id() const;
      int size() const;
      const GArray< Pair<int> >& inSolverIds() const;
      const GArray< Pair<int> >& outSolverIds() const;
   
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
   const GArray< Pair<int> >& Vertex::inSolverIds() const
   {  return inSolverIds_; }

   inline 
   const GArray< Pair<int> >& Vertex::outSolverIds() const
   {  return outSolverIds_; }

} 
#endif 
