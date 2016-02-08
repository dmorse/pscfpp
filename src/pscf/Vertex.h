#ifndef PSCF_VERTEX_H
#define PSCF_VERTEX_H

/*
* PFTS - Polymer Field Theory Simulator
*
* Copyright 2013, David Morse (morse012@.umn.edu)
* Distributed under the terms of the GNU General Public License.
*/

#include <util/containers/GArray.h>
#include <util/containers/Pair.h>

namespace Pscf
{ 

   class BlockDescriptor;
   using namespace Util;

   /**
   * A junction or chain end in a block polymer.
   *
   * \ingroup Pscf_Base_Module
   */
   class Vertex
   {
   public:

      Vertex();
      ~Vertex();
  
      /**
      * Set the integer identifier of this vertex.
      * 
      * \param id identifier
      */ 
      void setId(int id);

      /**
      * Add to the list of attached blocks.
      * 
      * \param block attached BlockDescriptor object
      */ 
      void addBlock(const BlockDescriptor& block);

      /**
      * Get the id of this vertex.
      */
      int id() const;

      /**
      * Get the number of attached blocks.
      */
      int size() const;

      /**
      * Get the block and direction of an incoming propagator.
      *
      * The first element of the integer pair is the block id,
      * and the second is a direction id which is 1 if this 
      * vertex is vertex 1 of the block, and 1 if this vertex
      * is vertex 0.
      *
      * \param i index of incoming propagator
      * \return Pair<int> containing block index, direction index
      */
      const Pair<int>& inPropagatorId(int i) const;

      /**
      * Get the block and direction of an outgoing propagator
      *
      * The first element of the integer pair is the block id,
      * and the second is a direction id which is 0 if this 
      * vertex is vertex 0 of the block, and 1 if this vertex
      * is vertex 1.
      *
      * \param i index of incoming propagator
      * \return Pair<int> containing block index, direction index
      */
      const Pair<int>& outPropagatorId(int i) const;
   
   private:
   
      GArray< Pair<int> > inPropagatorIds_;
      GArray< Pair<int> > outPropagatorIds_;
      int id_;
   
   };

   inline int Vertex::id() const
   {  return id_; }

   inline int Vertex::size() const
   {  return outPropagatorIds_.size(); }

   inline 
   const Pair<int>& Vertex::inPropagatorId(int i) const
   {  return inPropagatorIds_[i]; }

   inline 
   const Pair<int>& Vertex::outPropagatorId(int i) const
   {  return outPropagatorIds_[i]; }

} 
#endif 
