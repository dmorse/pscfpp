#ifndef PSCF_VERTEX_H
#define PSCF_VERTEX_H

/*
* PSCF - Polymer Self-Consistent Field Theory
*
* Copyright 2016 - 2022, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include <util/containers/GArray.h>
#include <util/containers/Pair.h>

namespace Pscf
{ 

   class Edge;
   using namespace Util;

   /**
   * A junction or chain end in a block polymer.
   *
   * A vertex represents a vertex within the graph associated with a 
   * linear or acyclic branched polymer. Each vertex is either a free
   * end of a block or a junction at which two or more blocks (or edges)
   * are connected. 
   *
   * \ingroup Pscf_Chem_Module
   */
   class Vertex
   {
   public:

      /**
      * Constructor.
      */
      Vertex();

      /**
      * Destructor.
      */
      ~Vertex();
  
      /**
      * Set the integer identifier of this vertex.
      * 
      * \param id identifier
      */ 
      void setId(int id);

      /**
      * Add block to the list of attached blocks.
      *
      * Preconditions: The id for this vertex must have been set, vertex
      * ids must have been set for the block, and the id of this vertex
      * must match one of the ids for the two vertices attached to the
      * block.
      * 
      * \param block attached Edge object
      */ 
      void addEdge(Edge const & block);

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
      * and the second is a direction id which is 0 if this 
      * vertex is vertex 1 of the block, and 1 if this vertex
      * is vertex 0.
      *
      * \param i index of incoming propagator
      * \return Pair<int> containing block index, direction index
      */
      Pair<int> const & inPropagatorId(int i) const;

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
      Pair<int> const & outPropagatorId(int i) const;
   
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
   Pair<int> const & Vertex::inPropagatorId(int i) const
   {  return inPropagatorIds_[i]; }

   inline 
   Pair<int> const & Vertex::outPropagatorId(int i) const
   {  return outPropagatorIds_[i]; }

} 
#endif 
