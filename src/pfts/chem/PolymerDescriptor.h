#ifndef PFTS_POLYMER_DESCRIPTOR_H
#define PFTS_POLYMER_DESCRIPTOR_H

/*
* PFTS - Polymer Field Theory Simulator
*
* Copyright 2013, David Morse (morse012@.umn.edu)
* Distributed under the terms of the GNU General Public License.
*/

#include <util/param/ParamComposite.h>   // base class
#include <util/containers/DArray.h>      // member template
#include <util/containers/Pair.h>        // member template
#include "Monomer.h"                     // member template argument
#include "Block.h"                       // member template argument
#include "Vertex.h"                      // member template argument

namespace Pfts{ 

   using namespace Util;

   /**
   * Descriptor for an acyclic block polymer.
   */
   class PolymerDescriptor : public ParamComposite
   {
   public:

      /**
      * Constructor.
      */
      PolymerDescriptor();

      /**
      * Read and initialize polymer structure.
      *
      * \param in input parameter stream
      */
      void readParameters(std::istream& in);

      /**
      * Number of blocks.
      */
      int nBlock() const;  //

      /**
      * Number of vertices (ends and/or junctions).
      */
      int  nVertex() const;

      /**
      * Volume per molecule, in units of reference volume.
      */
      double volume() const;

      /**
      * Get a specified Block.
      * 
      * \param id block index
      */
      const Block& blocks(int id) const;

      /**
      * Get a specified Vertex.
      * 
      * \param id vertex index
      */
      const Vertex& vertex(int id) const;

      /**
      * Computation plan.
      *
      * The links array is ordered in the order of the computation
      * plan, i.e., in an order in which the links must be solved
      * so that the intitial condition for each link is provided
      * by the solution of links that have already been solved.
      */
      const DArray< Pair<int> >& plan() const;

   private:

      DArray<Block> blocks_;
      DArray<Vertex> vertices_;
      DArray< Pair<int> > plan_;
      double volume_;
      int nBlock_;
      int nVertex_;

      void makePlan();

   };
  
   /*
   * Number of blocks.
   */
   inline int PolymerDescriptor::nBlock() const
   {  return nBlock_; }

   /*
   * Number of vertices (ends and/or junctions)
   */
   inline int PolymerDescriptor::nVertex() const
   {  return nVertex_; }

   /*
   * Volume per molecule, in units of reference volume.
   */
   inline double PolymerDescriptor::volume() const
   {  return volume_; }

   /*
   * Get a specified Block.
   */
   inline const Block& PolymerDescriptor::blocks(int id) const
   {  return blocks_[id]; }

   /*
   * Get a specified Vertex.
   */
   inline const Vertex& PolymerDescriptor::vertex(int id) const
   {  return vertices_[id]; }

} 
#endif 
