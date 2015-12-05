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
      int nVertex() const;

      /**
      * Number of solvers (twice nBlock).
      */
      int nSolver() const;  //

      /**
      * Volume per molecule, in units of reference volume.
      */
      double volume() const;

      /**
      * Get a specified Block.
      * 
      * \param id block index
      */
      const Block& block(int id) const;

      /**
      * Get a specified Vertex.
      * 
      * \param id vertex index
      */
      const Vertex& vertex(int id) const;

      /**
      * Solver identifier, indexed by order of computation.
      *
      * An array of propagator ids ordered in the order in which 
      * they should be computed, so that the intitial condition 
      * for each link is provided by the solution of links that 
      * have already been solved.
      */
      const Pair<int>& solverId(int i) const;

   protected:

      virtual void makePlan();

   private:

      /// Array of Block objects in this polymer.
      DArray<Block> blocks_;

      /// Array of Vertex objects in this polymer.
      DArray<Vertex> vertices_;

      /// Solver ids, indexed in order in which they should be invoked.
      DArray< Pair<int> > solverIds_;

      /// Number of blocks in this polymer
      int nBlock_;

      /// Number of vertices (ends or junctions) in this polymer
      int nVertex_;

      /// Number of solvers (two per block).
      int nSolver_;

      double volume_;

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
   * Number of solvers.
   */
   inline int PolymerDescriptor::nSolver() const
   {  return nSolver_; }

   /*
   * Volume per molecule, in units of reference volume.
   */
   inline double PolymerDescriptor::volume() const
   {  return volume_; }

   /*
   * Get a specified Block.
   */
   inline const Block& PolymerDescriptor::block(int id) const
   {  return blocks_[id]; }

   /*
   * Get a specified Vertex.
   */
   inline const Vertex& PolymerDescriptor::vertex(int id) const
   {  return vertices_[id]; }

   /*
   * Get a solver id, indexed in order of computation.
   */
   inline 
   const Pair<int>& PolymerDescriptor::solverId(int id) const
   {  return solverIds_[id]; }

}
#endif 
