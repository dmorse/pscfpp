#ifndef PSCF_BLOCK_DESCRIPTOR_H
#define PSCF_BLOCK_DESCRIPTOR_H

/*
* PSCF - Polymer Self-Consistent Field Theory
*
* Copyright 2016 - 2019, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include "PolymerType.h"
#include <util/containers/Pair.h>

#include <iostream>

namespace Pscf
{ 

   using namespace Util;

   /**
   * Description of a linear homopolymer block within a block polymer.
   *
   * This class defines the blockId, monomerId, length and vertexIds of 
   * a block within a block polymer. It serves as a base class for the
   * BlockTmpl class template, which is a template for classes that solve 
   * the modified diffusion equation for the two propagators associated 
   * with a block. VertexIds should be set for all blocks in a block 
   * polymer before the associated Vertex objects are initialized.
   *
   * \ref pscf_Block_page "Parameter File Format"
   * \ingroup Pscf_Chem_Module
   */
   class BlockDescriptor
   {
   public:

      /**
      * Constructor.
      */ 
      BlockDescriptor();
 
      /**
      * Destructor.
      */ 
      virtual ~BlockDescriptor();
 
      /**
      * Serialize to/from archive.
      *
      * \param ar input or output Archive
      * \param versionId archive format version index
      */ 
      template <class Archive>
      void serialize(Archive& ar, unsigned int versionId);

      /// \name Setters
      //@{
    
      /**
      * Set the id for this block.
      *
      * \param id integer index for this block
      */ 
      void setId(int id);
  
      /**
      * Set indices of associated vertices.
      *
      * \param vertexAId integer id of vertex A
      * \param vertexBId integer id of vertex B
      */ 
      void setVertexIds(int vertexAId, int vertexBId);
  
      /**
      * Set the monomer id.
      *
      * \param monomerId integer id of monomer type (>=0)
      */ 
      void setMonomerId(int monomerId);
  
      /**
      * Set the length of this block.
      *
      * The ``length" is steric volume / reference volume.
      *
      * \param length block length (number of monomers).
      */ 
      virtual void setLength(double length);
  
      /**
      * Set the polymer type.
      *
      * By convention, if the polymer type of a block with block index id is
      * PolymerType::Linear, then vertexId(0) = id and vertexId(1) = id + 1.
      * The PolymerType enumeration value for the block is used by the
      * inserter and extractor operators to define a shorter string
      * representation for blocks in linear polymers, for which the string
      * representation does not include values for vertex ids. Vertex id
      * values for blocks in a linear poiymer must be set explicitly by 
      * calling the setVertexIds function with consecutive values, as done 
      * in the PolymerTmpl::readParameters function.
      *
      * \param type  type of polymer (branched or linear)
      */ 
      void setPolymerType(PolymerType::Enum type);
  
      //@}
      /// \name Accessors (getters)
      //@{
    
      /**
      * Get the id of this block.
      */ 
      int id() const;
  
      /**
      * Get the monomer type id.
      */ 
      int monomerId() const;
  
      /**
      * Get the pair of associated vertex ids.
      */ 
      const Pair<int>& vertexIds() const;
  
      /**
      * Get id of an associated vertex.
      *
      * \param i index of vertex (0 or 1)
      */ 
      int vertexId(int i) const;
  
      /**
      * Get the length (number of monomers) in this block.
      */
      double length() const;

      /**
      * Get the type of the parent polymer (branched or linear).
      */
      PolymerType::Enum polymerType() const;

      //@}

   private:
  
      /// Identifier for this block, unique within the polymer.
      int id_;

      /// Identifier for the associated monomer type.
      int monomerId_;

      /// Length of this block = volume / monomer reference volume. 
      double length_;

      /// Indexes of associated vertices
      Pair<int> vertexIds_;

      /// Type of polymer (branched or linear)
      PolymerType::Enum polymerType_;

      friend 
      std::istream& operator >> (std::istream& in, BlockDescriptor &block);

      friend 
      std::ostream& operator << (std::ostream& out, 
                                 const BlockDescriptor &block);

   };

   /**
   * istream extractor for a BlockDescriptor.
   *
   * \param in  input stream
   * \param block  BlockDescriptor to be read from stream
   * \return modified input stream
   */
   std::istream& operator >> (std::istream& in, BlockDescriptor &block);

   /**
   * ostream inserter for a BlockDescriptor.
   *
   * \param out  output stream
   * \param block  BlockDescriptor to be written to stream
   * \return modified output stream
   */
   std::ostream& operator << (std::ostream& out, const BlockDescriptor &block);

   // Inline member functions

   /*
   * Get the id of this block.
   */ 
   inline int BlockDescriptor::id() const
   {  return id_; }

   /*
   * Get the monomer type id.
   */ 
   inline int BlockDescriptor::monomerId() const
   {  return monomerId_; }

   /*
   * Get the pair of associated vertex ids.
   */ 
   inline const Pair<int>& BlockDescriptor::vertexIds() const
   {  return vertexIds_; }

   /*
   * Get id of an associated vertex.
   */ 
   inline int BlockDescriptor::vertexId(int i) const
   {  return vertexIds_[i]; }

   /*
   * Get the length (number of monomers) in this block.
   */
   inline double BlockDescriptor::length() const
   {  return length_; }
    
   /*
   * Get the polymer type (branched or linear).
   */
   inline PolymerType::Enum BlockDescriptor::polymerType() const
   {  return polymerType_; }
    
   /*
   * Serialize to/from an archive.
   */
   template <class Archive>
   void BlockDescriptor::serialize(Archive& ar, unsigned int)
   {
      ar & id_;
      ar & monomerId_;
      ar & vertexIds_;
      ar & length_;
   }
    
} 
#endif 
