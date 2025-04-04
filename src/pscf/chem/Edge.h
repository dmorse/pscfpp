#ifndef PSCF_EDGE_H
#define PSCF_EDGE_H

/*
* PSCF - Polymer Self-Consistent Field Theory
*
* Copyright 2016 - 2022, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include "PolymerType.h"
#include "PolymerModel.h"
#include <util/containers/Pair.h>

#include <iostream>

namespace Pscf
{

   using namespace Util;

   /**
   * Descriptor for a block within an acyclic block polymer.
   *
   * This class defines the block id, monomerId, and vertexIds of a block 
   * within a block polymer. It also defines either the length of a block,
   * for a continuous thread model, or the number of beads in the block
   * for a discrete bead model. 
   *
   * This class serves as a base class for the BlockTmpl class template, 
   * which is a base class for a class named Block in each implementation 
   * level namespace. Each such implementation-levl Block class is thus 
   * a subclass of Edge.
   *
   * The class is designed to store either a value for length (the length
   * of the block) when PolymerModel::isThread(), or a value for nBead
   * (the number of beads in the block) when PolymerModel::isBead(). In
   * the case of a bead model, it also stores a pair of boolean flags 
   * to indicate whether the block owns either or both of the associated
   * vertices It is an error to try to set or get a value for nBead or
   * the vertex ownership flags when a thread model is in use, because
   * these variables are only meaningful for a bead model. Similarly,
   * it is an error to get or set a value for the block length when a
   * bead model is in use, because the length variable is only meaningful
   * in the context of a thread model. 
   *
   * Block objects associated with a polymer are normally stored in an
   * array named blocks. The id value of each block should be set to the
   * element index of the Block within that array.  VertexIds should be 
   * set for all blocks in a block polymer before the associated Vertex 
   * objects are initialized.
   *
   * \ref user_param_block_sec "Parameter File Format"
   * \ingroup Pscf_Chem_Module
   */
   class Edge
   {
   public:

      /**
      * Constructor.
      */
      Edge();

      /**
      * Destructor.
      */
      virtual ~Edge();

      /**
      * Serialize to/from archive.
      *
      * \param ar  input or output Archive
      * \param versionId  archive format version index
      */
      template <class Archive>
      void serialize(Archive& ar, unsigned int versionId);

      /// \name Setters
      //@{

      /**
      * Set the id for this block.
      *
      * \param id  integer index for this block
      */
      void setId(int id);

      /**
      * Set the monomer type id.
      *
      * \param monomerId  integer id of monomer type
      */
      void setMonomerId(int monomerId);

      /**
      * Set indices of associated vertices.
      *
      * \param vertexId0  integer id of vertex 0
      * \param vertexId1  integer id of vertex 1
      */
      void setVertexIds(int vertexId0, int vertexId1);

      /**
      * Set ownership of associated vertices (only valid for bead model).
      *
      * The concept of "ownership" of vertex beads is only meaningful in
      * the context of a bead model, in which we require that each vertex 
      * be owned by one of the associated blocks. The monomer type of the 
      * vertex is the same as the type of the block that owns it.
      *
      * Precondition: PolymerModel::isBead()
      *
      * \param own0  Does this block own vertex 0 ?
      * \param own1  Does this block own vertex 1 ?
      */
      void setVertexOwnership(bool own0, bool own1);

      /**
      * Set the number of beads in this block (only valid for bead model).
      *
      * Precondition: PolymerModel::isBead()
      *
      * \param nBead  number of beads (bead model)
      */
      virtual void setNBead(int nBead);

      /**
      * Set the length of this block (only valid for thread model).
      *
      * In the continuous thread model, the length of a block is given by
      * the ratio of block steric volume / monomer reference volume. 
      *
      * Precondition: PolymerModel::isThread()
      *
      * \param length  block length (thread model).
      */
      virtual void setLength(double length);

      /**
      * Set the polymer type (branched or linear).
      *
      * By convention, if the polymer type of a block with index id is
      * PolymerType::Linear, then vertexId(0) = id and vertexId(1) = id + 1.
      * In the case of a bead model, a standard convention also exists for
      * ownership of the two attached vertices.  The stream insertion and 
      * extraction operators for an Edge can thus use a shorter string
      * representation for linear polymers in which vertex ids and (when
      * relevant) owhership flags are omitted.
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
      * Does this block own an associated vertex (bead model).
      *
      * Precondition: PolymerModel::isBead()
      *
      * \param i index of vertex (0 or 1)
      */
      bool ownsVertex(int i) const;

      /**
      * Get the number of beads in this block, in the bead model.
      *
      * Precondition: PolymerModel::isBead()
      */
      double nBead() const;

      /**
      * Get the length of this block, in the thread model.
      *
      * Precondition: PolymerModel::isThread()
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

      /// Number of beads in block 
      /// Only valid for bead model, if PolymerModel::isBead()
      int nBead_;

      /// Length of this block = volume / monomer reference volume.
      /// Only valid for thread model, if PolymerModel::isThread()
      double length_;

      /// Indexes of associated vertices
      Pair<int> vertexIds_;

      /// Pair of bools to indicate ownership of associated vertices
      /// Only valid for bead model, if PolymerModel::isBead()
      Pair<bool> ownsVertex_;

      /// Type of polymer that contains this block (branched or linear)
      PolymerType::Enum polymerType_;

      friend
      std::istream& operator >> (std::istream& in, Edge &block);

      friend
      std::ostream& operator << (std::ostream& out,
                                 const Edge &block);

   };

   /**
   * Input stream extractor (>>) for a Edge.
   *
   * The polymerType must be set before an Edge can be read from a stream.
   * The block id must be set explicitly by calling setId, and are not 
   * read from an istream. Vertex id values for blocks in a linear polymer 
   * must be set explicitly by calling setVertexIds with consecutive 
   * values, as done in the function Pscf::PolymerTmpl::readParameters.
   *
   * Different text representations are used for bead and thread models,
   * and for linear and branched polymers:
   *
   * Thread model:
   *
   * In the thread model, if PolymerModel::isThread(), the text 
   * representation of an Edge for a branched polymer is:
   * \code
   *    monomerId length vertexId(0) vertexid(1)
   * \endcode
   * Here, length is a floating point number. For a linear polymer in the
   * thread model, the vertex ids are omitted, because the vertex ids for
   * block id of linear polymer must be equal to id and id + 1.
   *
   * Bead model:
   * 
   * In the bead model, if PolymerModel::isBead(), the text representation 
   * for a branched polymer is:
   * \code
   *    monomerId nBead vertexId(0) vertexid(1) ownsVertex(0) ownsVertex(1)
   * \endcode
   * Here, nBead is the integer number of beads in the chain, while 
   * ownsVertex(0) and ownsVertex(1) are boolean variables (0 or 1) that 
   * indicate whether the block owns vertex 0 and 1, respecively. The 
   * corresponding text representation for a bead-spring linear chain 
   * omits both the vertex ids and the ownsVertex flags, because values of 
   * these are set by convention.  In a linear bead-spring polymer, 
   * vertex ids for block number id are id and id+1, ownsVertex(1) == 1, 
   * while ownsVertex(0) == 1 for the first block (id = 0) and 
   * ownsVertex(0) = 0 for all other blocks (id > 0).
   *
   * \param in  input stream
   * \param block  Edge to be read from stream
   * \return  modified input stream
   */
   std::istream& operator >> (std::istream& in, Edge &block);

   /**
   * Output stream inserter (<<) for a Edge.
   *
   * Different text representations are used for linear and branched
   * polymers, as discussed in documentation for the stream extractor
   * (>>) operator. Vertex ids are output only for blocks of branched
   * polymers.
   *
   * \param out  output stream
   * \param block  Edge to be written to stream
   * \return modified output stream
   */
   std::ostream&
   operator << (std::ostream& out, const Edge &block);

   // Inline member functions

   /*
   * Get the id of this block.
   */
   inline int Edge::id() const
   {  return id_; }

   /*
   * Get the monomer type id.
   */
   inline int Edge::monomerId() const
   {  return monomerId_; }

   /*
   * Get the pair of associated vertex ids.
   */
   inline const Pair<int>& Edge::vertexIds() const
   {  return vertexIds_; }

   /*
   * Get id of an associated vertex.
   */
   inline int Edge::vertexId(int i) const
   {  return vertexIds_[i]; }

   /*
   * Get ownsVertex flag for an associated vertex (bead model only).
   */
   inline bool Edge::ownsVertex(int i) const
   {  
      UTIL_CHECK(PolymerModel::isBead());  
      return ownsVertex_[i]; 
   }

   /*
   * Get the number of beads in this block (bead model).
   */
   inline double Edge::nBead() const
   {
      UTIL_CHECK(PolymerModel::isBead());  
      return nBead_; 
   }

   /*
   * Get the length (number of monomers) in this block.
   */
   inline double Edge::length() const
   {
      UTIL_CHECK(PolymerModel::isThread());  
      return length_; 
   }

   /*
   * Get the polymer type (branched or linear).
   */
   inline PolymerType::Enum Edge::polymerType() const
   {  return polymerType_; }

   /*
   * Serialize to/from an archive.
   */
   template <class Archive>
   void Edge::serialize(Archive& ar, unsigned int)
   {
      ar & id_;
      ar & monomerId_;
      ar & vertexIds_;
      ar & length_;
   }

}
#endif
