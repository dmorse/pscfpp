/*
* PSCF - Polymer Self-Consistent Field Theory
*
* Copyright 2015 - 2025, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include "Edge.h"

namespace Pscf
{

   /*
   * Constructor.
   */
   Edge::Edge()
    : id_(-1),
      monomerId_(-1),
      nBead_(-1),
      length_(-1.0),
      vertexIds_(),
      ownsVertex_(),
      polymerType_(PolymerType::Branched)
   {
      // Initialize vertex ids to null value -1.
      vertexIds_[0] = -1;
      vertexIds_[1] = -1;
   }

   /*
   * Destructor (virtual)
   */
   Edge::~Edge()
   {}

   /*
   * Set the id for this block.
   */
   void Edge::setId(int id)
   {  id_ = id; }

   /*
   * Set the monomer type id.
   */
   void Edge::setMonomerId(int monomerId)
   {  monomerId_ = monomerId; }

   /*
   * Set the indices of the two associated vertices.
   */
   void Edge::setVertexIds(int vertexId0, int vertexId1)
   {
      vertexIds_[0] = vertexId0;
      vertexIds_[1] = vertexId1;
   }

   /*
   * Set ownership of the two associated vertices.
   */
   void Edge::setVertexOwnership(bool own0, bool own1)
   {
      ownsVertex_[0] = own0;
      ownsVertex_[1] = own1;
   }

   /*
   * Set the number of beads in the block.
   */
   void Edge::setNBead(int nBead)
   {
      UTIL_CHECK(PolymerModel::isBead());  
      nBead_ = nBead; 
   }

   /*
   * Set the length of this block.
   */
   void Edge::setLength(double length)
   {  
      UTIL_CHECK(PolymerModel::isThread());  
      length_ = length; 
   }

   /*
   * Set the type of the polymer containing this block.
   */
   void Edge::setPolymerType(PolymerType::Enum type)
   {  polymerType_ = type; }

   /*
   * Extract a Edge from an istream.
   */
   std::istream& operator >> (std::istream& in, Edge& block)
   {
      // Read monomer type id
      in >> block.monomerId_;

      // Read length or nBead
      if (PolymerModel::isThread()) {
         in >> block.length_;
      } else 
      if (PolymerModel::isBead()) {
         in >> block.nBead_;
      }

      // For branched polyimers, read topology information 
      if (block.polymerType_ == PolymerType::Branched) {
         in >> block.vertexIds_[0];
         in >> block.vertexIds_[1];
         if (PolymerModel::isBead()) {
            in >> block.ownsVertex_[0];
            in >> block.ownsVertex_[1];
         }
      }

      return in;
   }

   /*
   * Output a Edge to an ostream, without line breaks.
   */
   std::ostream& operator  << (std::ostream& out, 
                               Edge const & block)
   {
      // Write monomer type id
      out << "  " << block.monomerId_;
      out << "  ";

      // Write length or nBead
      if (PolymerModel::isThread()) {
         out.setf(std::ios::scientific);
         out.width(20);
         out.precision(12);
         out << block.length_;
      } else
      if (PolymerModel::isBead()) {
         out << block.nBead_;
      }

      // For branched polymer, write topology information
      if (block.polymerType_ == PolymerType::Branched) {
         out << "  " << block.vertexIds_[0];
         out << "  " << block.vertexIds_[1];
         if (PolymerModel::isBead()) {
            out << "  " << block.ownsVertex_[0];
            out << "  " << block.ownsVertex_[1];
         }
      }

      return out;
   }

}
