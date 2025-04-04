/*
* PSCF - Polymer Self-Consistent Field Theory
*
* Copyright 2016 - 2022, The Regents of the University of Minnesota
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
      length_(-1.0),
      vertexIds_(),
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
   * Set indices of associated vertices.
   */
   void Edge::setVertexIds(int vertexId0, int vertexId1)
   {
      vertexIds_[0] = vertexId0;
      vertexIds_[1] = vertexId1;
   }

   /*
   * Set the monomer id.
   */
   void Edge::setMonomerId(int monomerId)
   {  monomerId_ = monomerId; }

   /*
   * Set the length of this block.
   */
   void Edge::setLength(double length)
   {  length_ = length; }

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
      in >> block.monomerId_;
      in >> block.length_;
      if (block.polymerType_ == PolymerType::Branched) {
         in >> block.vertexIds_[0];
         in >> block.vertexIds_[1];
      }
      return in;
   }

   /*
   * Output a Edge to an ostream, without line breaks.
   */
   std::ostream& operator  << (std::ostream& out, 
                               Edge const & block)
   {
      out << "  " << block.monomerId_;
      out << "  ";
      out.setf(std::ios::scientific);
      out.width(20);
      out.precision(12);
      out << block.length_;
      if (block.polymerType_ == PolymerType::Branched) {
         out << "  " << block.vertexIds_[0];
         out << "  " << block.vertexIds_[1];
      }
      return out;
   }

}
