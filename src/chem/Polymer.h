#ifndef PSCF_CHEM_POLYMER_H
#define PSCF_CHEM_POLYMER_H

/*
* PFTS - Polymer Field Theory Simulator
*
* Copyright 2013, David Morse (morse012@.umn.edu)
* Distributed under the terms of the GNU General Public License.
*/

namespace Pfts{ 
namespace Chem{

   /**
   * Descriptor for an acyclic block polymer.
   */
   class Polymer 
   {

      /**
      * Initialize system parameters.
      */
      void readParameters(std::istream& in);

      Vertex& vertex(int id);
      PArray<Vertex>& ends();
      Array<Block>& blocks();

      /**
      * Computation plan.
      *
      * The links array is ordered in the order of the computation
      * plan, i.e., in an order in which the links must be solved
      * so that the intitial condition for each link is provided
      * by the solution of links that have already been solved.
      */
      Array<int>& plan();

      /**
      * Volume per molecule, in units of reference volume.
      */
      double volume();

      int  nVertex();

      int  nEnd();  //

      int  nBlock();  //

      void makePlan();

   private:

      int nVertex_;
      int nBlock_;
      DArray<Block> blocks_;
      DArray<Vertex> vertices_;
      DPArray<Vertex> ends_;
      DArray<int> plan_;

   };

} 
} 
#endif 
