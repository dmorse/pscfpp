/*
* PSCF++ - Polymer Self-Consistent Field Theory Simulator
*
* Copyright 2013, David Morse (morse012@.umn.edu)
* Distributed under the terms of the GNU General Public License.
*/

   /**
   * An acyclic block copolymer species.
   */
   class Polymer : public Species
   {

      /**
      * Initialize system parameters.
      */
      void readParameters(std::istream& in);

      /**
      * Calculate partition function.
      */
      void compute();

      Array<Vertex>&  vertices();
      PArray<Vertex>& ends();
      PArray<Vertex>& junctions();
      Array<Block>&   blocks();

      /**
      * Array of all links.
      *
      * The links array is ordered in the order of the computation
      * plan, i.e., in an order in which the links must be solved
      * so that the intitial condition for each link is provided
      * by the solution of links that have already been solved.
      */
      PArray<Link>& links();

      /**
      * Volume per molecule, in units of reference volume.
      */
      double volume();
      int    nVertex();
      int    nEnd();  //
      int    nBlock();  //

   protected:

      void makePlan();

   private:

      unsigned int    nVertex_;
      unsigned int    nBlock_;
      DArray<Vertex>  vertices_;
      DPArray<Vertex> ends_;
      DArray<Block>   blocks_;
      DPArray<Link>   links_;

   };
