
/*
* PSCF++ - Polymer Self-Consistent Field Theory Simulator
*
* Copyright 2013, David Morse (morse012@.umn.edu)
* Distributed under the terms of the GNU General Public License.
*/


   class Link
   {

   public:

      /**
      * Constructor.
      */
      void Link();

      /**
      * Create associations.
      */
      void init(Vertex& head, Vertex& tail, Block& block);

      /**
      * Create association with other links on head Vertex.
      */
      void init();

      /**
      * Solve modified diffusion equation.
      */
      virtual void compute() = 0;

      /**
      * Return partition function q(r) at tail vertex.
      */
      virtual Field& tailQFunction() = 0;

      const Block&  block() const;

      const Vertex& head() const;

      const Vertex& tail() const;

   private:

      Block*  blockPtr_;
      Vertex* headPtr_;
      Vertex* tailPtr_;
      DPArray<Link> sourceLinks_;
   };

