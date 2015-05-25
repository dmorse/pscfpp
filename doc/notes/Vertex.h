
/*
* PSCF++ - Polymer Self-Consistent Field Theory Simulator
*
* Copyright 2013, David Morse (morse012@.umn.edu)
* Distributed under the terms of the GNU General Public License.
*/

   /**
   * A junction or chain end in a block copolymer.
   */
   class Vertex
   {
   public:
   
      void addBlock(Block& block);
      int degree();
      Link& inLink(unsigned int);
      Link& outLink(unsigned int);
   
   private:
   
      int degree_;
      std::vector<Link*> inLinks_;
      std::vector<Link*> outLinks_;
   
   };
