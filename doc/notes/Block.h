
/*
* PSCF++ - Polymer Self-Consistent Field Theory Simulator
*
* Copyright 2013, David Morse (morse012@.umn.edu)
* Distributed under the terms of the GNU General Public License.
*/

   /**
   * A linear homopolymer block within a block copolymer.
   */
   class Block
   {
   
      Block();
   
      setVertices(Vertex& v1, Vertex& v2);
      setId(unsigned int id);
      setMonomer(Monomer& monomer);
      setLength(double length);
      computeNStep(double dS0);
   
      /**
      * Return unnormalized concentration field.
      */
      virtual Field& rho() = 0;

      unsigned int id() const;
      const Link& link1() const;
      const Link& link2() const;
      const Monomer& monomer() const;
      double length() const;
      double ds() const;
      unsigned int nStep() const;
    
   private:
   
      Link     link1_;
      Link     link2_;
      double   length_;
      double   ds_;
      Monomer* monomerPtr_
      int      nStep_;
      int      id_;

   };
