
/*
* PSCF++ - Polymer Self-Consistent Field Theory Simulator
*
* Copyright 2013, David Morse (morse012@.umn.edu)
* Distributed under the terms of the GNU General Public License.
*/


   /**
   * A scalar function of position.
   */
   class Field
   {
   public:
   
      Field();
     
      /**
      * Clear all field values.
      */ 
      void clear();
   
      Array<double>& grid();
      Array<double>& coeffs();
   
      const Array<double>& grid() const;
      const Array<double>& coeff() const;
   
   protected:
   
      virtual void computeBasis() = 0;
      virtual void computeGrid() = 0;
   
      bool hasGrid_;
      bool hasCoeff_;
   
   };
