System{
  Mixture{
     nMonomer  2
     monomers  A   1.0  
               B   1.0 
     nPolymer  1
     Polymer{
        type    linear
        nBlock  2
        blocks  0  0.5
                1  0.5
        phi     1.0
     }
     ds   0.01
  }


  ChiInteraction{
     chi  0   0   0.0
          1   0   15.0
          1   1   0.0
  }
   
unitCell Lamellar   1.3835952907
mesh  	 32
groupName P_-1

  AmIterator{
   maxItr 300
   epsilon 1e-14
   maxHist 5
   isFlexible 1
  }

}
