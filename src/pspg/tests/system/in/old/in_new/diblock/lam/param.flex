System{
  Mixture{
     nMonomer  2
     monomers  0   A   1.0  
               1   B   1.0 
     nPolymer  1
     Polymer{
        nBlock  2
        nVertex 3
        blocks  0  0  0  1  0.56
                1  1  1  2  0.44
        phi     1.0
     }
     ds   0.01
  }


  ChiInteraction{
     chi  0   0   0.0
          1   0   12.0
          1   1   0.0
  }
   
unitCell Lamellar   1.3835952906
mesh  	 40
groupName P_-1

  AmIterator{
   maxItr 100
   epsilon 1e-12
   maxHist 10
   isFlexible 1
  }

}
