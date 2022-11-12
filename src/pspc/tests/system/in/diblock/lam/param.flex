System{
  Mixture{
     nMonomer  2
     monomers[
               1.0  
               1.0 
     ]
     nPolymer  1
     Polymer{
        type    linear
        nBlock  2
        blocks[
                0  0.5
                1  0.5
        ]
        phi     1.0
     }
     ds   0.01
  }
  Interaction{
     chi(  
          1   0   15.0
     )
  }
  Domain{
     mesh        32
     lattice     Lamellar   
     groupName   P_-1
  }
  AmIterator{
     maxItr 300
     epsilon 1e-10
     maxHist 10
     isFlexible   1
  }
}


     unitCell Lamellar   1.4814442100
