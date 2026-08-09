System{
  Mixture{
    nMonomer   2
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
    epsilon   1.0e-10
    maxItr    300
    maxHist   10
    verbose   1
    isFlexible  1
  }
}

