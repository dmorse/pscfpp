System{
  Mixture{
    nMonomer  2
    monomers  2.0  
              1.0 
    nPolymer  1
    Polymer{
      type    linear
      nBlock  2
      blocks  0  0.25
              1  0.75
      phi     1.0
    }
    ds   0.01
  }
  Interaction{
    chi(  
         1   0   25.0
    ) 
  }
  Domain{
    mesh            32  32  32
    lattice         cubic  
    groupName       F_d_-3_m:1 
  }
  AmIteratorBasis{
    epsilon  1.0e-8
    maxItr   800
    maxHist  20
    verbose  1
    isFlexible   1
  }
}

