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
    chi  0   0   0.0
         1   0   25.0
         1   1   0.0
  }
  Domain{
    mesh              32  32  32
    lattice           Cubic  
    groupName         F_d_-3_m:1 
  }
  AmIterator{
    maxItr   800
    epsilon  1e-8
    maxHist  20
    isFlexible   1
  }
}

    unitCell  Cubic   6.9947873
