System{
  polymerModel  Bead
  Mixture{
     nMonomer  2
     monomers[
               0.1
               0.1
     ]
     nPolymer  1
     Polymer{
        nBlock  2
        blocks[
                0   50
                1   50
        ]
        phi     1.0
     }
  }
  Interaction{
     chi(
          1   0   0.15
     )
  }
  Domain{
    mesh         32
    lattice      lamellar 
    groupName    P_-1 
  }
  AmIterator{
     epsilon      1.0e-10
     maxItr       300
     maxHist      10
     verbose      1
     isFlexible   1
  }
}

    unitCell Lamellar   1.3835952906
