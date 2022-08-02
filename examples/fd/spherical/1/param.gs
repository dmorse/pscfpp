System{
  Mixture{
     nMonomer  2
     monomers[  
            1.0  
            1.0 
     ]
     nPolymer  2
     Polymer{
        type      linear
        nBlock    2
        blocks[
             0  0.125
             1  0.875
        ]
        mu        5.60908591992e+00   
     }
     Polymer{
        type      linear
        nBlock    1
        blocks[
            1  1.000
        ]
        mu        -3.6515433392e-03  
     }
     vMonomer     0.045787
     ds           0.005
  }
  ChiInteraction{
     chi(
         0  1    88.5
     )
  }
  Domain{
     mode      Spherical
     isShell           0
     xMax           4.00
     nx              401
  }
  NrIterator{
     epsilon   0.0000001
  }
  hasSweep     1
  MuSweep{
     ns                    200
     baseFileName          gs/
     homogeneousMode         1
     dMu              +0.40000
                       0.00000
  }
}

