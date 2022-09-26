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
  Interaction{
     chi(
         0  1    88.5
     )
  }
  Domain{
     mode      Spherical
     xMax           4.00
     nx              401
  }
  Iterator{
     epsilon   0.0000001
  }
  Sweep{
     ns                    200
     baseFileName          gs/
     homogeneousMode         1
     nParameter              1
     parameters[
        mu_polymer  0    +0.40000
     ]
  }
}

  MuSweep{
     ns                    200
     baseFileName          gs/
     homogeneousMode         1
     dMu              +0.40000
                       0.00000
  }

