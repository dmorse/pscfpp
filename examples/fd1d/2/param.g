System{
  Mixture{
     nMonomer  2
     monomers  0   A   1.0  
               1   B   1.0 
     nPolymer  2
     Polymer{
        nBlock    2
        nVertex   3
        blocks    0  0  0  1  0.125
                  1  1  1  2  0.875
        ensemble  open
        mu        5.74775182858
     }
     Polymer{
        nBlock    1
        nVertex   2
        blocks    0  1  0  1  1.000
        ensemble  open
        mu        -0.0042051012549
     }
     ds   0.005
  }
  ChiInteraction{
     chi   0  1    88.5
           0  0     0.0
           1  1     0.0
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
}

        phi     0.018223
        phi     0.981777
