 format  1  0

MONOMERS
N_monomer
                   3
kuhn
    6.0700000000E+00    6.0700000000E+00    6.0700000000E+00

CHAINS
N_chain
                   1
N_block
                   2
block_monomer
                   1                   2
block_length
    3.5000000000E-01    6.5000000000E-01

SOLVENTS
N_solvent
                   1
solvent_monomer
                   3
solvent_size
    2.0000000000E-02

COMPOSITION
ensemble
                   0
phi_chain
    5.0000000000E-01
phi_solvent
    5.0000000000E-01

INTERACTION
interaction_type
               'chi'
chi
    3.0000000000E+01
    3.0000000000E+01    1.0000000000E+02

UNIT_CELL
dim
                   1
crystal_system
          'lamellar'
N_cell_param
                   1
cell_param
   1.62132470000E+01

DISCRETIZATION
ngrid
                 100
chain_step
    1.0000000000E-02

BASIS
group_name
                '-1'

ITERATE
input_filename
          'ref/w.bf'
output_prefix
              'out/'
max_itr
                 100
error_max
    1.0000000000E-08
domain
                   T
itr_algo
                'NR'
N_cut
                  51

FIELD_TO_RGRID
input_filename
           'out/rho'
output_filename
         'out/rho.bf'

FINISH
