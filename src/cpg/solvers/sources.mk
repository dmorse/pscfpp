cpg_solvers_= 
  #cpg/solvers/Propagator.cu \
  #cpg/solvers/Block.cu \
  #cpg/solvers/Polymer.cu \
  #cpg/solvers/Solvent.cu \
  #cpg/solvers/Mixture.cu

cpg_solvers_OBJS=\
     $(addprefix $(BLD_DIR)/, $(cpg_solvers_:.cu=.o))

