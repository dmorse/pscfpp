rpg_solvers_= \
     rpg/solvers/Propagator.cu \
     rpg/solvers/Block.cu \
     rpg/solvers/Solvent.cu \
     rpg/solvers/Polymer.cu \
     rpg/solvers/Mixture.cu \

rpg_solvers_OBJS=\
     $(addprefix $(BLD_DIR)/, $(rpg_solvers_:.cu=.o))

