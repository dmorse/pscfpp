pspg_solvers_= \
     pspg/solvers/WaveList.cu \
     pspg/solvers/Propagator.cu \
     pspg/solvers/Block.cu \
     pspg/solvers/Solvent.cu \
     pspg/solvers/Polymer.cu \
     pspg/solvers/Mixture.cu \

pspg_solvers_SRCS=\
     $(addprefix $(SRC_DIR)/, $(pspg_solvers_))
pspg_solvers_OBJS=\
     $(addprefix $(BLD_DIR)/, $(pspg_solvers_:.cu=.o))

