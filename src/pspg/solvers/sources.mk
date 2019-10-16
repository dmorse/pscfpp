
pspg_solvers_=

pspg_solvers_SRCS=\
     $(addprefix $(SRC_DIR)/, $(pspg_solvers_))
pspg_solvers_OBJS=\
     $(addprefix $(BLD_DIR)/, $(pspg_solvers_:.cpp=.o))

