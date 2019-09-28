
pssp_gpu_solvers_=

pssp_gpu_solvers_SRCS=\
     $(addprefix $(SRC_DIR)/, $(pssp_gpu_solvers_))
pssp_gpu_solvers_OBJS=\
     $(addprefix $(BLD_DIR)/, $(pssp_gpu_solvers_:.cpp=.o))

