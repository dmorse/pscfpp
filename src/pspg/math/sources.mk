pspg_math_= \
     pspg/math/LinearAlgebra.cu \
     pspg/math/ParallelReductions.cu \
     pspg/math/KernelWrappers.cu \
     pspg/math/ThreadGrid.cu \

pspg_math_SRCS=\
     $(addprefix $(SRC_DIR)/, $(pspg_math_))
pspg_math_OBJS=\
     $(addprefix $(BLD_DIR)/, $(pspg_math_:.cu=.o))

