pscf_math_= \
  pscf/math/LuSolver.cpp \
  pscf/math/TridiagonalSolver.cpp \
  pscf/math/IntVec.cpp \
  pscf/math/Field.cpp

pscf_math_SRCS=\
     $(addprefix $(SRC_DIR)/, $(pscf_math_))
pscf_math_OBJS=\
     $(addprefix $(BLD_DIR)/, $(pscf_math_:.cpp=.o))

