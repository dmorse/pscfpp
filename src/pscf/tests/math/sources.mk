pscf_tests_math_=pscf/tests/math/Test.cc

pscf_tests_math_OBJS=\
     $(addprefix $(BLD_DIR)/, $(pscf_tests_math_:.cc=.o))

