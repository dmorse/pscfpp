pscf_tests_cpu_=pscf/tests/cpu/Test.cpp

pscf_tests_cpu_OBJS=\
     $(addprefix $(BLD_DIR)/, $(pscf_tests_cpu_:.cpp=.o))

pscf_tests_cpu_DEPS=\
     $(addprefix $(BLD_DIR)/, $(pscf_tests_cpu_:.cpp=.d))

