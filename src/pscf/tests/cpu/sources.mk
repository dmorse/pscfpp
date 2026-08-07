pscf_tests_cpu_=pscf/tests/cpu/Test.cu

pscf_tests_cpu_OBJS=\
     $(addprefix $(BLD_DIR)/, $(pscf_tests_cpu_:.cu=.o))

pscf_tests_cpu_DEPS=\
     $(addprefix $(BLD_DIR)/, $(pscf_tests_cpu_:.cu=.d))

