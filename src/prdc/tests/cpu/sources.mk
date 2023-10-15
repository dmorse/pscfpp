prdc_tests_cpu_=prdc/tests/cpu/Test.cc

prdc_tests_cpu_OBJS=\
     $(addprefix $(BLD_DIR)/, $(prdc_tests_cpu_:.cc=.o))

