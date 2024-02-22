prdc_tests_crystal_=prdc/tests/crystal/Test.cpp

prdc_tests_crystal_OBJS=\
     $(addprefix $(BLD_DIR)/, $(prdc_tests_crystal_:.cpp=.o))

