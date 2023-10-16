pspg_tests_field_=pspg/tests/field/Test.cu

pspg_tests_field_OBJS=\
     $(addprefix $(BLD_DIR)/, $(pspg_tests_field_:.cu=.o))

