pspg_tests_iterator_=pspg/tests/iterator/Test.ccu

pspg_tests_iterator_SRCS=\
     $(addprefix $(SRC_DIR)/, $(pspg_tests_iterator_))
pspg_tests_iterator_OBJS=\
     $(addprefix $(BLD_DIR)/, $(pspg_tests_iterator_:.ccu=.o))

