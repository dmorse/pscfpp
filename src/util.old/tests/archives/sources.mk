ifdef UTIL_MPI
util_tests_archives_=util/tests/archives/MpiTest.cc
else
util_tests_archives_=util/tests/archives/Test.cc
endif

util_tests_archives_SRCS=\
     $(addprefix $(SRC_DIR)/, $(util_tests_archives_))
util_tests_archives_OBJS=\
     $(addprefix $(BLD_DIR)/, $(util_tests_archives_:.cc=.o))

