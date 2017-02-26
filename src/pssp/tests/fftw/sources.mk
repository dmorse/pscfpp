pssp_tests_fftw_=pssp/tests/fftw/Test.cc

pssp_tests_fftw_SRCS=\
     $(addprefix $(SRC_DIR)/, $(pssp_tests_fftw_))
pssp_tests_fftw_OBJS=\
     $(addprefix $(BLD_DIR)/, $(pssp_tests_fftw_:.cc=.o))

