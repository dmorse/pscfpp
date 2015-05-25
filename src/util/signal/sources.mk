util_signal_=util/signal/Signal.cpp 

util_signal_SRCS=\
     $(addprefix $(SRC_DIR)/, $(util_signal_))
util_signal_OBJS=\
     $(addprefix $(BLD_DIR)/, $(util_signal_:.cpp=.o))
