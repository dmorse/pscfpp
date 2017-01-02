
util_format_=util/format/Bool.cpp \
    util/format/Dbl.cpp util/format/Format.cpp \
    util/format/Int.cpp util/format/Lng.cpp \
    util/format/Str.cpp util/format/write.cpp 

util_format_SRCS=$(addprefix $(SRC_DIR)/, $(util_format_))
util_format_OBJS=$(addprefix $(BLD_DIR)/, $(util_format_:.cpp=.o))

