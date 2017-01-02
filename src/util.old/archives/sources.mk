
util_archives_=\
    util/archives/MemoryOArchive.cpp \
    util/archives/MemoryIArchive.cpp \
    util/archives/MemoryCounter.cpp \
    util/archives/BinaryFileOArchive.cpp \
    util/archives/BinaryFileIArchive.cpp \
    util/archives/TextFileOArchive.cpp \
    util/archives/TextFileIArchive.cpp \
    util/archives/XdrFileOArchive.cpp \
    util/archives/XdrFileIArchive.cpp 

util_archives_SRCS=$(addprefix $(SRC_DIR)/, $(util_archives_))
util_archives_OBJS=$(addprefix $(BLD_DIR)/, $(util_archives_:.cpp=.o))

