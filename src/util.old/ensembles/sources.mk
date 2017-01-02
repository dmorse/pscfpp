util_ensembles_=\
    util/ensembles/EnergyEnsemble.cpp \
    util/ensembles/BoundaryEnsemble.cpp \
    util/ensembles/SpeciesEnsemble.cpp 

util_ensembles_SRCS=$(addprefix $(SRC_DIR)/, $(util_ensembles_))
util_ensembles_OBJS=$(addprefix $(BLD_DIR)/, $(util_ensembles_:.cpp=.o))

