include $(SRC_DIR)/pfts/chem/sources.mk
include $(SRC_DIR)/pfts/solvers/sources.mk

pfts_= $(pfts_chem_) $(pfts_solvers_) 

pfts_SRCS=\
     $(addprefix $(SRC_DIR)/, $(pfts_))
pfts_OBJS=\
     $(addprefix $(BLD_DIR)/, $(pfts_:.cpp=.o))

$(pfts_LIB): $(pfts_OBJS)
	$(AR) rcs $(pfts_LIB) $(pfts_OBJS)

