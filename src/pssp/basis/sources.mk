pssp_basis_= \
     pssp/basis/groupFile.cpp

pssp_basis_SRCS=\
     $(addprefix $(SRC_DIR)/, $(pssp_basis_))
pssp_basis_OBJS=\
     $(addprefix $(BLD_DIR)/, $(pssp_basis_:.cpp=.o))

# Specialized rule for groupFile.o
$(BLD_DIR)/pssp/basis/groupFile.o: $(SRC_DIR)/pssp/basis/groupFile.cpp
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) $(INCLUDES) $(DEFINES) -D DAT_DIR=$(DAT_DIR) -c -o $@ $<
ifdef MAKEDEP
	$(MAKEDEP) $(INCLUDES) $(DEFINES) -D DAT_DIR=$(DAT_DIR) $(MAKE_DEPS) -S$(SRC_DIR) -B$(BLD_DIR) $<
endif

