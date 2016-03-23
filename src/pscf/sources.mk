include $(SRC_DIR)/pscf/chem/sources.mk
include $(SRC_DIR)/pscf/inter/sources.mk
include $(SRC_DIR)/pscf/math/sources.mk

pscf_= \
  $(pscf_chem_) $(pscf_inter_) $(pscf_math_) 

pscf_SRCS=\
     $(addprefix $(SRC_DIR)/, $(pscf_))
pscf_OBJS=\
     $(addprefix $(BLD_DIR)/, $(pscf_:.cpp=.o))

$(pscf_LIB): $(pscf_OBJS)
	$(AR) rcs $(pscf_LIB) $(pscf_OBJS)

