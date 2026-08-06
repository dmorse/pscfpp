rp_field_OBJS_=
rp_field_DEPS_=

ifdef PSCF_CPP
  rp_field_CPP= \
     rp/field/CFields_c.cpp \
     rp/field/WFields_c.cpp \
     rp/field/Mask_c.cpp \
     rp/field/FieldIo_c.cpp \
     rp/field/Domain_c.cpp
  rp_field_OBJS+= \
       $(addprefix $(BLD_DIR)/, $(rp_field_CPP:.cpp=.o))
  rp_field_DEPS+= \
       $(addprefix $(BLD_DIR)/, $(rp_field_CPP:.cpp=.d))
endif

ifdef PSCF_CUDA
  rp_field_CU= \
     rp/field/CFields_u.cu \
     rp/field/WFields_u.cu \
     rp/field/Mask_u.cu \
     rp/field/FieldIo_u.cu \
     rp/field/Domain_u.cu
  rp_field_OBJS+= \
     $(addprefix $(BLD_DIR)/, $(rp_field_CU:.cu=.o))
  rp_field_DEPS+= \
       $(addprefix $(BLD_DIR)/, $(rp_field_CU:.cu=.d))
endif
