rp_field_OBJS_=
rp_field_DEPS_=

ifdef PSCF_CPP
  rp_field_CPP= \
     rp/field/CFields.cpp \
     rp/field/WFields.cpp \
     rp/field/Mask.cpp \
     rp/field/FieldIo.cpp \
     rp/field/Domain.cpp
  rp_field_OBJS+= \
       $(addprefix $(BLD_DIR)/, $(rp_field_CPP:.cpp=.o))
  rp_field_DEPS+= \
       $(addprefix $(BLD_DIR)/, $(rp_field_CPP:.cpp=.d))
endif

ifdef PSCF_CUDA
  rp_field_CU= \
     rp/field/CFields.cu \
     rp/field/WFields.cu \
     rp/field/Mask.cu \
     rp/field/FieldIo.cu \
     rp/field/Domain.cu
  rp_field_OBJS+= \
     $(addprefix $(BLD_DIR)/, $(rp_field_CU:.cu=.ou))
  rp_field_DEPS+= \
     $(addprefix $(BLD_DIR)/, $(rp_field_CU:.cu=.du))
endif
