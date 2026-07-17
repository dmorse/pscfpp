rp_field_CPP= \
  rp/field/FieldIo.cpp \
  rp/field/Domain.cpp \
  rp/field/WFields.cpp \
  rp/field/CFields.cpp \
  rp/field/Mask.cpp 

rp_field_CUDA= \
  rp/field/FieldIo.cu \
  rp/field/Domain.cu \
  rp/field/WFields.cu \
  rp/field/CFields.cu \
  rp/field/Mask.cu 

rp_field_OBJS_=
ifdef PSCF_CPP
rp_field_OBJS+= \
     $(addprefix $(BLD_DIR)/, $(rp_field_CPP:.cpp=.o))
endif
ifdef PSCF_CUDA
rp_field_OBJS+= \
     $(addprefix $(BLD_DIR)/, $(rp_field_CUDA:.cu=.o))
endif

