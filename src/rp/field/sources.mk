rp_field_OBJS_=

ifdef PSCF_CPP
rp_field_CPP= \
   rp/field/CFields_c.cpp \
   rp/field/WFields_c.cpp \
   rp/field/Mask_c.cpp \
   rp/field/FieldIo_c.cpp \
   rp/field/Domain_c.cpp
rp_field_OBJS+= \
     $(addprefix $(BLD_DIR)/, $(rp_field_CPP:.cpp=.o))
endif

ifdef PSCF_CUDA
rp_field_CUDA= \
   rp/field/CFields_u.cu \
   rp/field/WFields_u.cu \
   rp/field/Mask_u.cu \
   rp/field/FieldIo_u.cu \
   rp/field/Domain_u.cu
rp_field_OBJS+= \
   $(addprefix $(BLD_DIR)/, $(rp_field_CUDA:.cu=.o))
endif
