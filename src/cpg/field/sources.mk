cpg_field_= 
  #cpg/field/FieldIo.cu \
  #cpg/field/Domain.cu \
  #cpg/field/WFields.cu \
  #cpg/field/CFields.cu

cpg_field_OBJS=\
     $(addprefix $(BLD_DIR)/, $(cpg_field_:.cu=.o))

