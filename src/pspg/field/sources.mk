pspg_field_= \
  pspg/field/FFTBatched.cu \
  pspg/field/FieldIo.cu \
  pspg/field/WFieldContainer.cu \
  pspg/field/CFieldContainer.cu \
  pspg/field/Domain.cu 

pspg_field_OBJS=\
	  $(addprefix $(BLD_DIR)/, $(pspg_field_:.cu=.o))

