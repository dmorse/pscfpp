pspg_field_= \
  pspg/field/FFTBatched.cu \
  pspg/field/FieldIo.cu \
  pspg/field/WFieldContainer.cu \
  pspg/field/CFieldContainer.cu \
  pspg/field/RFieldComparison.cu \
  pspg/field/KFieldComparison.cu \
  pspg/field/Domain.cu 

pspg_field_SRCS=\
	  $(addprefix $(SRC_DIR)/, $(pspg_field_))
pspg_field_OBJS=\
	  $(addprefix $(BLD_DIR)/, $(pspg_field_:.cu=.o))

