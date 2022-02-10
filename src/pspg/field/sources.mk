pspg_field_= \
  pspg/field/FFT.cu \
  pspg/field/FFTBatched.cu \
  pspg/field/FieldIo.cu \
  pspg/field/BFieldComparison.cu \
  pspg/field/RFieldComparison.cu \
  pspg/field/KFieldComparison.cu \
  pspg/field/Domain.cu

pspg_field_SRCS=\
	  $(addprefix $(SRC_DIR)/, $(pspg_field_))
pspg_field_OBJS=\
	  $(addprefix $(BLD_DIR)/, $(pspg_field_:.cu=.o))

