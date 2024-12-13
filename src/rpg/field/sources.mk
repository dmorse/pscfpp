rpg_field_= \
  rpg/field/FieldIo.cu \
  rpg/field/WFieldContainer.cu \
  rpg/field/CFieldContainer.cu \
  rpg/field/Domain.cu 

rpg_field_OBJS=\
	  $(addprefix $(BLD_DIR)/, $(rpg_field_:.cu=.o))

