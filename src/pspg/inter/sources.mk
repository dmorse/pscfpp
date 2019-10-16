pspg_inter_= \
  pspg/inter/ChiInteraction.cu 

pspg_inter_SRCS=\
     $(addprefix $(SRC_DIR)/, $(pspg_inter_))
pspg_inter_OBJS=\
     $(addprefix $(BLD_DIR)/, $(pspg_inter_:.cu=.o))

