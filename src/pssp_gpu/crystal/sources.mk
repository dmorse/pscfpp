pspg_crystal_= \
  pspg/crystal/UnitCell1.cu \
  pspg/crystal/UnitCell2.cu \
  pspg/crystal/UnitCell3.cu \
  pspg/crystal/shiftToMinimum.cu

pspg_crystal_SRCS=\
     $(addprefix $(SRC_DIR)/, $(pspg_crystal_))
pspg_crystal_OBJS=\
     $(addprefix $(BLD_DIR)/, $(pspg_crystal_:.cu=.o))

