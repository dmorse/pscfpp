rpg_system_= \
   rpg/system/System.cu \
   rpg/system/SystemConstRef.cu

rpg_system_OBJS=\
     $(addprefix $(BLD_DIR)/, $(rpg_system_:.cu=.o))

