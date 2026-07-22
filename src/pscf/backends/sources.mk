pscf_backends_CPP =\
   pscf/backends/CPT.cpp

pscf_backends_OBJS=\
     $(addprefix $(BLD_DIR)/, $(pscf_backends_CPP:.cpp=.o))

ifdef PSCF_CUDA
pscf_backends_CU = pscf/backends/CUT.cu 

pscf_backends_OBJS+=\
     $(addprefix $(BLD_DIR)/, $(pscf_backends_CU:.cu=.o))
endif
