
rp_fts_trajectory_OBJS=
rp_fts_trajectory_DEPS=

ifdef PSCF_CPP
  rp_fts_trajectory_CPP= \
    rp/fts/trajectory/TrajectoryReader.cpp \
    rp/fts/trajectory/RGridTrajectoryReader.cpp \
    rp/fts/trajectory/TrajectoryReaderFactory.cpp
  rp_fts_trajectory_OBJS+=\
       $(addprefix $(BLD_DIR)/, $(rp_fts_trajectory_CPP:.cpp=.o))
  rp_fts_trajectory_DEPS+=\
       $(addprefix $(BLD_DIR)/, $(rp_fts_trajectory_CPP:.cpp=.d))
endif

ifdef PSCF_CUDA
  rp_fts_trajectory_CU= \
    rp/fts/trajectory/TrajectoryReader.cu \
    rp/fts/trajectory/RGridTrajectoryReader.cu \
    rp/fts/trajectory/TrajectoryReaderFactory.cu
  rp_fts_trajectory_OBJS+=\
      $(addprefix $(BLD_DIR)/, $(rp_fts_trajectory_CU:.cu=.ou))
  rp_fts_trajectory_DEPS+=\
      $(addprefix $(BLD_DIR)/, $(rp_fts_trajectory_CU:.cu=.du))
endif
