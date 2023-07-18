pspc_simulate_analyzer_= \
  pspc/simulate/analyzer/Analyzer.cpp \
  pspc/simulate/analyzer/AverageListAnalyzer.cpp \
  pspc/simulate/analyzer/AnalyzerManager.cpp \
  pspc/simulate/analyzer/AnalyzerFactory.cpp \
  pspc/simulate/analyzer/TrajectoryWriter.cpp \
  pspc/simulate/analyzer/McHamiltonianAnalyzer.cpp 
  
pspc_simulate_analyzer_SRCS=\
     $(addprefix $(SRC_DIR)/, $(pspc_simulate_analyzer_))
pspc_simulate_analyzer_OBJS=\
     $(addprefix $(BLD_DIR)/, $(pspc_simulate_analyzer_:.cpp=.o))

