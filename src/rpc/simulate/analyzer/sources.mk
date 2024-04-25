rpc_simulate_analyzer_= \
  rpc/simulate/analyzer/Analyzer.cpp \
  rpc/simulate/analyzer/AverageListAnalyzer.cpp \
  rpc/simulate/analyzer/AnalyzerManager.cpp \
  rpc/simulate/analyzer/AnalyzerFactory.cpp \
  rpc/simulate/analyzer/TrajectoryWriter.cpp \
  rpc/simulate/analyzer/HamiltonianAnalyzer.cpp \
  rpc/simulate/analyzer/HamiltonianAutoCorr.cpp \
  rpc/simulate/analyzer/BinaryStructureFactorGrid.cpp \
  rpc/simulate/analyzer/StepLogger.cpp \
  rpc/simulate/analyzer/LinearResponseAnalyzer.cpp \
  
rpc_simulate_analyzer_OBJS=\
     $(addprefix $(BLD_DIR)/, $(rpc_simulate_analyzer_:.cpp=.o))

