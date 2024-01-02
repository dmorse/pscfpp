pspc_simulate_analyzer_= \
  pspc/simulate/analyzer/Analyzer.cpp \
  pspc/simulate/analyzer/AverageListAnalyzer.cpp \
  pspc/simulate/analyzer/AnalyzerManager.cpp \
  pspc/simulate/analyzer/AnalyzerFactory.cpp \
  pspc/simulate/analyzer/TrajectoryWriter.cpp \
  pspc/simulate/analyzer/HamiltonianAnalyzer.cpp \
  pspc/simulate/analyzer/HamiltonianAutoCorr.cpp \
  pspc/simulate/analyzer/BinaryStructureFactorGrid.cpp \
  pspc/simulate/analyzer/StepLogger.cpp 
  
pspc_simulate_analyzer_OBJS=\
     $(addprefix $(BLD_DIR)/, $(pspc_simulate_analyzer_:.cpp=.o))

