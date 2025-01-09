rpc_fts_analyzer_= \
  rpc/fts/analyzer/Analyzer.cpp \
  rpc/fts/analyzer/AverageListAnalyzer.cpp \
  rpc/fts/analyzer/AnalyzerManager.cpp \
  rpc/fts/analyzer/AnalyzerFactory.cpp \
  rpc/fts/analyzer/TrajectoryWriter.cpp \
  rpc/fts/analyzer/ConcentrationWriter.cpp \
  rpc/fts/analyzer/HamiltonianAnalyzer.cpp \
  rpc/fts/analyzer/HamiltonianAutoCorr.cpp \
  rpc/fts/analyzer/BinaryStructureFactorGrid.cpp \
  rpc/fts/analyzer/StepLogger.cpp \
  rpc/fts/analyzer/ThermoDerivativeAnalyzer.cpp \
  rpc/fts/analyzer/PerturbationDerivative.cpp \
  rpc/fts/analyzer/ChiDerivative.cpp \
  rpc/fts/analyzer/ConcentrationDerivative.cpp \
  rpc/fts/analyzer/MaxOrderParameter.cpp \
  rpc/fts/analyzer/FourthOrderParameter.cpp
  
  
rpc_fts_analyzer_OBJS=\
     $(addprefix $(BLD_DIR)/, $(rpc_fts_analyzer_:.cpp=.o))

