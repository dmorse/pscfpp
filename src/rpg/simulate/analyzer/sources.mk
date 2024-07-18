rpg_simulate_analyzer_= \
  rpg/simulate/analyzer/Analyzer.cu \
  rpg/simulate/analyzer/AverageListAnalyzer.cu \
  rpg/simulate/analyzer/AnalyzerManager.cu \
  rpg/simulate/analyzer/AnalyzerFactory.cu \
  rpg/simulate/analyzer/TrajectoryWriter.cu \
  rpg/simulate/analyzer/HamiltonianAnalyzer.cu \
  rpg/simulate/analyzer/BinaryStructureFactorGrid.cu \
  rpg/simulate/analyzer/StepLogger.cu \
  rpg/simulate/analyzer/ThermoDerivativeAnalyzer.cu \
  rpg/simulate/analyzer/PerturbationDerivative.cu \
  rpg/simulate/analyzer/ChiDerivative.cu \
  rpg/simulate/analyzer/ConcentrationDerivative.cu
  
rpg_simulate_analyzer_OBJS=\
     $(addprefix $(BLD_DIR)/, $(rpg_simulate_analyzer_:.cu=.o))

