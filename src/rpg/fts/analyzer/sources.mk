rpg_fts_analyzer_= \
  rpg/fts/analyzer/Analyzer.cu \
  rpg/fts/analyzer/AverageListAnalyzer.cu \
  rpg/fts/analyzer/AnalyzerManager.cu \
  rpg/fts/analyzer/AnalyzerFactory.cu \
  rpg/fts/analyzer/TrajectoryWriter.cu \
  rpg/fts/analyzer/ConcentrationWriter.cu \
  rpg/fts/analyzer/HamiltonianAnalyzer.cu \
  rpg/fts/analyzer/BinaryStructureFactorGrid.cu \
  rpg/fts/analyzer/StepLogger.cu \
  rpg/fts/analyzer/ThermoDerivativeAnalyzer.cu \
  rpg/fts/analyzer/PerturbationDerivative.cu \
  rpg/fts/analyzer/ChiDerivative.cu \
  rpg/fts/analyzer/ConcentrationDerivative.cu \
  rpg/fts/analyzer/MaxOrderParameter.cu \
  rpg/fts/analyzer/FourthOrderParameter.cu
  
rpg_fts_analyzer_OBJS=\
     $(addprefix $(BLD_DIR)/, $(rpg_fts_analyzer_:.cu=.o))

