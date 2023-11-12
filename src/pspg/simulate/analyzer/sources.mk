pspg_simulate_analyzer_= \
  pspg/simulate/analyzer/Analyzer.cu \
  pspg/simulate/analyzer/AverageListAnalyzer.cu \
  pspg/simulate/analyzer/AnalyzerManager.cu \
  pspg/simulate/analyzer/AnalyzerFactory.cu \
  pspg/simulate/analyzer/TrajectoryWriter.cu \
  pspg/simulate/analyzer/HamiltonianAnalyzer.cu \
  pspg/simulate/analyzer/BinaryStructureFactorGrid.cu
  
pspg_simulate_analyzer_OBJS=\
     $(addprefix $(BLD_DIR)/, $(pspg_simulate_analyzer_:.cu=.o))

