
rp_fts_analyzer_OBJS=
rp_fts_analyzer_DEPS=

ifdef PSCF_CPP
  rp_fts_analyzer_CPP= \
    rp/fts/analyzer/Analyzer.cpp \
    rp/fts/analyzer/AnalyzerManager.cpp \
    rp/fts/analyzer/AverageAnalyzer.cpp \
    rp/fts/analyzer/AverageListAnalyzer.cpp \
    rp/fts/analyzer/BinaryChiDerivative.cpp \
    rp/fts/analyzer/BinaryStructureFactor.cpp \
    rp/fts/analyzer/ConcentrationDerivative.cpp \
    rp/fts/analyzer/CubicLengthDerivative.cpp \
    rp/fts/analyzer/ConcentrationWriter.cpp \
    rp/fts/analyzer/HamiltonianAnalyzer.cpp \
    rp/fts/analyzer/PerturbationDerivative.cpp \
    rp/fts/analyzer/StepLogger.cpp \
    rp/fts/analyzer/TrajectoryWriter.cpp \
    rp/fts/analyzer/FourthOrderParameter.cpp \
    rp/fts/analyzer/MaxOrderParameter.cpp \
    rp/fts/analyzer/AnalyzerFactory.cpp
  rp_fts_analyzer_OBJS+=\
       $(addprefix $(BLD_DIR)/, $(rp_fts_analyzer_CPP:.cpp=.o))
  rp_fts_analyzer_DEPS+=\
       $(addprefix $(BLD_DIR)/, $(rp_fts_analyzer_CPP:.cpp=.d))
endif

ifdef PSCF_CUDA
  rp_fts_analyzer_CU= \
    rp/fts/analyzer/Analyzer.cu \
    rp/fts/analyzer/AnalyzerManager.cu \
    rp/fts/analyzer/AverageAnalyzer.cu \
    rp/fts/analyzer/AverageListAnalyzer.cu \
    rp/fts/analyzer/BinaryChiDerivative.cu \
    rp/fts/analyzer/BinaryStructureFactor.cu \
    rp/fts/analyzer/ConcentrationDerivative.cu \
    rp/fts/analyzer/CubicLengthDerivative.cu \
    rp/fts/analyzer/ConcentrationWriter.cu \
    rp/fts/analyzer/HamiltonianAnalyzer.cu \
    rp/fts/analyzer/PerturbationDerivative.cu \
    rp/fts/analyzer/StepLogger.cu \
    rp/fts/analyzer/TrajectoryWriter.cu \
    rp/fts/analyzer/FourthOrderParameter.cu \
    rp/fts/analyzer/MaxOrderParameter.cu \
    rp/fts/analyzer/AnalyzerFactory.cu
  rp_fts_analyzer_OBJS+=\
      $(addprefix $(BLD_DIR)/, $(rp_fts_analyzer_CU:.cu=.ou))
  rp_fts_analyzer_DEPS+=\
      $(addprefix $(BLD_DIR)/, $(rp_fts_analyzer_CU:.cu=.du))
endif
