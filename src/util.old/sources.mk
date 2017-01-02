include $(SRC_DIR)/util/misc/sources.mk
include $(SRC_DIR)/util/format/sources.mk
include $(SRC_DIR)/util/containers/sources.mk
include $(SRC_DIR)/util/mpi/sources.mk
include $(SRC_DIR)/util/signal/sources.mk
include $(SRC_DIR)/util/param/sources.mk
include $(SRC_DIR)/util/math/sources.mk
include $(SRC_DIR)/util/space/sources.mk
include $(SRC_DIR)/util/random/sources.mk
include $(SRC_DIR)/util/boundary/sources.mk
include $(SRC_DIR)/util/crystal/sources.mk
include $(SRC_DIR)/util/ensembles/sources.mk
include $(SRC_DIR)/util/accumulators/sources.mk
include $(SRC_DIR)/util/archives/sources.mk
include $(SRC_DIR)/util/xmltag/sources.mk

util_=$(util_misc_) $(util_format_) \
    $(util_containers_) $(util_mpi_) \
    $(util_signal_) $(util_param_) $(util_math_) \
    $(util_space_) $(util_random_) $(util_boundary_) \
    $(util_crystal_) $(util_ensembles_) \
    $(util_accumulators_) $(util_archives_) \
    $(util_xmltag_)

util_SRCS=\
     $(addprefix $(SRC_DIR)/, $(util_))
util_OBJS=\
     $(addprefix $(BLD_DIR)/, $(util_:.cpp=.o))

$(util_LIB): $(util_OBJS)
	$(AR) rcs $(util_LIB) $(util_OBJS)

