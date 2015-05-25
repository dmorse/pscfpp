
util_xmltag_=\
    util/xmltag/XmlBase.cpp \
    util/xmltag/XmlXmlTag.cpp \
    util/xmltag/XmlAttribute.cpp \
    util/xmltag/XmlStartTag.cpp \
    util/xmltag/XmlEndTag.cpp 


util_xmltag_SRCS=$(addprefix $(SRC_DIR)/, $(util_xmltag_))
util_xmltag_OBJS=$(addprefix $(BLD_DIR)/, $(util_xmltag_:.cpp=.o))

