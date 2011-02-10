#include $(MK)/header.mk

TARGETS := libHemeLbNet.$(LIBEXT)
SRCS := net.cc \
        PhasedBroadcast.cc

$(TARGETS)_DEPS := $(subst .cc,.$(OBJEXT), $(SRCS))

INCLUDES_$(d) := $(INCLUDES_$(parent))

include $(MK)/footer.mk
