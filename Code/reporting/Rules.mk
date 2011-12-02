include $(MK)/header.mk

TARGETS := libHemeLbReporting.$(LIBEXT)
SRCS := Timers.cc Reporter.cc

$(TARGETS)_DEPS := $(subst .cc,.$(OBJEXT), $(SRCS))

INCLUDES_$(d) := $(INCLUDES_$(parent))

include $(MK)/footer.mk
