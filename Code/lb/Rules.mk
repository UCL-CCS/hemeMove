include $(MK)/header.mk

TARGETS := libHemeLbLb.$(LIBEXT)

SUBDIRS := collisions

$(TARGETS)_DEPS = lb.o \
                  io.o \
                  StabilityTester.o \
                  EntropyTester.o \
                  SimulationState.o \
                  BoundaryComms.o \
                  BoundaryValues.o \
		  $(foreach sd,$(SUBDIRS_$(d)),$(call subtree_tgts,$(sd)))

INCLUDES_$(d) := $(INCLUDES_$(parent))

include $(MK)/footer.mk
