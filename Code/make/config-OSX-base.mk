include $(MK)/config-default.mk

HEMELB_CFG_ON_OSX := true
HEMELB_CFG_ON_BSD := true

HEMELB_DEFS += HEMELB_CFG_ON_BSD \
	HEMELB_CFG_ON_OSX

HEMELB_CXXFLAGS += -pthread

