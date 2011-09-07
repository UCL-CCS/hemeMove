CXX := mpic++
CC := mpic++

EXE := hemelb
UNITTESTS := unitTests

HEMELB_DEBUG_LEVEL := 0
HEMELB_STEERING_LIB := $(or $(HEMELB_STEERING_LIB),basic)

PMETIS_INCLUDE_DIR := $(TOP)/parmetis/include $(TOP)/parmetis/metis/include
PMETIS_LIBRARY_DIR := $(TOP)/parmetis/build/Linux-x86_64/libmetis/ $(TOP)/parmetis/build/Linux-x86_64/libparmetis/
HEMELB_CFLAGS :=
HEMELB_CXXFLAGS := -g -pedantic -Wall -Wextra -Wno-unused
HEMELB_DEFS := TIXML_USE_STL
HEMELB_INCLUDEPATHS = $(TOP) $(PMETIS_INCLUDE_DIR)
HEMELB_LIBPATHS = $(PMETIS_LIBRARY_DIR)

$(EXE)_LIBS = -lparmetis -lmetis
$(UNITTESTS)_LIBS = -lparmetis -lmetis -lcppunit
$(UNITTESTS)_CXXFLAGS := $(EXE)_CXXFLAGS
$(UNITTESTS)_INCLUDEPATHS := $(EXE)_INCLUDEPATHS
$(UNITTESTS)_LIBPATHS := $(EXE)_LIBPATHS
