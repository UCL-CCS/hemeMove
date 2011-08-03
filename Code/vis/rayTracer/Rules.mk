include $(MK)/header.mk

SRCS :=	BlockIterator.cc \
	ClusterBuilder.cc \
	RayTracer.cc \
	RectangularIterator.cc \
	SiteIterator.cc \

INCLUDES_$(d) := $(INCLUDES_$(parent))

include $(MK)/footer.mk
