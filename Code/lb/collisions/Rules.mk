include $(MK)/header.mk

SRCS := Visitor.cc \
	Collision.cc \
        InletOutletCollision.cc \
        InletOutletWallCollision.cc \
        MidFluidCollision.cc \
        WallCollision.cc \
	StreamAndCollide.cc \
	PostStep.cc \

INCLUDES_$(d) := $(INCLUDES_$(parent))

include $(MK)/footer.mk
