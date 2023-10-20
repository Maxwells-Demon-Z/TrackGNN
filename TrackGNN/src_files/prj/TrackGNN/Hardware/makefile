#
# this file was created by a computer. trust it.
#

# compiler tools
XILINX_VITIS ?= /tools/Xilinx/Vitis/2021.1
XILINX_XRT ?= /opt/xilinx/xrt

HOST_CXX ?= g++
RM = rm -f
RMDIR = rm -rf

VITIS_PLATFORM = xilinx_u280_xdma_201920_3
VITIS_PLATFORM_PATH = $(VITIS_PLATFORM)

# host compiler global settings
CXXFLAGS += -std=c++1y -DVITIS_PLATFORM=$(VITIS_PLATFORM) -D__USE_XOPEN2K8 -I$(XILINX_XRT)/include/ -I/tools/Xilinx/Vitis_HLS/2021.1/include/ -I"../../TrackGNN_kernels/src" -O2 -g -Wall -c -fmessage-length=0
LDFLAGS += -lxilinxopencl -lpthread -lrt -lstdc++ -L$(XILINX_XRT)/lib/ -Wl,-rpath-link,$(XILINX_XRT)/lib

#
# host files
#

HOST_OBJECTS += src/host.o
HOST_OBJECTS += src/xcl2.o

HOST_EXE = TrackGNN

BUILD_SUBDIRS += src/

#
# primary build targets
#

.PHONY: all clean
all:  $(HOST_EXE)

clean:
	-$(RM) $(HOST_EXE) $(HOST_OBJECTS)

.PHONY: incremental
incremental: all


nothing:

#
# host rules
#

src/host.o: ../src/host.cpp /usr/include/c++/7/fstream /usr/include/c++/7/vector
	-@mkdir -p $(@D)
	$(HOST_CXX) $(CXXFLAGS) -o "$@" "$<"

src/xcl2.o: ../src/xcl2.cpp ../src/xcl2.hpp
	-@mkdir -p $(@D)
	$(HOST_CXX) $(CXXFLAGS) -o "$@" "$<"

$(HOST_EXE): $(HOST_OBJECTS)
	$(HOST_CXX) -o "$@" $(+) $(LDFLAGS)

