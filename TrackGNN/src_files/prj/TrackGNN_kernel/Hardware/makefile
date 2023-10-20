#
# this file was created by a computer. trust it.
#

# compiler tools
XILINX_VITIS ?= /tools/Xilinx/Vitis/2021.1
XILINX_XRT ?= /opt/xilinx/xrt
XILINX_VIVADO ?= /tools/Xilinx/Vivado/2021.1
XILINX_HLS ?= $(XILINX_VITIS)/Vivado_HLS

VPP ?= ${XILINX_VITIS}/bin/v++
RM = rm -f
RMDIR = rm -rf

VITIS_PLATFORM = xilinx_u280_xdma_201920_3
VITIS_PLATFORM_PATH = $(VITIS_PLATFORM)

TARGET = hw
VPP_OPTS = --target $(TARGET)

#
# Hw kernel files
#


BUILD_SUBDIRS += build
HW_KERNEL_OBJS += build/kernel_compute_graph.xo
SRC_KERNEL_OBJS += build/kernel_compute_graph.xo

ALL_MESSAGE_FILES = $(subst .xo,.o,.mdb,$(SRC_KERNEL_OBJS))

#
# primary build targets
#

.PHONY: all clean
all: $(HW_KERNEL_OBJS)


clean:
	-$(RM) $(SRC_KERNEL_OBJS) $(ALL_MESSAGE_FILES) 
	-$(RMDIR) $(BUILD_SUBDIRS)
	-$(RMDIR) .Xil

.PHONY: incremental
incremental: all


nothing:

#
# Hw Kernels compile;
#

build/kernel_compute_graph.xo: ../src/kernel.cpp ../src/config.h ../src/data_files/input_nw_weight.txt ../src/data_files/input_nw_bias.txt ../src/data_files/node_nw_weight_0_0.txt ../src/data_files/node_nw_weight_0_1.txt ../src/data_files/node_nw_weight_0_2.txt ../src/data_files/node_nw_weight_1.txt ../src/data_files/node_nw_weight_2.txt ../src/data_files/node_nw_weight_3.txt ../src/data_files/node_nw_bias_0.txt ../src/data_files/node_nw_bias_1.txt ../src/data_files/node_nw_bias_2.txt ../src/data_files/node_nw_bias_3.txt ../src/data_files/edge_nw_weight_0_0.txt ../src/data_files/edge_nw_weight_0_1.txt ../src/data_files/edge_nw_weight_1.txt ../src/data_files/edge_nw_weight_2.txt ../src/data_files/edge_nw_weight_3.txt ../src/data_files/edge_nw_bias_0.txt ../src/data_files/edge_nw_bias_1.txt ../src/data_files/edge_nw_bias_2.txt ../src/data_files/edge_nw_bias_3.txt ../src/data_files/activation_sigm_result_rom.txt ../src/data_files/activation_sigm_node_nw_rom.txt ../src/data_files/activation_sigm_edge_nw_rom.txt ../src/data_files/activation_tanh_node_nw_rom.txt ../src/data_files/activation_tanh_edge_nw_rom.txt kernel_compute_graph-compile.cfg
	-@mkdir -p $(@D)
	-@$(RM) $@
	$(VPP) $(VPP_OPTS) --compile -I"$(<D)" --config kernel_compute_graph-compile.cfg -o"$@" "$<" 

