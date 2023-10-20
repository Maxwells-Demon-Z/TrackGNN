
# Configuration file template.
TEMP = """
# ifndef __CONFIG_H_ 
# define __CONFIG_H_ 

// Mode.
//# define SIM
//# define EMU

// Includes.
# include <array>
# include "ap_fixed.h"
# include "ap_int.h"
# include "hls_stream.h"

// Ceiling division function.
template <typename T>
static constexpr T ceildiv(T dividend, T divisor) {
    # pragma HLS INLINE
    return (dividend + divisor - 1) / divisor; }

// Prime data types.
# define IDX_WID {TEMP_IDX_WID}
# define DT0_INT_WID {TEMP_DT0_INT_WID}
# define DT1_INT_WID {TEMP_DT1_INT_WID}
# define DT_FRAC_WID {TEMP_DT_FRAC_WID}
typedef ap_uint<1> BOL_TYPE;
typedef ap_int<IDX_WID> IDX_TYPE;
typedef ap_fixed<DT0_INT_WID+DT_FRAC_WID, DT0_INT_WID, AP_RND, AP_SAT> DT0_TYPE;
typedef ap_fixed<DT1_INT_WID+DT_FRAC_WID, DT1_INT_WID, AP_RND, AP_SAT> DT1_TYPE;

// Maximum numbers.
constexpr int GRAPH_NUM_MAX = 1024;
constexpr int NODE_NUM_MAX = {TEMP_NODE_NUM_MAX};
constexpr int EDGE_NUM_MAX = {TEMP_EDGE_NUM_MAX};
// Dimensions.
constexpr int NODE_FEA_DIM = {TEMP_NODE_FEA_DIM};
constexpr int NODE_EMB_DIM = {TEMP_NODE_EMB_DIM};
// Layers.
constexpr int NUM_OF_LAYERS = {TEMP_NUM_OF_LAYERS};
// Parallel parameters.
constexpr int NODE_UNI_PARA = {TEMP_NODE_UNI_PARA};
constexpr int EDGE_UNI_PARA = {TEMP_EDGE_UNI_PARA};
constexpr int DIM_PARA = {TEMP_DIM_PARA};
// Times of calculation.
constexpr int NODE_UNI_COUNT = ceildiv(NODE_NUM_MAX, NODE_UNI_PARA);
constexpr int EDGE_UNI_COUNT = ceildiv(EDGE_NUM_MAX, EDGE_UNI_PARA);
constexpr int DIM_COUNT = ceildiv(NODE_EMB_DIM, DIM_PARA);
// Input data size.
//   Maximum bit width:
//     64 bytes (512 bits).
//   Input node feature:
//     {TEMP_INPUT_0} bytes ({TEMP_INPUT_1} bits -> {TEMP_INPUT_2} bits) per data.
//     {TEMP_INPUT_3} ({TEMP_INPUT_4} -> {TEMP_INPUT_3}) data per group.
//     64/({TEMP_INPUT_0}*{TEMP_INPUT_3}) = {TEMP_INPUT_5} groups per cycle.
//     {TEMP_NODE_NUM_MAX}/{TEMP_INPUT_5} = {TEMP_INPUT_6} cycles totally.
//   Input adjacent list:
//     {TEMP_INPUT_10} bytes ({TEMP_INPUT_11} bits -> {TEMP_INPUT_12} bits) per data.
//     2 data per group.
//     64/({TEMP_INPUT_10}*2) = {TEMP_INPUT_13} groups per cycle.
//     {TEMP_EDGE_NUM_MAX}/{TEMP_INPUT_13} = {TEMP_INPUT_14} cycles totally.
constexpr int INPUT_FEA_DIM_SIZE = {TEMP_INPUT_3};
constexpr int INPUT_ADJ_DIM_SIZE = 2;
constexpr int INPUT_FEA_UNI_SIZE = {TEMP_INPUT_5};
constexpr int INPUT_ADJ_UNI_SIZE = {TEMP_INPUT_13};
constexpr int INPUT_FEA_COUNT = {TEMP_INPUT_6};
constexpr int INPUT_ADJ_COUNT = {TEMP_INPUT_14};
// Output data size.
//   Maximum bit width:
//     64 bytes (512 bits).
//   Output result:
//     {TEMP_OUTPUT_0} bytes ({TEMP_OUTPUT_1} bits -> {TEMP_OUTPUT_2} bits) per data.
//     {TEMP_EDGE_UNI_PARA} data per group.
//     1 group per cycle.
//     {TEMP_EDGE_NUM_MAX}/{TEMP_EDGE_UNI_PARA} = {TEMP_OUTPUT_3} cycles totally.
constexpr int OUTPUT_RESULT_SIZE = EDGE_UNI_PARA;
constexpr int OUTPUT_RESULT_COUNT = EDGE_UNI_COUNT;
// FIFO parameter.
constexpr int FIFO_DEPTH = {TEMP_FIFO_DEPTH};

// Activation functions.
//# define ACT_SIGM
//# define ACT_TANH
# define NODE_MLP_LAYER 4
# define EDGE_MLP_LAYER 4
# define SIGM_SIZE {TEMP_SIGM_SIZE}
# define TANH_SIZE {TEMP_TANH_SIZE}
# define SIGM_TAIL {TEMP_SIGM_TAIL}
# define TANH_TAIL {TEMP_TANH_TAIL}

// Secondary data types.
typedef std::array<DT0_TYPE, DIM_PARA> DT0_PARA;
typedef std::array<DT1_TYPE, DIM_PARA> DT1_PARA;
typedef std::array<DT0_TYPE, NODE_EMB_DIM> DT0_VECTOR;
typedef std::array<DT1_TYPE, NODE_EMB_DIM> DT1_VECTOR;
struct EDGE_TYPE{
	DT1_PARA node_emb_src;
	DT1_PARA node_emb_dst; };
typedef std::array<DT0_TYPE, INPUT_FEA_UNI_SIZE*INPUT_FEA_DIM_SIZE> INPUT_FEA_TYPE;
typedef std::array<IDX_TYPE, INPUT_ADJ_UNI_SIZE*INPUT_ADJ_DIM_SIZE> INPUT_ADJ_TYPE;
typedef std::array<DT1_TYPE, EDGE_UNI_PARA> OUTPUT_RESULT_TYPE;

// Kernel.
extern "C" {
void kernel_compute_graph(
		// Numbers.
		IDX_TYPE num_of_graphs,
		IDX_TYPE num_of_nodes[GRAPH_NUM_MAX],
		IDX_TYPE num_of_edges[GRAPH_NUM_MAX],
		// Input data.
		INPUT_FEA_TYPE node_feature[GRAPH_NUM_MAX][INPUT_FEA_COUNT],
		INPUT_ADJ_TYPE adj_list[GRAPH_NUM_MAX][INPUT_ADJ_COUNT],
		// Output data.
		OUTPUT_RESULT_TYPE result[GRAPH_NUM_MAX][OUTPUT_RESULT_COUNT]);
}

# endif 

"""

# Main function.
def generate_config(
        sim_or_emu,
        act_func,
        node_num_max,
        edge_num_max,
        dt_frac_wid,
        max_data_of_wandb,
        max_data_of_cpt,
        node_fea_dim,
        node_emb_dim,
        node_uni_para,
        edge_uni_para,
        dim_para,
        sigm_tail,
        sigm_size,
        tanh_tail,
        tanh_size,
        num_of_layers,
        fifo_depth):
    
    # Print the information.
    print('INFO: Start generating config file.')

    # Get the template.
    temp = TEMP
    
    # Simulate or emulate.
    if sim_or_emu:
        temp = temp.replace('//# define SIM', '# define SIM')
    else:
        temp = temp.replace('//# define EMU', '# define EMU')

    # Activation function.
    if act_func == 'sigmoid':
        temp = temp.replace('//# define ACT_SIGM', '# define ACT_SIGM')
    elif act_func == 'tanh':
        temp = temp.replace('//# define ACT_TANH', '# define ACT_TANH')
    else:
        raise Exception('Invalid activation function.')
    
    # Node and edge numbers and index width.
    num_of_nodes = 1
    num_of_edges = 1
    width_of_node = 1
    width_of_edge = 1
    while num_of_nodes < node_num_max:
        num_of_nodes *= 2 
        width_of_node += 1
    while num_of_edges < edge_num_max:
        num_of_edges *= 2
        width_of_edge += 1
    width_of_idx = max(width_of_node, width_of_edge) + 1
    temp = temp.replace('{TEMP_IDX_WID}', str(width_of_idx))
    temp = temp.replace('{TEMP_NODE_NUM_MAX}', str(num_of_nodes))
    temp = temp.replace('{TEMP_EDGE_NUM_MAX}', str(num_of_edges))
    
    # Integer widths.
    data_of_cpt = 1.0
    data_of_wandb = 1.0
    width_of_dt0_int = 1
    width_of_dt1_int = 1
    while data_of_cpt < max_data_of_cpt:
        data_of_cpt *= 2
        width_of_dt0_int += 1
    while data_of_wandb < max_data_of_wandb:
        data_of_wandb *= 2
        width_of_dt1_int += 1
    width_of_dt0_int += 2
    width_of_dt1_int += 2
    temp = temp.replace('{TEMP_DT0_INT_WID}', str(width_of_dt0_int))
    temp = temp.replace('{TEMP_DT1_INT_WID}', str(width_of_dt1_int))

    # Fraction width.
    temp = temp.replace('{TEMP_DT_FRAC_WID}', str(dt_frac_wid))

    # Feature and embedding dimensions.
    temp = temp.replace('{TEMP_NODE_FEA_DIM}', str(node_fea_dim))
    temp = temp.replace('{TEMP_NODE_EMB_DIM}', str(node_emb_dim))

    # Parallel parameters.
    temp = temp.replace('{TEMP_NODE_UNI_PARA}', str(node_uni_para))
    temp = temp.replace('{TEMP_EDGE_UNI_PARA}', str(edge_uni_para))
    temp = temp.replace('{TEMP_DIM_PARA}', str(dim_para))

    # Activation function parameters.
    temp = temp.replace('{TEMP_SIGM_SIZE}', str(sigm_size))
    temp = temp.replace('{TEMP_TANH_SIZE}', str(tanh_size))
    temp = temp.replace('{TEMP_SIGM_TAIL}', str(sigm_tail))
    temp = temp.replace('{TEMP_TANH_TAIL}', str(tanh_tail))

    # Number of layers.
    temp = temp.replace('{TEMP_NUM_OF_LAYERS}', str(num_of_layers))

    # FIFO depth.
    temp = temp.replace('{TEMP_FIFO_DEPTH}', str(fifo_depth))

    # Input data size.
    dt0_width = width_of_dt0_int + dt_frac_wid
    bytes = 1
    while 8*bytes < dt0_width:
        bytes *= 2
    data = 1
    while data < node_fea_dim:
        data *= 2
    group = 64 // (bytes*data)
    cycle = num_of_nodes // group
    temp = temp.replace('{TEMP_INPUT_0}', str(bytes))
    temp = temp.replace('{TEMP_INPUT_1}', str(dt0_width))
    temp = temp.replace('{TEMP_INPUT_2}', str(8*bytes))
    temp = temp.replace('{TEMP_INPUT_3}', str(data))
    temp = temp.replace('{TEMP_INPUT_4}', str(node_fea_dim))
    temp = temp.replace('{TEMP_INPUT_5}', str(group))
    temp = temp.replace('{TEMP_INPUT_6}', str(cycle))
    idx_width = width_of_idx
    bytes = 1
    while 8*bytes < idx_width:
        bytes *= 2
    group = 64 // (bytes*2)
    cycle = num_of_edges // group
    temp = temp.replace('{TEMP_INPUT_10}', str(bytes))
    temp = temp.replace('{TEMP_INPUT_11}', str(idx_width))
    temp = temp.replace('{TEMP_INPUT_12}', str(8*bytes))
    temp = temp.replace('{TEMP_INPUT_13}', str(group))
    temp = temp.replace('{TEMP_INPUT_14}', str(cycle))

    # Output data size.
    dt1_width = width_of_dt1_int + dt_frac_wid
    bytes = 1
    while 8*bytes < dt1_width:
        bytes *= 2
    cycle = num_of_edges // edge_uni_para
    temp = temp.replace('{TEMP_OUTPUT_0}', str(bytes))
    temp = temp.replace('{TEMP_OUTPUT_1}', str(dt1_width))
    temp = temp.replace('{TEMP_OUTPUT_2}', str(8*bytes))
    temp = temp.replace('{TEMP_OUTPUT_3}', str(cycle))

    # Write the config file.
    with open('./prj/config.h', 'w') as f:
        f.write(temp)
    
    # Print the information.
    print('INFO: Finish generating config file.')

