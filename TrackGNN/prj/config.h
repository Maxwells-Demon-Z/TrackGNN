
# ifndef __CONFIG_H_ 
# define __CONFIG_H_ 

// Mode.
# define SIM
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
# define IDX_WID 12
# define DT0_INT_WID 8
# define DT1_INT_WID 6
# define DT_FRAC_WID 8
typedef ap_uint<1> BOL_TYPE;
typedef ap_int<IDX_WID> IDX_TYPE;
typedef ap_fixed<DT0_INT_WID+DT_FRAC_WID, DT0_INT_WID, AP_RND, AP_SAT> DT0_TYPE;
typedef ap_fixed<DT1_INT_WID+DT_FRAC_WID, DT1_INT_WID, AP_RND, AP_SAT> DT1_TYPE;

// Maximum numbers.
constexpr int GRAPH_NUM_MAX = 1024;
constexpr int NODE_NUM_MAX = 512;
constexpr int EDGE_NUM_MAX = 1024;
// Dimensions.
constexpr int NODE_FEA_DIM = 5;
constexpr int NODE_EMB_DIM = 8;
// Layers.
constexpr int NUM_OF_LAYERS = 1;
// Parallel parameters.
constexpr int NODE_UNI_PARA = 4;
constexpr int EDGE_UNI_PARA = 8;
constexpr int DIM_PARA = 2;
// Times of calculation.
constexpr int NODE_UNI_COUNT = ceildiv(NODE_NUM_MAX, NODE_UNI_PARA);
constexpr int EDGE_UNI_COUNT = ceildiv(EDGE_NUM_MAX, EDGE_UNI_PARA);
constexpr int DIM_COUNT = ceildiv(NODE_EMB_DIM, DIM_PARA);
// Input data size.
//   Maximum bit width:
//     64 bytes (512 bits).
//   Input node feature:
//     2 bytes (16 bits -> 16 bits) per data.
//     8 (5 -> 8) data per group.
//     64/(2*8) = 4 groups per cycle.
//     512/4 = 128 cycles totally.
//   Input adjacent list:
//     2 bytes (12 bits -> 16 bits) per data.
//     2 data per group.
//     64/(2*2) = 16 groups per cycle.
//     1024/16 = 64 cycles totally.
constexpr int INPUT_FEA_DIM_SIZE = 8;
constexpr int INPUT_ADJ_DIM_SIZE = 2;
constexpr int INPUT_FEA_UNI_SIZE = 4;
constexpr int INPUT_ADJ_UNI_SIZE = 16;
constexpr int INPUT_FEA_COUNT = 128;
constexpr int INPUT_ADJ_COUNT = 64;
// Output data size.
//   Maximum bit width:
//     64 bytes (512 bits).
//   Output result:
//     2 bytes (14 bits -> 16 bits) per data.
//     8 data per group.
//     1 group per cycle.
//     1024/8 = 128 cycles totally.
constexpr int OUTPUT_RESULT_SIZE = EDGE_UNI_PARA;
constexpr int OUTPUT_RESULT_COUNT = EDGE_UNI_COUNT;
// FIFO parameter.
constexpr int FIFO_DEPTH = 8;

// Activation functions.
//# define ACT_SIGM
# define ACT_TANH
# define NODE_MLP_LAYER 4
# define EDGE_MLP_LAYER 4
# define SIGM_SIZE 26
# define TANH_SIZE 29
# define SIGM_TAIL 6
# define TANH_TAIL 5

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

