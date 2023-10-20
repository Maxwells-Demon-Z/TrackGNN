
// **************************** Readme **************************** //

/*

	The project is to design an accelerator to process data loaded
	from colliders, which can determine the authenticity of tracks.

	The design of accelerator bases on the article of FlowGNN, and
	refers to the open-source code - https://arxiv.org/abs/2204.13103
	and https://github.com/sharc-lab/FlowGNN.

	The algorithm of processing bases on the article of sPhenix, and
	refers to the open-source code - https://arxiv.org/abs/1810.06111
	and https://bitbucket.org/dtyu/trigger-detection-pipeline.

	If there is any issue or question, please contact Hanqing Zhang
	by the email hanqing.zhang@zju.edu.cn.

*/



// **************************** Includes **************************** //

# include "config.h"



// **************************** Global Variables **************************** //

// Input data.
DT0_TYPE g_node_feature[NODE_NUM_MAX][NODE_FEA_DIM];
IDX_TYPE g_adj_list_adapter[EDGE_UNI_PARA][EDGE_NUM_MAX][2];
IDX_TYPE g_adj_list_scatter[EDGE_UNI_PARA][EDGE_NUM_MAX][2];
IDX_TYPE g_determine[NODE_NUM_MAX];

// Ping-pong data.
DT0_TYPE g_msg_src_ping[2][EDGE_UNI_PARA][NODE_NUM_MAX][NODE_EMB_DIM];
DT0_TYPE g_msg_src_pong[2][EDGE_UNI_PARA][NODE_NUM_MAX][NODE_EMB_DIM];
DT0_TYPE g_msg_dst_ping[2][EDGE_UNI_PARA][NODE_NUM_MAX][NODE_EMB_DIM];
DT0_TYPE g_msg_dst_pong[2][EDGE_UNI_PARA][NODE_NUM_MAX][NODE_EMB_DIM];
IDX_TYPE g_msg_src_ping_flag[EDGE_UNI_PARA][NODE_NUM_MAX];
IDX_TYPE g_msg_src_pong_flag[EDGE_UNI_PARA][NODE_NUM_MAX];
IDX_TYPE g_msg_dst_ping_flag[EDGE_UNI_PARA][NODE_NUM_MAX];
IDX_TYPE g_msg_dst_pong_flag[EDGE_UNI_PARA][NODE_NUM_MAX];
DT1_TYPE g_node_emb_ping[NODE_NUM_MAX][NODE_EMB_DIM];
DT1_TYPE g_node_emb_pong[NODE_NUM_MAX][NODE_EMB_DIM];



// **************************** Weights and Biases **************************** //

// Input network.
DT1_TYPE g_input_nw_weight[NODE_EMB_DIM][NODE_FEA_DIM] = {
		# include "data_files/input_nw_weight.txt"
		};
DT1_TYPE g_input_nw_bias[NODE_EMB_DIM] = {
		# include "data_files/input_nw_bias.txt"
		};
// Node network.
DT1_TYPE g_node_nw_weight_0_0[NODE_EMB_DIM][NODE_EMB_DIM] = {
		# include "data_files/node_nw_weight_0_0.txt"
		};
DT1_TYPE g_node_nw_weight_0_1[NODE_EMB_DIM][NODE_EMB_DIM] = {
		# include "data_files/node_nw_weight_0_1.txt"
		};
DT1_TYPE g_node_nw_weight_0_2[NODE_EMB_DIM][NODE_EMB_DIM] = {
		# include "data_files/node_nw_weight_0_2.txt"
		};
DT1_TYPE g_node_nw_weight_1[NODE_EMB_DIM][NODE_EMB_DIM] = {
		# include "data_files/node_nw_weight_1.txt"
		};
DT1_TYPE g_node_nw_weight_2[NODE_EMB_DIM][NODE_EMB_DIM] = {
		# include "data_files/node_nw_weight_2.txt"
		};
DT1_TYPE g_node_nw_weight_3[NODE_EMB_DIM][NODE_EMB_DIM] = {
		# include "data_files/node_nw_weight_3.txt"
		};
DT1_TYPE g_node_nw_bias_0[NODE_EMB_DIM] = {
		# include "data_files/node_nw_bias_0.txt"
		};
DT1_TYPE g_node_nw_bias_1[NODE_EMB_DIM] = {
		# include "data_files/node_nw_bias_1.txt"
		};
DT1_TYPE g_node_nw_bias_2[NODE_EMB_DIM] = {
		# include "data_files/node_nw_bias_2.txt"
		};
DT1_TYPE g_node_nw_bias_3[NODE_EMB_DIM] = {
		# include "data_files/node_nw_bias_3.txt"
		};
// Edge network.
DT1_TYPE g_edge_nw_weight_0_0[NODE_EMB_DIM][NODE_EMB_DIM] = {
		# include "data_files/edge_nw_weight_0_0.txt"
		};
DT1_TYPE g_edge_nw_weight_0_1[NODE_EMB_DIM][NODE_EMB_DIM] = {
		# include "data_files/edge_nw_weight_0_1.txt"
		};
DT1_TYPE g_edge_nw_weight_1[NODE_EMB_DIM][NODE_EMB_DIM] = {
		# include "data_files/edge_nw_weight_1.txt"
		};
DT1_TYPE g_edge_nw_weight_2[NODE_EMB_DIM][NODE_EMB_DIM] = {
		# include "data_files/edge_nw_weight_2.txt"
		};
DT1_TYPE g_edge_nw_weight_3[NODE_EMB_DIM] = {
		# include "data_files/edge_nw_weight_3.txt"
		};
DT1_TYPE g_edge_nw_bias_0[NODE_EMB_DIM] = {
		# include "data_files/edge_nw_bias_0.txt"
		};
DT1_TYPE g_edge_nw_bias_1[NODE_EMB_DIM] = {
		# include "data_files/edge_nw_bias_1.txt"
		};
DT1_TYPE g_edge_nw_bias_2[NODE_EMB_DIM] = {
		# include "data_files/edge_nw_bias_2.txt"
		};
DT1_TYPE g_edge_nw_bias_3 =
		# include "data_files/edge_nw_bias_3.txt"
		;



// **************************** Look-Up Tables **************************** //

// Look-up tables in activation functions.
DT1_TYPE activation_sigm_result_rom[EDGE_UNI_PARA][SIGM_SIZE] = {
		# include "data_files/activation_sigm_result_rom.txt"
		};
DT1_TYPE activation_sigm_node_nw_rom[NODE_MLP_LAYER][NODE_UNI_PARA][DIM_PARA][SIGM_SIZE] = {
		# include "data_files/activation_sigm_node_nw_rom.txt"
		};
DT1_TYPE activation_sigm_edge_nw_rom[EDGE_MLP_LAYER][EDGE_UNI_PARA][DIM_PARA][SIGM_SIZE] = {
		# include "data_files/activation_sigm_edge_nw_rom.txt"
		};
DT1_TYPE activation_tanh_node_nw_rom[NODE_MLP_LAYER][NODE_UNI_PARA][DIM_PARA][TANH_SIZE] = {
		# include "data_files/activation_tanh_node_nw_rom.txt"
		};
DT1_TYPE activation_tanh_edge_nw_rom[EDGE_MLP_LAYER][EDGE_UNI_PARA][DIM_PARA][TANH_SIZE] = {
		# include "data_files/activation_tanh_edge_nw_rom.txt"
		};



// **************************** Load **************************** //

/*

	This file is to load input data, including node features and
	adjacent list, and then rename indexes of nodes in order to make
	the edge-processing network work as soon as possible.

	Version: 0.1
		Author: Hanqing Zhang
		Email : hanqing.zhang@zju.edu.cn
		Date  : June-30-2023
		Commit: Add function "load".

	Version: 0.2
		Author: Hanqing Zhang
		Email : hanqing.zhang@zju.edu.cn
		Date  : July-4-2023
		Commit: Add "merge_sort" function to sort the adjacent list.

	Version: 0.3
		Author: Hanqing Zhang
		Email : hanqing.zhang@zju.edu.cn
		Date  : July-7-2023
		Commit: Divide the array "adj_list" into two arrays because
				modules "adapter" and "msg_scatter" need this array
				simultaneously.

	Version: 0.4
		Author: Hanqing Zhang
		Email : hanqing.zhang@zju.edu.cn
		Date  : July-12-2023
		Commit: Change the "merge_sort" function because it can work
				by O(n) instead of O(nlogn).

	Version: 0.5
		Author: Hanqing Zhang
		Email : hanqing.zhang@zju.edu.cn
		Date  : July-20-2023
		Commit: Modify the "merge_sort" function to fixed the long-
				time-synthesis problem.

	Version: 0.6
		Author: Hanqing Zhang
		Email : hanqing.zhang@zju.edu.cn
		Date  : August-5-2023
		Commit: Delete "sort" function and add "rename" function,
				because the unit can work without sorting.

	Version: 0.7
		Author: Hanqing Zhang
		Email : hanqing.zhang@zju.edu.cn
		Date  : September-2-2023
		Commit: Merge "rename" and "load_graph" functions to reduce
				invoking time.
	
	Version: 0.8
		Author: Hanqing Zhang
		Email : hanqing.zhang@zju.edu.cn
		Date  : September-29-2023
		Commit: Update the number of nodes to omit isolated nodes.

	Function:
		void load_graph(
			IDX_TYPE num_of_nodes,
			IDX_TYPE num_of_edges,
			INPUT_FEA_TYPE node_feature[INPUT_FEA_COUNT],
			INPUT_ADJ_TYPE adj_list[INPUT_ADJ_COUNT]);
		Input arguments:
			num_of_nodes: the number of nodes from s_axilite interface.
			num_of_edges: the number of edges from s_axilite interface.
			node_feature: data of node features from m_axi interface.
			adj_list: data of adjacent list from m_axi interface.

*/

static IDX_TYPE load_graph(
		IDX_TYPE num_of_nodes,
		IDX_TYPE num_of_edges,
		INPUT_FEA_TYPE node_feature[INPUT_FEA_COUNT],
		INPUT_ADJ_TYPE adj_list[INPUT_ADJ_COUNT])
{
	# pragma HLS INLINE off

	// Declarations of arrays loading input data.
	DT0_TYPE node_feature_cache[NODE_NUM_MAX][NODE_FEA_DIM];
	IDX_TYPE adj_list_cache[EDGE_NUM_MAX][2];
	# pragma HLS ARRAY_PARTITION variable=node_feature_cache dim=1 type=cyclic factor=INPUT_FEA_UNI_SIZE
	# pragma HLS ARRAY_PARTITION variable=node_feature_cache dim=2 type=complete
	# pragma HLS ARRAY_PARTITION variable=adj_list_cache dim=1 type=cyclic factor=INPUT_ADJ_UNI_SIZE
	# pragma HLS ARRAY_PARTITION variable=adj_list_cache dim=2 type=complete

	// Declarations of arrays used to rename indexes.
	IDX_TYPE map_list[NODE_NUM_MAX];
	IDX_TYPE edge_idx_base;
	IDX_TYPE node_idx_max = (IDX_TYPE)0;
	IDX_TYPE node_idx_now = (IDX_TYPE)0;
	IDX_TYPE node_idx_last = (IDX_TYPE)-1;
	IDX_TYPE node_idx_next = (IDX_TYPE)-1;

	// Load node features.
	IDX_TYPE iter_fea_num = ceildiv(num_of_nodes, (IDX_TYPE)INPUT_FEA_UNI_SIZE);
	for (	IDX_TYPE iter_fea_cnt=(IDX_TYPE)0;
			iter_fea_cnt<iter_fea_num;
			iter_fea_cnt+=(IDX_TYPE)1)
	{
		# pragma HLS LOOP_TRIPCOUNT min=INPUT_FEA_COUNT max=INPUT_FEA_COUNT avg=INPUT_FEA_COUNT
		# pragma HLS PIPELINE II=1

		INPUT_FEA_TYPE cache = node_feature[iter_fea_cnt];

		for (IDX_TYPE uni=0; uni<INPUT_FEA_UNI_SIZE; uni++) {
			# pragma HLS UNROLL
			IDX_TYPE dst_idx = iter_fea_cnt * (IDX_TYPE)INPUT_FEA_UNI_SIZE + uni;
			for (IDX_TYPE dim=0; dim<NODE_FEA_DIM; dim++) {
				# pragma HLS UNROLL
				IDX_TYPE src_idx = uni * (IDX_TYPE)INPUT_FEA_DIM_SIZE + dim;
				node_feature_cache[dst_idx][dim] = cache[src_idx]; }}
	}

	// Load adjacent list.
	IDX_TYPE iter_adj_num = ceildiv(num_of_edges, (IDX_TYPE)INPUT_ADJ_UNI_SIZE);
	for (	IDX_TYPE iter_adj_cnt=(IDX_TYPE)0;
			iter_adj_cnt<iter_adj_num;
			iter_adj_cnt+=(IDX_TYPE)1)
	{
		# pragma HLS LOOP_TRIPCOUNT min=INPUT_ADJ_COUNT max=INPUT_ADJ_COUNT avg=INPUT_ADJ_COUNT
		# pragma HLS PIPELINE II=1

		INPUT_ADJ_TYPE cache = adj_list[iter_adj_cnt];

		for (IDX_TYPE uni=0; uni<INPUT_ADJ_UNI_SIZE; uni++) {
			# pragma HLS UNROLL
			IDX_TYPE dst_idx = iter_adj_cnt * (IDX_TYPE)INPUT_ADJ_UNI_SIZE + uni;
			for (IDX_TYPE dim=0; dim<2; dim++) {
				# pragma HLS UNROLL
				IDX_TYPE src_idx = uni * (IDX_TYPE)INPUT_ADJ_DIM_SIZE + dim;
				adj_list_cache[dst_idx][dim] = cache[src_idx]; }}
	}

	// Initialize arrays.
	for (	IDX_TYPE node_idx=(IDX_TYPE)0;
			node_idx<num_of_nodes;
			node_idx+=(IDX_TYPE)1)
	{
		# pragma HLS LOOP_TRIPCOUNT min=NODE_NUM_MAX max=NODE_NUM_MAX avg=NODE_NUM_MAX
		# pragma HLS PIPELINE II=1
		map_list[node_idx] = (IDX_TYPE)-1;
	}
	for (	IDX_TYPE node_idx_base=(IDX_TYPE)0;
			node_idx_base<num_of_nodes;
			node_idx_base+=(IDX_TYPE)NODE_UNI_PARA)
	{
		# pragma HLS LOOP_TRIPCOUNT min=NODE_UNI_COUNT max=NODE_UNI_COUNT avg=NODE_UNI_COUNT
		# pragma HLS PIPELINE II=1
		g_determine[node_idx_base] = (IDX_TYPE)-1;
	}
	for (	IDX_TYPE node_idx_base=(IDX_TYPE)0;
			node_idx_base<num_of_nodes;
			node_idx_base+=(IDX_TYPE)NODE_UNI_PARA)
	{
		# pragma HLS LOOP_TRIPCOUNT min=NODE_UNI_COUNT max=NODE_UNI_COUNT avg=NODE_UNI_COUNT
		# pragma HLS PIPELINE II=1

		for (IDX_TYPE node_idx_offset=0; node_idx_offset<NODE_UNI_PARA; node_idx_offset++) {
			# pragma HLS UNROLL
			IDX_TYPE node_idx = node_idx_base + node_idx_offset;
			BOL_TYPE node_idx_overflow = (node_idx >= num_of_nodes);
			if (!node_idx_overflow)
			for (IDX_TYPE uni=0; uni<EDGE_UNI_PARA; uni++) {
				# pragma HLS UNROLL
				g_msg_src_ping_flag[uni][node_idx] = (IDX_TYPE)-1;
				g_msg_src_pong_flag[uni][node_idx] = (IDX_TYPE)-1;
				g_msg_dst_ping_flag[uni][node_idx] = (IDX_TYPE)-1;
				g_msg_dst_pong_flag[uni][node_idx] = (IDX_TYPE)-1; }}
	}

	// Iterate to rename each node in the adjacent list.
	IDX_TYPE iter_cnt = (IDX_TYPE)0;
	IDX_TYPE iter_num = num_of_edges * (IDX_TYPE)2 - (IDX_TYPE)1;
	while (true)
	{
		# pragma HLS LOOP_TRIPCOUNT min=(EDGE_NUM_MAX*2) max=(EDGE_NUM_MAX*2) avg=(EDGE_NUM_MAX*2)
		# pragma HLS PIPELINE II=2

		// Get indexes.
		IDX_TYPE edge_idx = iter_cnt / (IDX_TYPE)2;
		IDX_TYPE io = iter_cnt % (IDX_TYPE)2;

		// Flags of the first and the last loops.
		BOL_TYPE first =
				(edge_idx % (IDX_TYPE)EDGE_UNI_PARA == (IDX_TYPE)0) &&
				(io == (IDX_TYPE)0);
		BOL_TYPE last =
				((edge_idx % (IDX_TYPE)EDGE_UNI_PARA == (IDX_TYPE)EDGE_UNI_PARA - (IDX_TYPE)1) ||
				(edge_idx == num_of_edges - (IDX_TYPE)1)) &&
				(io == (IDX_TYPE)1);

		// If it is the first loop.
		// Initialize the maximum node index.
		// Save the base edge index.
		if (first) {
			node_idx_max = (IDX_TYPE)-1;
			edge_idx_base = edge_idx; }

		// Get the old node index in the adjacent list.
		// Get the new index.
		// Whether the node has been chosen.
		// If true, use the new index.
		// If false, rename it.
		// Get the max index in an edge group.
		IDX_TYPE node_idx_old = adj_list_cache[edge_idx][io];
		IDX_TYPE node_idx_new = map_list[node_idx_old];
		BOL_TYPE occupied = (node_idx_new != (IDX_TYPE)-1);
		IDX_TYPE node_idx = occupied? node_idx_new: node_idx_now;
		node_idx_max = (node_idx_max > node_idx)? node_idx_max: node_idx;
		node_idx_now += occupied? (IDX_TYPE)0: (IDX_TYPE)1;
		map_list[node_idx_old] = node_idx;

		// New adjacent list.
		// New node feature.
		for (IDX_TYPE uni=0; uni<EDGE_UNI_PARA; uni++) {
			# pragma HLS UNROLL
			g_adj_list_adapter[uni][edge_idx][io] = node_idx;
			g_adj_list_scatter[uni][edge_idx][io] = node_idx; }
		for (IDX_TYPE dim=0; dim<NODE_FEA_DIM; dim++) {
			# pragma HLS UNROLL
			g_node_feature[node_idx][dim] = node_feature_cache[node_idx_old][dim]; }

		// If it is the last loop.
		// Update the data in the determination array.
		if (last) {
			// Get the base node index.
			IDX_TYPE node_idx_base = node_idx_max - node_idx_max % (IDX_TYPE)NODE_UNI_PARA;
			// If the index is not greater than the last one.
			if (node_idx_base <= node_idx_last) {
				if (node_idx_next < num_of_nodes) {
				g_determine[node_idx_next] = edge_idx_base;
				node_idx_last = node_idx_next;
				node_idx_next = node_idx_next + (IDX_TYPE)NODE_UNI_PARA; }}
			// If the index is greater than the last one.
			else {
				g_determine[node_idx_base] = edge_idx_base;
				node_idx_last = node_idx_base;
				node_idx_next = node_idx_base + (IDX_TYPE)NODE_UNI_PARA; }}

		// Update iteration variable.
		if (iter_cnt < iter_num) {
			iter_cnt += (IDX_TYPE)1; }
		else break;
	}

	// New number of nodes.
	return node_idx_now;
}



// **************************** Activation Functions **************************** //

/*

	This part is to implement activation functions, including Sigmoid 
	and Tanh, which are implemented by the method of look-up table.

	Version: 0.1
		Author: Hanqing Zhang
		Email : hanqing.zhang@zju.edu.cn
		Date  : July-2-2023
		Commit: Add ROMs declarations and functions.

	Version: 0.2
		Author: Hanqing Zhang
		Email : hanqing.zhang@zju.edu.cn
		Date  : July-3-2023
		Commit: Fix the bug of index overflow.

	Function:
		DT1_TYPE act_sigm_result(
			IDX_TYPE uni_idx,
			DT0_TYPE in);
		Input arguments:
			uni_idx: index of the edge-processing unit.
			in: input data.
		Return value type:
			DT1_TYPE.

	Function:
		DT1_TYPE act_sigm_node_nw(
			IDX_TYPE layer_idx,
			IDX_TYPE uni_idx,
			IDX_TYPE dim_idx,
			DT0_TYPE in);
		Input arguments:
			layer_idx: index of the MLP layer.
			uni_idx: index of the node-processing unit.
			dim_idx: index of the dimension.
			in: input data.
		Return value type:
			DT1_TYPE.

	Function:
		DT1_TYPE act_sigm_edge_nw(
			IDX_TYPE layer_idx,
			IDX_TYPE uni_idx,
			IDX_TYPE dim_idx,
			DT0_TYPE in);
		Input arguments:
			layer_idx: index of the MLP layer.
			uni_idx: index of the edge-processing unit.
			dim_idx: index of the dimension.
			in: input data.
		Return value type:
			DT1_TYPE.

	Function:
		DT1_TYPE act_tanh_node_nw(
			IDX_TYPE layer_idx,
			IDX_TYPE uni_idx,
			IDX_TYPE dim_idx,
			DT0_TYPE in);
		Input arguments:
			layer_idx: index of the MLP layer.
			uni_idx: index of the node-processing unit.
			dim_idx: index of the dimension.
			in: input data.
		Return value type:
			DT1_TYPE.

	Function:
		DT1_TYPE act_tanh_edge_nw(
			IDX_TYPE layer_idx,
			IDX_TYPE uni_idx,
			IDX_TYPE dim_idx,
			DT0_TYPE in);
		Input arguments:
			layer_idx: index of the MLP layer.
			uni_idx: index of the edge-processing unit.
			dim_idx: index of the dimension.
			in: input data.
		Return value type:
			DT1_TYPE.

*/

static DT1_TYPE act_sigm_result(
		IDX_TYPE uni_idx,
		DT0_TYPE in)
{
	# pragma HLS INLINE
	# pragma HLS ARRAY_PARTITION variable=activation_sigm_result_rom dim=1 type=complete

	// Whether the input data is negative.
	// If it is negative, convert it into positive.
	BOL_TYPE sign;
	DT0_TYPE pos_in;
	# pragma HLS AGGREGATE variable=pos_in
	if (in > (DT0_TYPE)0.0) {
		sign = false;
		pos_in = in; }
	else {
		sign = true;
		pos_in = -in; }

	// Divide bits of the input data into two parts.
	// The front part is the index of the look-up table.
	// The back part is used to interpolate.
	typedef ap_uint<DT0_INT_WID+DT_FRAC_WID-SIGM_TAIL-1> BASE_TYPE;
	typedef ap_uint<SIGM_TAIL+1> OFFSET_TYPE;
	BASE_TYPE base = 0;
	OFFSET_TYPE offset_0 = 0;
	OFFSET_TYPE offset_1 = 0;
	# pragma HLS AGGREGATE variable=base
	# pragma HLS AGGREGATE variable=offset_0
	# pragma HLS AGGREGATE variable=offset_1
	for (IDX_TYPE i=0; i<DT0_INT_WID+DT_FRAC_WID-SIGM_TAIL-1; i++)
		# pragma HLS UNROLL
		base[i] = pos_in[i+SIGM_TAIL];
	for (IDX_TYPE i=0; i<SIGM_TAIL; i++)
		# pragma HLS UNROLL
		offset_1[i] = pos_in[i];
	offset_0 = (1 << SIGM_TAIL) - offset_1;

	// If the index is overflowed.
	// Output the minimum or the maximum data.
	if (base >= SIGM_SIZE - 1)
	{
		if (sign)
			return (DT1_TYPE)0.0;
		else
			return (DT1_TYPE)1.0;
	}

	// If the index is not overflowed.
	// Calculate by linear interpolation algorithm.
	else
	{
		// Y(X) = (Y0 * X1 + Y1 * X0) / (X0 + X1).
		// Get the two data according to the index.
		typedef ap_fixed<SIGM_TAIL+DT_FRAC_WID+2, SIGM_TAIL+2> EXP_TYPE;
		EXP_TYPE y_0 = activation_sigm_result_rom[uni_idx][base    ];
		EXP_TYPE y_1 = activation_sigm_result_rom[uni_idx][base + 1];

		// Calculate (Y0 * X1) and (Y1 * X0).
		// Data types of X0 and X1 are integer.
		// Data types of Y0 and Y1 are fixed-point.
		// Multiply Y0 and Y1 by each bit of X0 and X1.
		EXP_TYPE sub_0[SIGM_TAIL+1];
		EXP_TYPE sub_1[SIGM_TAIL+1];
		# pragma HLS ARRAY_PARTITION variable=sub_0 dim=1 type=complete
		# pragma HLS ARRAY_PARTITION variable=sub_1 dim=1 type=complete
		for (IDX_TYPE i=0; i<SIGM_TAIL+1; i++)
			# pragma HLS UNROLL
			if (offset_0[i] == 1) sub_0[i] = (y_0 << i);
			else sub_0[i] = 0.0;
		for (IDX_TYPE i=0; i<SIGM_TAIL+1; i++)
			# pragma HLS UNROLL
			if (offset_1[i] == 1) sub_1[i] = (y_1 << i);
			else sub_1[i] = 0.0;

		// Calculate (Y0 * X1 + Y1 * X0) / (X0 + X1).
		// Sum fixed-point numbers together.
		EXP_TYPE sum = 0.0;
		for (IDX_TYPE i=0; i<SIGM_TAIL+1; i++) {
			# pragma HLS UNROLL
			sum += sub_0[i];
			sum += sub_1[i]; }
		DT1_TYPE result = (sum >> SIGM_TAIL);

		// Output data.
		if (sign)
			return ((DT1_TYPE)1.0 - result);
		else
			return (                result);
	}
}

static DT1_TYPE act_sigm_node_nw(
		IDX_TYPE layer_idx,
		IDX_TYPE uni_idx,
		IDX_TYPE dim_idx,
		DT0_TYPE in)
{
	# pragma HLS INLINE
	# pragma HLS ARRAY_PARTITION variable=activation_sigm_node_nw_rom dim=1 type=complete
	# pragma HLS ARRAY_PARTITION variable=activation_sigm_node_nw_rom dim=2 type=complete
	# pragma HLS ARRAY_PARTITION variable=activation_sigm_node_nw_rom dim=3 type=complete

	// Whether the input data is negative.
	// If it is negative, convert it into positive.
	BOL_TYPE sign;
	DT0_TYPE pos_in;
	# pragma HLS AGGREGATE variable=pos_in
	if (in > (DT0_TYPE)0.0) {
		sign = false;
		pos_in = in; }
	else {
		sign = true;
		pos_in = -in; }

	// Divide bits of the input data into two parts.
	// The front part is the index of the look-up table.
	// The back part is used to interpolate.
	typedef ap_uint<DT0_INT_WID+DT_FRAC_WID-SIGM_TAIL-1> BASE_TYPE;
	typedef ap_uint<SIGM_TAIL+1> OFFSET_TYPE;
	BASE_TYPE base = 0;
	OFFSET_TYPE offset_0 = 0;
	OFFSET_TYPE offset_1 = 0;
	# pragma HLS AGGREGATE variable=base
	# pragma HLS AGGREGATE variable=offset_0
	# pragma HLS AGGREGATE variable=offset_1
	for (IDX_TYPE i=0; i<DT0_INT_WID+DT_FRAC_WID-SIGM_TAIL-1; i++)
		# pragma HLS UNROLL
		base[i] = pos_in[i+SIGM_TAIL];
	for (IDX_TYPE i=0; i<SIGM_TAIL; i++)
		# pragma HLS UNROLL
		offset_1[i] = pos_in[i];
	offset_0 = (1 << SIGM_TAIL) - offset_1;

	// If the index is overflowed.
	// Output the minimum or the maximum data.
	if (base >= SIGM_SIZE - 1)
	{
		if (sign)
			return (DT1_TYPE)0.0;
		else
			return (DT1_TYPE)1.0;
	}

	// If the index is not overflowed.
	// Calculate by linear interpolation algorithm.
	else
	{
		// Y(X) = (Y0 * X1 + Y1 * X0) / (X0 + X1).
		// Get the two data according to the index.
		typedef ap_fixed<SIGM_TAIL+DT_FRAC_WID+2, SIGM_TAIL+2> EXP_TYPE;
		EXP_TYPE y_0 = activation_sigm_node_nw_rom[layer_idx][uni_idx][dim_idx][base    ];
		EXP_TYPE y_1 = activation_sigm_node_nw_rom[layer_idx][uni_idx][dim_idx][base + 1];

		// Calculate (Y0 * X1) and (Y1 * X0).
		// Data types of X0 and X1 are integer.
		// Data types of Y0 and Y1 are fixed point.
		// Multiply Y0 and Y1 by each bit of X0 and X1.
		EXP_TYPE sub_0[SIGM_TAIL+1];
		EXP_TYPE sub_1[SIGM_TAIL+1];
		# pragma HLS ARRAY_PARTITION variable=sub_0 dim=1 type=complete
		# pragma HLS ARRAY_PARTITION variable=sub_1 dim=1 type=complete
		for (IDX_TYPE i=0; i<SIGM_TAIL+1; i++)
			# pragma HLS UNROLL
			if (offset_0[i] == 1) sub_0[i] = (y_0 << i);
			else sub_0[i] = 0.0;
		for (IDX_TYPE i=0; i<SIGM_TAIL+1; i++)
			# pragma HLS UNROLL
			if (offset_1[i] == 1) sub_1[i] = (y_1 << i);
			else sub_1[i] = 0.0;

		// Calculate (Y0 * X1 + Y1 * X0) / (X0 + X1).
		// Sum fixed-point numbers together.
		EXP_TYPE sum = 0.0;
		for (IDX_TYPE i=0; i<SIGM_TAIL+1; i++) {
			# pragma HLS UNROLL
			sum += sub_0[i];
			sum += sub_1[i]; }
		DT1_TYPE result = (sum >> SIGM_TAIL);

		// Output data.
		if (sign)
			return ((DT1_TYPE)1.0 - result);
		else
			return (                result);
	}
}

static DT1_TYPE act_sigm_edge_nw(
		IDX_TYPE layer_idx,
		IDX_TYPE uni_idx,
		IDX_TYPE dim_idx,
		DT0_TYPE in)
{
	# pragma HLS INLINE
	# pragma HLS ARRAY_PARTITION variable=activation_sigm_edge_nw_rom dim=1 type=complete
	# pragma HLS ARRAY_PARTITION variable=activation_sigm_edge_nw_rom dim=2 type=complete
	# pragma HLS ARRAY_PARTITION variable=activation_sigm_edge_nw_rom dim=3 type=complete

	// Whether the input data is negative.
	// If it is negative, convert it into positive.
	BOL_TYPE sign;
	DT0_TYPE pos_in;
	# pragma HLS AGGREGATE variable=pos_in
	if (in > (DT0_TYPE)0.0) {
		sign = false;
		pos_in = in; }
	else {
		sign = true;
		pos_in = -in; }

	// Divide bits of the input data into two parts.
	// The front part is the index of the look-up table.
	// The back part is used to interpolate.
	typedef ap_uint<DT0_INT_WID+DT_FRAC_WID-SIGM_TAIL-1> BASE_TYPE;
	typedef ap_uint<SIGM_TAIL+1> OFFSET_TYPE;
	BASE_TYPE base = 0;
	OFFSET_TYPE offset_0 = 0;
	OFFSET_TYPE offset_1 = 0;
	# pragma HLS AGGREGATE variable=base
	# pragma HLS AGGREGATE variable=offset_0
	# pragma HLS AGGREGATE variable=offset_1
	for (IDX_TYPE i=0; i<DT0_INT_WID+DT_FRAC_WID-SIGM_TAIL-1; i++)
		# pragma HLS UNROLL
		base[i] = pos_in[i+SIGM_TAIL];
	for (IDX_TYPE i=0; i<SIGM_TAIL; i++)
		# pragma HLS UNROLL
		offset_1[i] = pos_in[i];
	offset_0 = (1 << SIGM_TAIL) - offset_1;

	// If the index is overflowed.
	// Output the minimum or the maximum data.
	if (base >= SIGM_SIZE - 1)
	{
		if (sign)
			return (DT1_TYPE)0.0;
		else
			return (DT1_TYPE)1.0;
	}

	// If the index is not overflowed.
	// Calculate by linear interpolation algorithm.
	else
	{
		// Y(X) = (Y0 * X1 + Y1 * X0) / (X0 + X1).
		// Get the two data according to the index.
		typedef ap_fixed<SIGM_TAIL+DT_FRAC_WID+2, SIGM_TAIL+2> EXP_TYPE;
		EXP_TYPE y_0 = activation_sigm_edge_nw_rom[layer_idx][uni_idx][dim_idx][base    ];
		EXP_TYPE y_1 = activation_sigm_edge_nw_rom[layer_idx][uni_idx][dim_idx][base + 1];

		// Calculate (Y0 * X1) and (Y1 * X0).
		// Data types of X0 and X1 are integer.
		// Data types of Y0 and Y1 are fixed point.
		// Multiply Y0 and Y1 by each bit of X0 and X1.
		EXP_TYPE sub_0[SIGM_TAIL+1];
		EXP_TYPE sub_1[SIGM_TAIL+1];
		# pragma HLS ARRAY_PARTITION variable=sub_0 dim=1 type=complete
		# pragma HLS ARRAY_PARTITION variable=sub_1 dim=1 type=complete
		for (IDX_TYPE i=0; i<SIGM_TAIL+1; i++)
			# pragma HLS UNROLL
			if (offset_0[i] == 1) sub_0[i] = (y_0 << i);
			else sub_0[i] = 0.0;
		for (IDX_TYPE i=0; i<SIGM_TAIL+1; i++)
			# pragma HLS UNROLL
			if (offset_1[i] == 1) sub_1[i] = (y_1 << i);
			else sub_1[i] = 0.0;

		// Calculate (Y0 * X1 + Y1 * X0) / (X0 + X1).
		// Sum fixed-point numbers together.
		EXP_TYPE sum = 0.0;
		for (IDX_TYPE i=0; i<SIGM_TAIL+1; i++) {
			# pragma HLS UNROLL
			sum += sub_0[i];
			sum += sub_1[i]; }
		DT1_TYPE result = (sum >> SIGM_TAIL);

		// Output data.
		if (sign)
			return ((DT1_TYPE)1.0 - result);
		else
			return (                result);
	}
}

static DT1_TYPE act_tanh_node_nw(
		IDX_TYPE layer_idx,
		IDX_TYPE uni_idx,
		IDX_TYPE dim_idx,
		DT0_TYPE in)
{
	# pragma HLS INLINE
	# pragma HLS ARRAY_PARTITION variable=activation_tanh_node_nw_rom dim=1 type=complete
	# pragma HLS ARRAY_PARTITION variable=activation_tanh_node_nw_rom dim=2 type=complete
	# pragma HLS ARRAY_PARTITION variable=activation_tanh_node_nw_rom dim=3 type=complete

	// Whether the input data is negative.
	// If it is negative, convert it into positive.
	BOL_TYPE sign;
	DT0_TYPE pos_in;
	# pragma HLS AGGREGATE variable=pos_in
	if (in > (DT0_TYPE)0.0) {
		sign = false;
		pos_in = in; }
	else {
		sign = true;
		pos_in = -in; }

	// Divide bits of the input data into two parts.
	// The front part is the index of the look-up table.
	// The back part is used to interpolate.
	typedef ap_uint<DT0_INT_WID+DT_FRAC_WID-TANH_TAIL-1> BASE_TYPE;
	typedef ap_uint<TANH_TAIL+1> OFFSET_TYPE;
	BASE_TYPE base = 0;
	OFFSET_TYPE offset_0 = 0;
	OFFSET_TYPE offset_1 = 0;
	# pragma HLS AGGREGATE variable=base
	# pragma HLS AGGREGATE variable=offset_0
	# pragma HLS AGGREGATE variable=offset_1
	for (IDX_TYPE i=0; i<DT0_INT_WID+DT_FRAC_WID-TANH_TAIL-1; i++)
		# pragma HLS UNROLL
		base[i] = pos_in[i+TANH_TAIL];
	for (IDX_TYPE i=0; i<TANH_TAIL; i++)
		# pragma HLS UNROLL
		offset_1[i] = pos_in[i];
	offset_0 = (1 << TANH_TAIL) - offset_1;

	// If the index is overflowed.
	// Output the minimum or the maximum data.
	if (base >= TANH_SIZE - 1)
	{
		if (sign)
			return (DT1_TYPE)(-1.0);
		else
			return (DT1_TYPE)( 1.0);
	}

	// If the index is not overflowed.
	// Calculate by linear interpolation algorithm.
	else
	{
		// Y(X) = (Y0 * X1 + Y1 * X0) / (X0 + X1).
		// Get the two data according to the index.
		typedef ap_fixed<TANH_TAIL+DT_FRAC_WID+2, TANH_TAIL+2> EXP_TYPE;
		EXP_TYPE y_0 = activation_tanh_node_nw_rom[layer_idx][uni_idx][dim_idx][base    ];
		EXP_TYPE y_1 = activation_tanh_node_nw_rom[layer_idx][uni_idx][dim_idx][base + 1];

		// Calculate (Y0 * X1) and (Y1 * X0).
		// Data types of X0 and X1 are integer.
		// Data types of Y0 and Y1 are fixed point.
		// Multiply Y0 and Y1 by each bit of X0 and X1.
		EXP_TYPE sub_0[TANH_TAIL+1];
		EXP_TYPE sub_1[TANH_TAIL+1];
		# pragma HLS ARRAY_PARTITION variable=sub_0 dim=1 type=complete
		# pragma HLS ARRAY_PARTITION variable=sub_1 dim=1 type=complete
		for (IDX_TYPE i=0; i<TANH_TAIL+1; i++)
			# pragma HLS UNROLL
			if (offset_0[i] == 1) sub_0[i] = (y_0 << i);
			else sub_0[i] = 0.0;
		for (IDX_TYPE i=0; i<TANH_TAIL+1; i++)
			# pragma HLS UNROLL
			if (offset_1[i] == 1) sub_1[i] = (y_1 << i);
			else sub_1[i] = 0.0;

		// Calculate (Y0 * X1 + Y1 * X0) / (X0 + X1).
		// Sum fixed-point numbers together.
		EXP_TYPE sum = 0.0;
		for (IDX_TYPE i=0; i<TANH_TAIL+1; i++) {
			# pragma HLS UNROLL
			sum += sub_0[i];
			sum += sub_1[i]; }
		DT1_TYPE result = (sum >> TANH_TAIL);

		// Output data.
		if (sign)
			return ( - result);
		else
			return (   result);
	}
}

static DT1_TYPE act_tanh_edge_nw(
		IDX_TYPE layer_idx,
		IDX_TYPE uni_idx,
		IDX_TYPE dim_idx,
		DT0_TYPE in)
{
	# pragma HLS INLINE
	# pragma HLS ARRAY_PARTITION variable=activation_tanh_edge_nw_rom dim=1 type=complete
	# pragma HLS ARRAY_PARTITION variable=activation_tanh_edge_nw_rom dim=2 type=complete
	# pragma HLS ARRAY_PARTITION variable=activation_tanh_edge_nw_rom dim=3 type=complete

	// Whether the input data is negative.
	// If it is negative, convert it into positive.
	BOL_TYPE sign;
	DT0_TYPE pos_in;
	# pragma HLS AGGREGATE variable=pos_in
	if (in > (DT0_TYPE)0.0) {
		sign = false;
		pos_in = in; }
	else {
		sign = true;
		pos_in = -in; }

	// Divide bits of the input data into two parts.
	// The front part is the index of the look-up table.
	// The back part is used to interpolate.
	typedef ap_uint<DT0_INT_WID+DT_FRAC_WID-TANH_TAIL-1> BASE_TYPE;
	typedef ap_uint<TANH_TAIL+1> OFFSET_TYPE;
	BASE_TYPE base = 0;
	OFFSET_TYPE offset_0 = 0;
	OFFSET_TYPE offset_1 = 0;
	# pragma HLS AGGREGATE variable=base
	# pragma HLS AGGREGATE variable=offset_0
	# pragma HLS AGGREGATE variable=offset_1
	for (IDX_TYPE i=0; i<DT0_INT_WID+DT_FRAC_WID-TANH_TAIL-1; i++)
		# pragma HLS UNROLL
		base[i] = pos_in[i+TANH_TAIL];
	for (IDX_TYPE i=0; i<TANH_TAIL; i++)
		# pragma HLS UNROLL
		offset_1[i] = pos_in[i];
	offset_0 = (1 << TANH_TAIL) - offset_1;

	// If the index is overflowed.
	// Output the minimum or the maximum data.
	if (base >= TANH_SIZE - 1)
	{
		if (sign)
			return (DT1_TYPE)(-1.0);
		else
			return (DT1_TYPE)( 1.0);
	}

	// If the index is not overflowed.
	// Calculate by linear interpolation algorithm.
	else
	{
		// Y(X) = (Y0 * X1 + Y1 * X0) / (X0 + X1).
		// Get the two data according to the index.
		typedef ap_fixed<TANH_TAIL+DT_FRAC_WID+2, TANH_TAIL+2> EXP_TYPE;
		EXP_TYPE y_0 = activation_tanh_edge_nw_rom[layer_idx][uni_idx][dim_idx][base    ];
		EXP_TYPE y_1 = activation_tanh_edge_nw_rom[layer_idx][uni_idx][dim_idx][base + 1];

		// Calculate (Y0 * X1) and (Y1 * X0).
		// Data types of X0 and X1 are integer.
		// Data types of Y0 and Y1 are fixed point.
		// Multiply Y0 and Y1 by each bit of X0 and X1.
		EXP_TYPE sub_0[TANH_TAIL+1];
		EXP_TYPE sub_1[TANH_TAIL+1];
		# pragma HLS ARRAY_PARTITION variable=sub_0 dim=1 type=complete
		# pragma HLS ARRAY_PARTITION variable=sub_1 dim=1 type=complete
		for (IDX_TYPE i=0; i<TANH_TAIL+1; i++)
			# pragma HLS UNROLL
			if (offset_0[i] == 1) sub_0[i] = (y_0 << i);
			else sub_0[i] = 0.0;
		for (IDX_TYPE i=0; i<TANH_TAIL+1; i++)
			# pragma HLS UNROLL
			if (offset_1[i] == 1) sub_1[i] = (y_1 << i);
			else sub_1[i] = 0.0;

		// Calculate (Y0 * X1 + Y1 * X0) / (X0 + X1).
		// Sum fixed-point numbers together.
		EXP_TYPE sum = 0.0;
		for (IDX_TYPE i=0; i<TANH_TAIL+1; i++) {
			# pragma HLS UNROLL
			sum += sub_0[i];
			sum += sub_1[i]; }
		DT1_TYPE result = (sum >> TANH_TAIL);

		// Output data.
		if (sign)
			return ( - result);
		else
			return (   result);
	}
}



// **************************** Input Network **************************** //

/*

	This part is to implement the input network. For a certain node,
	it gets all the features, multiply them with weights, and sum
	them together to form one dimension of the node embeddings. It
	then sends the dimension to the adapter, and calculate the next
	dimension.

	Version: 0.1
		Author: Hanqing Zhang
		Email : hanqing.zhang@zju.edu.cn
		Date  : June-29-2023
		Commit: Add function "input_nw".

	Version: 0.2
		Author: Hanqing Zhang
		Email : hanqing.zhang@zju.edu.cn
		Date  : July-2-2023
		Commit: Modify the output-data type to hls::stream.

	Version: 0.3
		Author: Hanqing Zhang
		Email : hanqing.zhang@zju.edu.cn
		Date  : July-13-2023
		Commit: Modify the unit because the original one cannot be
				synthesized.

	Version: 0.4
		Author: Hanqing Zhang
		Email : hanqing.zhang@zju.edu.cn
		Date  : July-23-2023
		Commit: Add an extra output FIFO for the bypass channel.

	Version: 0.5
		Author: Hanqing Zhang
		Email : hanqing.zhang@zju.edu.cn
		Date  : October-16-2023
		Commit: Fixed the bug that the last layer cannot be 
				activated.

	Function:
		void input_nw(
			BOL_TYPE last_layer,
			IDX_TYPE num_of_nodes,
			hls::stream<DT1_PARA> mlp_out_0[NODE_UNI_PARA],
			hls::stream<DT1_PARA> mlp_out_1[NODE_UNI_PARA]);
		Input arguments:
			last_layer: flag of the last layer.
			num_of_nodes: the number of nodes.
		Output arguments:
			mlp_out_0: node embeddings sent into the adapter.
			mlp_out_1: node embeddings sent into the bypass channel.

*/

static void input_nw(
		BOL_TYPE last_layer,
		IDX_TYPE num_of_nodes,
		hls::stream<DT1_PARA> mlp_out_0[NODE_UNI_PARA],
		hls::stream<DT1_PARA> mlp_out_1[NODE_UNI_PARA])
{
	# pragma HLS INLINE off
	# pragma HLS ARRAY_PARTITION variable=g_input_nw_weight dim=1 type=cyclic factor=DIM_PARA
	# pragma HLS ARRAY_PARTITION variable=g_input_nw_weight dim=2 type=complete
	# pragma HLS ARRAY_PARTITION variable=g_input_nw_bias dim=1 type=cyclic factor=DIM_PARA

	// Process NODE_UNI_PARA nodes each time.
	// Process DIM_PARA dimensions each time.
	IDX_TYPE iter_uni_num = ceildiv(num_of_nodes, (IDX_TYPE)NODE_UNI_PARA);
	IDX_TYPE iter_dim_num = ceildiv((IDX_TYPE)NODE_EMB_DIM, (IDX_TYPE)DIM_PARA);
	for (	IDX_TYPE iter_uni_cnt=(IDX_TYPE)0,
			node_idx_base=(IDX_TYPE)0;
			iter_uni_cnt<iter_uni_num;
			iter_uni_cnt+=(IDX_TYPE)1,
			node_idx_base+=(IDX_TYPE)NODE_UNI_PARA)
	# pragma HLS LOOP_TRIPCOUNT min=NODE_UNI_COUNT max=NODE_UNI_COUNT avg=NODE_UNI_COUNT
	for (	IDX_TYPE iter_dim_cnt=(IDX_TYPE)0,
			dim_out_base=(IDX_TYPE)0;
			iter_dim_cnt<iter_dim_num;
			iter_dim_cnt+=(IDX_TYPE)1,
			dim_out_base+=(IDX_TYPE)DIM_PARA)
	{
		# pragma HLS LOOP_TRIPCOUNT min=DIM_COUNT max=DIM_COUNT avg=DIM_COUNT
		# pragma HLS PIPELINE II=1 style=frp

		// Registers that stores products and outputs of the linear layer.
		DT0_TYPE linear_out[NODE_UNI_PARA][DIM_PARA];
		DT0_TYPE product[NODE_UNI_PARA][DIM_PARA][NODE_FEA_DIM];
		# pragma HLS ARRAY_PARTITION variable=linear_out dim=1 type=complete
		# pragma HLS ARRAY_PARTITION variable=linear_out dim=2 type=complete
		# pragma HLS ARRAY_PARTITION variable=product dim=1 type=complete
		# pragma HLS ARRAY_PARTITION variable=product dim=2 type=complete
		# pragma HLS ARRAY_PARTITION variable=product dim=3 type=complete

		// Get the bias of a certain output dimension.
		for (IDX_TYPE dim_out_offset=0; dim_out_offset<DIM_PARA; dim_out_offset++) {
			# pragma HLS UNROLL
			// Get the index of the output dimension.
			IDX_TYPE dim_out = dim_out_base + dim_out_offset;
			BOL_TYPE dim_out_overflow = (dim_out >= (IDX_TYPE)NODE_EMB_DIM);
			// Get the bias.
			DT1_TYPE bias = dim_out_overflow? (DT1_TYPE)0.0: g_input_nw_bias[dim_out];
			for (IDX_TYPE node_idx_offset=0; node_idx_offset<NODE_UNI_PARA; node_idx_offset++)
			# pragma HLS UNROLL
			linear_out[node_idx_offset][dim_out_offset] = bias; }

		// Parallel calculate the products of features and weights.
		for (IDX_TYPE dim_in=0; dim_in<NODE_FEA_DIM; dim_in++)
		# pragma HLS UNROLL
		for (IDX_TYPE dim_out_offset=0; dim_out_offset<DIM_PARA; dim_out_offset++) {
			# pragma HLS UNROLL
			// Get the index of the output dimension.
			IDX_TYPE dim_out = dim_out_base + dim_out_offset;
			BOL_TYPE dim_out_overflow = (dim_out >= (IDX_TYPE)NODE_EMB_DIM);
			// Get the weight.
			DT1_TYPE weight = dim_out_overflow? (DT1_TYPE)0.0: g_input_nw_weight[dim_out][dim_in];
			// Calculate products.
			for (IDX_TYPE node_idx_offset=0; node_idx_offset<NODE_UNI_PARA; node_idx_offset++) {
				# pragma HLS UNROLL
				// Get the index of the node.
				IDX_TYPE node_idx = node_idx_base + node_idx_offset;
				BOL_TYPE node_idx_overflow = (node_idx >= num_of_nodes);
				// Get the feature of the input node.
				DT0_TYPE input = node_idx_overflow? (DT0_TYPE)0.0: g_node_feature[node_idx][dim_in];
				// Calculate and store the product.
				product[node_idx_offset][dim_out_offset][dim_in] = weight * input; }}

		// Sum NODE_FEA_DIM products together.
		for (IDX_TYPE dim_in=0; dim_in<NODE_FEA_DIM; dim_in++)
		# pragma HLS UNROLL
		for (IDX_TYPE dim_out_offset=0; dim_out_offset<DIM_PARA; dim_out_offset++)
		# pragma HLS UNROLL
		for (IDX_TYPE node_idx_offset=0; node_idx_offset<NODE_UNI_PARA; node_idx_offset++)
			# pragma HLS UNROLL
			linear_out[node_idx_offset][dim_out_offset] += product[node_idx_offset][dim_out_offset][dim_in];

		// Go through the activation function and write it into FIFOs.
		for (IDX_TYPE node_idx_offset=0; node_idx_offset<NODE_UNI_PARA; node_idx_offset++) {
			# pragma HLS UNROLL
			// Data output by the activation function.
			DT1_PARA act_out;
			# pragma HLS AGGREGATE variable=act_out
			// Pass the data through the activation function.
			for (IDX_TYPE dim_out_offset=0; dim_out_offset<DIM_PARA; dim_out_offset++) {
				# pragma HLS UNROLL
				// Activation function.
				DT0_TYPE act_in = linear_out[node_idx_offset][dim_out_offset];
				# ifdef ACT_TANH
				act_out[dim_out_offset] = act_tanh_node_nw(0, node_idx_offset, dim_out_offset, act_in); }
				# endif
				# ifdef ACT_SIGM
				act_out[dim_out_offset] = act_sigm_node_nw(0, node_idx_offset, dim_out_offset, act_in); }
				# endif
			// Write data into FIFOs.
			if (!last_layer)
			mlp_out_1[node_idx_offset].write(act_out);
			mlp_out_0[node_idx_offset].write(act_out); }
	}
}



// **************************** Node-Processing Network **************************** //

/*

	This part is to implement node-processing network. It contains
	4 linear layers and 4 activation layers. For processing a node,
	the linear layer first fetches one input dimension, products it
	with corresponding weights, and sends them into all the output
	dimensions parallel. After iterating all the input dimension,
	the linear layer send all the output dimensions to the activation
	layer parallel. Once the activation layer receive the data. It
	will divide all the dimensions separately. Then, it'll make each
	data pass through the activation function, and send it into the
	next layer serially.

	Version: 0.1
		Author: Hanqing Zhang
		Email : hanqing.zhang@zju.edu.cn
		Date  : June-30-2023
		Commit: Add function "node_nw".

	Version: 0.2
		Author: Hanqing Zhang
		Email : hanqing.zhang@zju.edu.cn
		Date  : July-2-2023
		Commit: Modify the output-data type.

	Version: 0.3
		Author: Hanqing Zhang
		Email : hanqing.zhang@zju.edu.cn
		Data  : July-17-2023
		Commit: Modify for the DATAFLOW pragma.

	Version: 0.4
		Author: Hanqing Zhang
		Email : hanqing.zhang@zju.edu.cn
		Data  : September-20-2023
		Commit: Modify for the dimension parallel.

	Version: 0.5
		Author: Hanqing Zhang
		Email : hanqing.zhang@zju.edu.cn
		Data  : October-10-2023
		Commit: Fix the bug that the last layer cannot be activated.

	Version: 0.6
		Author: Hanqing Zhang
		Email : hanqing.zhang@zju.edu.cn
		Data  : October-19-2023
		Commit: Modify the depth of the FIFOs.

	Function:
		void node_nw(
			BOL_TYPE last_layer,
			IDX_TYPE num_of_nodes,
			DT1_TYPE node_emb[NODE_NUM_MAX][NODE_EMB_DIM],
			DT0_TYPE msg_dst[2][EDGE_UNI_PARA][NODE_NUM_MAX][NODE_EMB_DIM],
			DT0_TYPE msg_src[2][EDGE_UNI_PARA][NODE_NUM_MAX][NODE_EMB_DIM],
			IDX_TYPE msg_dst_flag[EDGE_UNI_PARA][NODE_NUM_MAX],
			IDX_TYPE msg_src_flag[EDGE_UNI_PARA][NODE_NUM_MAX],
			hls::stream<DT1_PARA> mlp_out_0[NODE_UNI_PARA],
			hls::stream<DT1_PARA> mlp_out_1[NODE_UNI_PARA]);
		Input arguments:
			last_layer: flag of the last layer.
			num_of_nodes: the number of nodes.
			node_emb: node embedding array.
			msg_dst: destination message array.
			msg_src: source message array.
			msg_dst_flag: flags of where the destination messages stored.
			msg_src_flag: flags of where the source messages stored.
		Output arguments:
			mlp_out_0: New node embeddings to the adapter.
			mlp_out_1: Bypass data to the node scattering module.

*/

static void node_linear_0(
		IDX_TYPE num_of_nodes,
		DT0_TYPE msg_dst[2][EDGE_UNI_PARA][NODE_NUM_MAX][NODE_EMB_DIM],
		DT0_TYPE msg_src[2][EDGE_UNI_PARA][NODE_NUM_MAX][NODE_EMB_DIM],
		IDX_TYPE msg_dst_flag[EDGE_UNI_PARA][NODE_NUM_MAX],
		IDX_TYPE msg_src_flag[EDGE_UNI_PARA][NODE_NUM_MAX],
		DT1_TYPE node_emb[NODE_NUM_MAX][NODE_EMB_DIM],
		hls::stream<DT0_VECTOR> linear_out[NODE_UNI_PARA],
		hls::stream<DT1_PARA> node_emb_old_out[NODE_UNI_PARA])
{
	# pragma HLS INLINE off
	# pragma HLS ARRAY_PARTITION variable=g_node_nw_weight_0_0 dim=1 type=complete
	# pragma HLS ARRAY_PARTITION variable=g_node_nw_weight_0_1 dim=1 type=complete
	# pragma HLS ARRAY_PARTITION variable=g_node_nw_weight_0_2 dim=1 type=complete
	# pragma HLS ARRAY_PARTITION variable=g_node_nw_weight_0_0 dim=2 type=cyclic factor=DIM_PARA
	# pragma HLS ARRAY_PARTITION variable=g_node_nw_weight_0_1 dim=2 type=cyclic factor=DIM_PARA
	# pragma HLS ARRAY_PARTITION variable=g_node_nw_weight_0_2 dim=2 type=cyclic factor=DIM_PARA
	# pragma HLS ARRAY_PARTITION variable=g_node_nw_bias_0 dim=1 type=complete

	// Temporary array of linear outputs.
	DT0_VECTOR linear_out_cache[NODE_UNI_PARA];
	# pragma HLS AGGREGATE variable=linear_out_cache
	# pragma HLS ARRAY_PARTITION dim=1 type=complete variable=linear_out_cache

	// Temporary array of flags.
	IDX_TYPE msg_dst_flag_cache[EDGE_UNI_PARA][NODE_UNI_PARA];
	IDX_TYPE msg_src_flag_cache[EDGE_UNI_PARA][NODE_UNI_PARA];
	# pragma HLS ARRAY_PARTITION dim=1 type=complete variable=msg_dst_flag_cache
	# pragma HLS ARRAY_PARTITION dim=2 type=complete variable=msg_dst_flag_cache
	# pragma HLS ARRAY_PARTITION dim=1 type=complete variable=msg_src_flag_cache
	# pragma HLS ARRAY_PARTITION dim=2 type=complete variable=msg_src_flag_cache

	// Process NODE_UNI_PARA nodes each time.
	// Process DIM_PARA dimensions each time.
	IDX_TYPE iter_uni_num = ceildiv(num_of_nodes, (IDX_TYPE)NODE_UNI_PARA);
	IDX_TYPE iter_dim_num = ceildiv((IDX_TYPE)NODE_EMB_DIM, (IDX_TYPE)DIM_PARA);
	for (	IDX_TYPE iter_uni_cnt=(IDX_TYPE)0,
			node_idx_base=(IDX_TYPE)0;
			iter_uni_cnt<iter_uni_num;
			iter_uni_cnt+=(IDX_TYPE)1,
			node_idx_base+=(IDX_TYPE)NODE_UNI_PARA)
	# pragma HLS LOOP_TRIPCOUNT min=NODE_UNI_COUNT max=NODE_UNI_COUNT avg=NODE_UNI_COUNT
	for (	IDX_TYPE iter_dim_cnt=(IDX_TYPE)0,
			dim_in_base=(IDX_TYPE)0;
			iter_dim_cnt<iter_dim_num;
			iter_dim_cnt+=(IDX_TYPE)1,
			dim_in_base+=(IDX_TYPE)DIM_PARA)
	{
		# pragma HLS LOOP_TRIPCOUNT min=DIM_COUNT max=DIM_COUNT avg=DIM_COUNT
		# pragma HLS PIPELINE II=1 style=frp
		# pragma HLS DEPENDENCE variable=msg_dst_flag inter true distance=DIM_COUNT
		# pragma HLS DEPENDENCE variable=msg_src_flag inter true distance=DIM_COUNT

		// Whether it is the first or the last iteration.
		BOL_TYPE first = (dim_in_base == (IDX_TYPE)0);
		BOL_TYPE last = (dim_in_base + (IDX_TYPE)DIM_PARA >= (IDX_TYPE)NODE_EMB_DIM);

		// If it is the first iteration.
		if (first) {
			// Get biases of all the output dimensions.
			for (IDX_TYPE dim_out=0; dim_out<NODE_EMB_DIM; dim_out++) {
				# pragma HLS UNROLL
				DT1_TYPE bias = g_node_nw_bias_0[dim_out];
				for (IDX_TYPE node_idx_offset=0; node_idx_offset<NODE_UNI_PARA; node_idx_offset++)
					# pragma HLS UNROLL
					linear_out_cache[node_idx_offset][dim_out] = bias; }

			// Get and update flags.
			for (IDX_TYPE node_idx_offset=0; node_idx_offset<NODE_UNI_PARA; node_idx_offset++) {
				# pragma HLS UNROLL
				// Get the certain node index.
				IDX_TYPE node_idx = node_idx_base + node_idx_offset;
				BOL_TYPE node_idx_overflow = (node_idx >= num_of_nodes);
				// Save flags to temporary arrays.
				if (!node_idx_overflow)
				for (IDX_TYPE edge_idx_offset=0; edge_idx_offset<EDGE_UNI_PARA; edge_idx_offset++) {
					# pragma HLS UNROLL
					msg_dst_flag_cache[edge_idx_offset][node_idx_offset] = msg_dst_flag[edge_idx_offset][node_idx];
					msg_src_flag_cache[edge_idx_offset][node_idx_offset] = msg_src_flag[edge_idx_offset][node_idx];
					msg_dst_flag[edge_idx_offset][node_idx] = (IDX_TYPE)-1;
					msg_src_flag[edge_idx_offset][node_idx] = (IDX_TYPE)-1; }
				else
				for (IDX_TYPE edge_idx_offset=0; edge_idx_offset<EDGE_UNI_PARA; edge_idx_offset++) {
					# pragma HLS UNROLL
					msg_dst_flag_cache[edge_idx_offset][node_idx_offset] = (IDX_TYPE)-1;
					msg_src_flag_cache[edge_idx_offset][node_idx_offset] = (IDX_TYPE)-1; }}}

		// Process NODE_UNI_PARA nodes each time.
		for (IDX_TYPE node_idx_offset=0; node_idx_offset<NODE_UNI_PARA; node_idx_offset++) {
			# pragma HLS UNROLL
			// Get the index of the certain node.
			IDX_TYPE node_idx = node_idx_base + node_idx_offset;
			BOL_TYPE node_idx_overflow = (node_idx >= num_of_nodes);
			// Temporary array of input messages.
			DT0_PARA input_msg_dst;
			DT0_PARA input_msg_src;
			# pragma HLS AGGREGATE variable=input_msg_dst
			# pragma HLS AGGREGATE variable=input_msg_src
			for (IDX_TYPE dim_in_offset=0; dim_in_offset<DIM_PARA; dim_in_offset++) {
				# pragma HLS UNROLL
				input_msg_dst[dim_in_offset] = (DT0_TYPE)0.0;
				input_msg_src[dim_in_offset] = (DT0_TYPE)0.0; }

			// Fetch messages.
			if (!node_idx_overflow)
			for (IDX_TYPE dim_in_offset=0; dim_in_offset<DIM_PARA; dim_in_offset++) {
				# pragma HLS UNROLL
				// Get the index of the input dimension.
				IDX_TYPE dim_in = dim_in_base + dim_in_offset;
				BOL_TYPE dim_in_overflow = (dim_in >= (IDX_TYPE)NODE_EMB_DIM);
				// Fetch messages according to flags.
				for (IDX_TYPE edge_idx_offset=0; edge_idx_offset<EDGE_UNI_PARA; edge_idx_offset++) {
					# pragma HLS UNROLL
					IDX_TYPE dst_flag = msg_dst_flag_cache[edge_idx_offset][node_idx_offset];
					IDX_TYPE src_flag = msg_src_flag_cache[edge_idx_offset][node_idx_offset];
					input_msg_dst[dim_in_offset] +=
							((dst_flag == (IDX_TYPE)-1) || dim_in_overflow)? (DT0_TYPE)0.0:
							((dst_flag == (IDX_TYPE)0)?
							msg_dst[0][edge_idx_offset][node_idx][dim_in]:
							msg_dst[1][edge_idx_offset][node_idx][dim_in]);
					input_msg_src[dim_in_offset] +=
							((src_flag == (IDX_TYPE)-1) || dim_in_overflow)? (DT0_TYPE)0.0:
							((src_flag == (IDX_TYPE)0)?
							msg_src[0][edge_idx_offset][node_idx][dim_in]:
							msg_src[1][edge_idx_offset][node_idx][dim_in]); }}

			// Fetch old embedding.
			DT1_PARA input_node_emb;
			for (IDX_TYPE dim_in_offset=0; dim_in_offset<DIM_PARA; dim_in_offset++) {
				# pragma HLS UNROLL
				// Get the index of the input dimension.
				IDX_TYPE dim_in = dim_in_base + dim_in_offset;
				BOL_TYPE dim_in_overflow = (dim_in >= (IDX_TYPE)NODE_EMB_DIM);
				// Fetch old embedding.
				input_node_emb[dim_in_offset] = node_idx_overflow?
						(DT1_TYPE)0.0: node_emb[node_idx][dim_in]; }
			// Write old embedding into bypass channel.
			node_emb_old_out[node_idx_offset].write(input_node_emb);

			// Parallel process all the output dimension.
			for (IDX_TYPE dim_out=0; dim_out<NODE_EMB_DIM; dim_out++) {
				# pragma HLS UNROLL
				// Temporary variables of products.
				DT0_TYPE product_msg_dst = (DT0_TYPE)0.0;
				DT0_TYPE product_msg_src = (DT0_TYPE)0.0;
				DT0_TYPE product_node_emb = (DT0_TYPE)0.0;
				// Parallel process all the input dimension.
				for (IDX_TYPE dim_in_offset=0; dim_in_offset<DIM_PARA; dim_in_offset++) {
					# pragma HLS UNROLL
					// Get the index of the input dimension.
					IDX_TYPE dim_in = dim_in_base + dim_in_offset;
					BOL_TYPE dim_in_overflow = (dim_in >= (IDX_TYPE)NODE_EMB_DIM);
					// Calculate products.
					DT0_TYPE product_msg_dst_cache = dim_in_overflow? (DT0_TYPE)0.0:
							(DT0_TYPE)(input_msg_dst[dim_in_offset] * g_node_nw_weight_0_0[dim_out][dim_in]);
					DT0_TYPE product_msg_src_cache = dim_in_overflow? (DT0_TYPE)0.0:
							(DT0_TYPE)(input_msg_src[dim_in_offset] * g_node_nw_weight_0_1[dim_out][dim_in]);
					DT0_TYPE product_node_emb_cache = dim_in_overflow? (DT0_TYPE)0.0:
							(DT0_TYPE)(input_node_emb[dim_in_offset] * g_node_nw_weight_0_2[dim_out][dim_in]);
					// Sum them together.
					product_msg_dst += product_msg_dst_cache;
					product_msg_src += product_msg_src_cache;
					product_node_emb += product_node_emb_cache; }
				// Get the old data and sum them together.
				DT0_TYPE data_old = linear_out_cache[node_idx_offset][dim_out];
				DT0_TYPE data_new = product_msg_dst + product_msg_src + product_node_emb + data_old;
				linear_out_cache[node_idx_offset][dim_out] = data_new; }}

		// If it is the last iteration.
		if (last)
		for (IDX_TYPE node_idx_offset=0; node_idx_offset<NODE_UNI_PARA; node_idx_offset++)
			# pragma HLS UNROLL
			linear_out[node_idx_offset].write(linear_out_cache[node_idx_offset]);
	}
}

static void node_linear_1(
		IDX_TYPE num_of_nodes,
		hls::stream<DT1_PARA> linear_in[NODE_UNI_PARA],
		hls::stream<DT0_VECTOR> linear_out[NODE_UNI_PARA],
		hls::stream<DT1_PARA> node_emb_old_in[NODE_UNI_PARA],
		hls::stream<DT1_PARA> node_emb_old_out[NODE_UNI_PARA])
{
	# pragma HLS INLINE off
	# pragma HLS ARRAY_PARTITION variable=g_node_nw_weight_1 dim=1 type=complete
	# pragma HLS ARRAY_PARTITION variable=g_node_nw_weight_1 dim=2 type=cyclic factor=DIM_PARA
	# pragma HLS ARRAY_PARTITION variable=g_node_nw_bias_1 dim=1 type=complete

	// Temporary array of linear outputs.
	DT0_VECTOR linear_out_cache[NODE_UNI_PARA];
	# pragma HLS AGGREGATE variable=linear_out_cache
	# pragma HLS ARRAY_PARTITION dim=1 type=complete variable=linear_out_cache

	// Process NODE_UNI_PARA nodes each time.
	// Process DIM_PARA dimensions each time.
	IDX_TYPE iter_uni_num = ceildiv(num_of_nodes, (IDX_TYPE)NODE_UNI_PARA);
	IDX_TYPE iter_dim_num = ceildiv((IDX_TYPE)NODE_EMB_DIM, (IDX_TYPE)DIM_PARA);
	for (	IDX_TYPE iter_uni_cnt=(IDX_TYPE)0,
			node_idx_base=(IDX_TYPE)0;
			iter_uni_cnt<iter_uni_num;
			iter_uni_cnt+=(IDX_TYPE)1,
			node_idx_base+=(IDX_TYPE)NODE_UNI_PARA)
	# pragma HLS LOOP_TRIPCOUNT min=NODE_UNI_COUNT max=NODE_UNI_COUNT avg=NODE_UNI_COUNT
	for (	IDX_TYPE iter_dim_cnt=(IDX_TYPE)0,
			dim_in_base=(IDX_TYPE)0;
			iter_dim_cnt<iter_dim_num;
			iter_dim_cnt+=(IDX_TYPE)1,
			dim_in_base+=(IDX_TYPE)DIM_PARA)
	{
		# pragma HLS LOOP_TRIPCOUNT min=DIM_COUNT max=DIM_COUNT avg=DIM_COUNT
		# pragma HLS PIPELINE II=1 style=frp

		// Whether it is the first or the last iteration.
		BOL_TYPE first = (dim_in_base == (IDX_TYPE)0);
		BOL_TYPE last = (dim_in_base + (IDX_TYPE)DIM_PARA >= (IDX_TYPE)NODE_EMB_DIM);

		// Get biases of all the output dimensions.
		if (first)
		for (IDX_TYPE dim_out=0; dim_out<NODE_EMB_DIM; dim_out++) {
			# pragma HLS UNROLL
			DT1_TYPE bias = g_node_nw_bias_1[dim_out];
			for (IDX_TYPE node_idx_offset=0; node_idx_offset<NODE_UNI_PARA; node_idx_offset++)
				# pragma HLS UNROLL
				linear_out_cache[node_idx_offset][dim_out] = bias; }

		// Process NODE_UNI_PARA nodes each time.
		for (IDX_TYPE node_idx_offset=0; node_idx_offset<NODE_UNI_PARA; node_idx_offset++) {
			# pragma HLS UNROLL

			// Fetch input data.
			DT1_PARA input;
			# pragma HLS AGGREGATE variable=input
			linear_in[node_idx_offset].read(input);

			// Bypass channel.
			DT1_PARA emb;
			# pragma HLS AGGREGATE variable=emb
			node_emb_old_in[node_idx_offset].read(emb);
			node_emb_old_out[node_idx_offset].write(emb);

			// Get the weight and calculate product.
			for (IDX_TYPE dim_out=0; dim_out<NODE_EMB_DIM; dim_out++) {
				# pragma HLS UNROLL
				// Calculate the sum of products.
				DT0_TYPE sum_of_product = (DT0_TYPE)0.0;
				for (IDX_TYPE dim_in_offset=0; dim_in_offset<DIM_PARA; dim_in_offset++) {
					# pragma HLS UNROLL
					IDX_TYPE dim_in = dim_in_base + dim_in_offset;
					IDX_TYPE dim_in_overflow = (dim_in >= (IDX_TYPE)NODE_EMB_DIM);
					DT1_TYPE product_cache = dim_in_overflow? (DT1_TYPE)0.0:
							(DT1_TYPE)(input[dim_in_offset] * g_node_nw_weight_1[dim_out][dim_in]);
					sum_of_product += product_cache; }
				// Get the old data and sum them together.
				DT0_TYPE data_old = linear_out_cache[node_idx_offset][dim_out];
				DT0_TYPE data_new = sum_of_product + data_old;
				linear_out_cache[node_idx_offset][dim_out] = data_new; }}

		// Send the data into the FIFO.
		if (last)
		for (IDX_TYPE node_idx_offset=0; node_idx_offset<NODE_UNI_PARA; node_idx_offset++)
			# pragma HLS UNROLL
			linear_out[node_idx_offset].write(linear_out_cache[node_idx_offset]);
	}
}

static void node_linear_2(
		IDX_TYPE num_of_nodes,
		hls::stream<DT1_PARA> linear_in[NODE_UNI_PARA],
		hls::stream<DT0_VECTOR> linear_out[NODE_UNI_PARA],
		hls::stream<DT1_PARA> node_emb_old_in[NODE_UNI_PARA],
		hls::stream<DT1_PARA> node_emb_old_out[NODE_UNI_PARA])
{
	# pragma HLS INLINE off
	# pragma HLS ARRAY_PARTITION variable=g_node_nw_weight_2 dim=1 type=complete
	# pragma HLS ARRAY_PARTITION variable=g_node_nw_weight_2 dim=2 type=cyclic factor=DIM_PARA
	# pragma HLS ARRAY_PARTITION variable=g_node_nw_bias_2 dim=1 type=complete

	// Temporary array of linear outputs.
	DT0_VECTOR linear_out_cache[NODE_UNI_PARA];
	# pragma HLS AGGREGATE variable=linear_out_cache
	# pragma HLS ARRAY_PARTITION dim=1 type=complete variable=linear_out_cache

	// Process NODE_UNI_PARA nodes each time.
	// Process DIM_PARA dimensions each time.
	IDX_TYPE iter_uni_num = ceildiv(num_of_nodes, (IDX_TYPE)NODE_UNI_PARA);
	IDX_TYPE iter_dim_num = ceildiv((IDX_TYPE)NODE_EMB_DIM, (IDX_TYPE)DIM_PARA);
	for (	IDX_TYPE iter_uni_cnt=(IDX_TYPE)0,
			node_idx_base=(IDX_TYPE)0;
			iter_uni_cnt<iter_uni_num;
			iter_uni_cnt+=(IDX_TYPE)1,
			node_idx_base+=(IDX_TYPE)NODE_UNI_PARA)
	# pragma HLS LOOP_TRIPCOUNT min=NODE_UNI_COUNT max=NODE_UNI_COUNT avg=NODE_UNI_COUNT
	for (	IDX_TYPE iter_dim_cnt=(IDX_TYPE)0,
			dim_in_base=(IDX_TYPE)0;
			iter_dim_cnt<iter_dim_num;
			iter_dim_cnt+=(IDX_TYPE)1,
			dim_in_base+=(IDX_TYPE)DIM_PARA)
	{
		# pragma HLS LOOP_TRIPCOUNT min=DIM_COUNT max=DIM_COUNT avg=DIM_COUNT
		# pragma HLS PIPELINE II=1 style=frp

		// Whether it is the first or the last iteration.
		BOL_TYPE first = (dim_in_base == (IDX_TYPE)0);
		BOL_TYPE last = (dim_in_base + (IDX_TYPE)DIM_PARA >= (IDX_TYPE)NODE_EMB_DIM);

		// Get biases of all the output dimensions.
		if (first)
		for (IDX_TYPE dim_out=0; dim_out<NODE_EMB_DIM; dim_out++) {
			# pragma HLS UNROLL
			DT1_TYPE bias = g_node_nw_bias_2[dim_out];
			for (IDX_TYPE node_idx_offset=0; node_idx_offset<NODE_UNI_PARA; node_idx_offset++)
				# pragma HLS UNROLL
				linear_out_cache[node_idx_offset][dim_out] = bias; }

		// Process NODE_UNI_PARA nodes each time.
		for (IDX_TYPE node_idx_offset=0; node_idx_offset<NODE_UNI_PARA; node_idx_offset++) {
			# pragma HLS UNROLL

			// Fetch input data.
			DT1_PARA input;
			# pragma HLS AGGREGATE variable=input
			linear_in[node_idx_offset].read(input);

			// Bypass channel.
			DT1_PARA emb;
			# pragma HLS AGGREGATE variable=emb
			node_emb_old_in[node_idx_offset].read(emb);
			node_emb_old_out[node_idx_offset].write(emb);

			// Get the weight and calculate product.
			for (IDX_TYPE dim_out=0; dim_out<NODE_EMB_DIM; dim_out++) {
				# pragma HLS UNROLL
				// Calculate the sum of products.
				DT0_TYPE sum_of_product = (DT0_TYPE)0.0;
				for (IDX_TYPE dim_in_offset=0; dim_in_offset<DIM_PARA; dim_in_offset++) {
					# pragma HLS UNROLL
					IDX_TYPE dim_in = dim_in_base + dim_in_offset;
					IDX_TYPE dim_in_overflow = (dim_in >= (IDX_TYPE)NODE_EMB_DIM);
					DT1_TYPE product_cache = dim_in_overflow? (DT1_TYPE)0.0:
							(DT1_TYPE)(input[dim_in_offset] * g_node_nw_weight_2[dim_out][dim_in]);
					sum_of_product += product_cache; }
				// Get the old data and sum them together.
				DT0_TYPE data_old = linear_out_cache[node_idx_offset][dim_out];
				DT0_TYPE data_new = sum_of_product + data_old;
				linear_out_cache[node_idx_offset][dim_out] = data_new; }}

		// Send the data into the FIFO.
		if (last)
		for (IDX_TYPE node_idx_offset=0; node_idx_offset<NODE_UNI_PARA; node_idx_offset++)
			# pragma HLS UNROLL
			linear_out[node_idx_offset].write(linear_out_cache[node_idx_offset]);
	}
}

static void node_linear_3(
		IDX_TYPE num_of_nodes,
		hls::stream<DT1_PARA> linear_in[NODE_UNI_PARA],
		hls::stream<DT0_VECTOR> linear_out[NODE_UNI_PARA],
		hls::stream<DT1_PARA> node_emb_old_in[NODE_UNI_PARA],
		hls::stream<DT1_PARA> node_emb_old_out[NODE_UNI_PARA])
{
	# pragma HLS INLINE off
	# pragma HLS ARRAY_PARTITION variable=g_node_nw_weight_3 dim=1 type=complete
	# pragma HLS ARRAY_PARTITION variable=g_node_nw_weight_3 dim=2 type=cyclic factor=DIM_PARA
	# pragma HLS ARRAY_PARTITION variable=g_node_nw_bias_3 dim=1 type=complete

	// Temporary array of linear outputs.
	DT0_VECTOR linear_out_cache[NODE_UNI_PARA];
	# pragma HLS AGGREGATE variable=linear_out_cache
	# pragma HLS ARRAY_PARTITION dim=1 type=complete variable=linear_out_cache

	// Process NODE_UNI_PARA nodes each time.
	// Process DIM_PARA dimensions each time.
	IDX_TYPE iter_uni_num = ceildiv(num_of_nodes, (IDX_TYPE)NODE_UNI_PARA);
	IDX_TYPE iter_dim_num = ceildiv((IDX_TYPE)NODE_EMB_DIM, (IDX_TYPE)DIM_PARA);
	for (	IDX_TYPE iter_uni_cnt=(IDX_TYPE)0,
			node_idx_base=(IDX_TYPE)0;
			iter_uni_cnt<iter_uni_num;
			iter_uni_cnt+=(IDX_TYPE)1,
			node_idx_base+=(IDX_TYPE)NODE_UNI_PARA)
	# pragma HLS LOOP_TRIPCOUNT min=NODE_UNI_COUNT max=NODE_UNI_COUNT avg=NODE_UNI_COUNT
	for (	IDX_TYPE iter_dim_cnt=(IDX_TYPE)0,
			dim_in_base=(IDX_TYPE)0;
			iter_dim_cnt<iter_dim_num;
			iter_dim_cnt+=(IDX_TYPE)1,
			dim_in_base+=(IDX_TYPE)DIM_PARA)
	{
		# pragma HLS LOOP_TRIPCOUNT min=DIM_COUNT max=DIM_COUNT avg=DIM_COUNT
		# pragma HLS PIPELINE II=1 style=frp

		// Whether it is the first or the last iteration.
		BOL_TYPE first = (dim_in_base == (IDX_TYPE)0);
		BOL_TYPE last = (dim_in_base + (IDX_TYPE)DIM_PARA >= (IDX_TYPE)NODE_EMB_DIM);

		// Get biases of all the output dimensions.
		if (first)
		for (IDX_TYPE dim_out=0; dim_out<NODE_EMB_DIM; dim_out++) {
			# pragma HLS UNROLL
			DT1_TYPE bias = g_node_nw_bias_3[dim_out];
			for (IDX_TYPE node_idx_offset=0; node_idx_offset<NODE_UNI_PARA; node_idx_offset++)
				# pragma HLS UNROLL
				linear_out_cache[node_idx_offset][dim_out] = bias; }

		// Process NODE_UNI_PARA nodes each time.
		for (IDX_TYPE node_idx_offset=0; node_idx_offset<NODE_UNI_PARA; node_idx_offset++) {
			# pragma HLS UNROLL

			// Fetch input data.
			DT1_PARA input;
			# pragma HLS AGGREGATE variable=input
			linear_in[node_idx_offset].read(input);

			// Bypass channel.
			DT1_PARA emb;
			# pragma HLS AGGREGATE variable=emb
			node_emb_old_in[node_idx_offset].read(emb);
			node_emb_old_out[node_idx_offset].write(emb);

			// Get the weight and calculate product.
			for (IDX_TYPE dim_out=0; dim_out<NODE_EMB_DIM; dim_out++) {
				# pragma HLS UNROLL
				// Calculate the sum of products.
				DT0_TYPE sum_of_product = (DT0_TYPE)0.0;
				for (IDX_TYPE dim_in_offset=0; dim_in_offset<DIM_PARA; dim_in_offset++) {
					# pragma HLS UNROLL
					IDX_TYPE dim_in = dim_in_base + dim_in_offset;
					IDX_TYPE dim_in_overflow = (dim_in >= (IDX_TYPE)NODE_EMB_DIM);
					DT1_TYPE product_cache = dim_in_overflow? (DT1_TYPE)0.0:
							(DT1_TYPE)(input[dim_in_offset] * g_node_nw_weight_3[dim_out][dim_in]);
					sum_of_product += product_cache; }
				// Get the old data and sum them together.
				DT0_TYPE data_old = linear_out_cache[node_idx_offset][dim_out];
				DT0_TYPE data_new = sum_of_product + data_old;
				linear_out_cache[node_idx_offset][dim_out] = data_new; }}

		// Send the data into the FIFO.
		if (last)
		for (IDX_TYPE node_idx_offset=0; node_idx_offset<NODE_UNI_PARA; node_idx_offset++)
			# pragma HLS UNROLL
			linear_out[node_idx_offset].write(linear_out_cache[node_idx_offset]);
	}
}

static void node_activation_0(
		IDX_TYPE num_of_nodes,
		hls::stream<DT0_VECTOR> activation_in[NODE_UNI_PARA],
		hls::stream<DT1_PARA> activation_out[NODE_UNI_PARA],
		hls::stream<DT1_PARA> node_emb_old_in[NODE_UNI_PARA],
		hls::stream<DT1_PARA> node_emb_old_out[NODE_UNI_PARA])
{
	# pragma HLS INLINE off

	// Temporary array of activation inputs.
	DT0_TYPE act_in[NODE_UNI_PARA][NODE_EMB_DIM];
	# pragma HLS ARRAY_PARTITION variable=act_in dim=1 type=complete
	# pragma HLS ARRAY_PARTITION variable=act_in dim=2 type=complete

	// Process NODE_UNI_PARA nodes each time.
	// Process DIM_PARA dimensions each time.
	IDX_TYPE iter_uni_num = ceildiv(num_of_nodes, (IDX_TYPE)NODE_UNI_PARA);
	IDX_TYPE iter_dim_num = ceildiv((IDX_TYPE)NODE_EMB_DIM, (IDX_TYPE)DIM_PARA);
	for (	IDX_TYPE iter_uni_cnt=(IDX_TYPE)0,
			node_idx_base=(IDX_TYPE)0;
			iter_uni_cnt<iter_uni_num;
			iter_uni_cnt+=(IDX_TYPE)1,
			node_idx_base+=(IDX_TYPE)NODE_UNI_PARA)
	# pragma HLS LOOP_TRIPCOUNT min=NODE_UNI_COUNT max=NODE_UNI_COUNT avg=NODE_UNI_COUNT
	for (	IDX_TYPE iter_dim_cnt=(IDX_TYPE)0,
			dim_out_base=(IDX_TYPE)0;
			iter_dim_cnt<iter_dim_num;
			iter_dim_cnt+=(IDX_TYPE)1,
			dim_out_base+=(IDX_TYPE)DIM_PARA)
	{
		# pragma HLS LOOP_TRIPCOUNT min=DIM_COUNT max=DIM_COUNT avg=DIM_COUNT
		# pragma HLS PIPELINE II=1 style=frp

		// Whether it is the first iteration.
		BOL_TYPE first = (dim_out_base == (IDX_TYPE)0);

		// Read new data from the last linear layer.
		if (first)
		for (IDX_TYPE node_idx_offset=0; node_idx_offset<NODE_UNI_PARA; node_idx_offset++) {
			# pragma HLS UNROLL
			DT0_VECTOR act_in_cache;
			# pragma HLS AGGREGATE variable=act_in_cache
			activation_in[node_idx_offset].read(act_in_cache);
			for (IDX_TYPE dim_in=0; dim_in<NODE_EMB_DIM; dim_in++)
				# pragma HLS UNROLL
				act_in[node_idx_offset][dim_in] = act_in_cache[dim_in]; }

		// Pass through the activation function.
		for (IDX_TYPE node_idx_offset=0; node_idx_offset<NODE_UNI_PARA; node_idx_offset++) {
			# pragma HLS UNROLL
			// Temporary array of activation outputs.
			DT1_PARA act_out;
			# pragma HLS AGGREGATE variable=act_out
			// Process DIM_PARA dimensions each time.
			for (IDX_TYPE dim_out_offset=0; dim_out_offset<DIM_PARA; dim_out_offset++) {
				# pragma HLS UNROLL
				// Get the index of the output dimension.
				IDX_TYPE dim_out = dim_out_base + dim_out_offset;
				BOL_TYPE dim_out_overflow = (dim_out >= (IDX_TYPE)NODE_EMB_DIM);
				// Pass through the activation function.
				# ifdef ACT_TANH
				act_out[dim_out_offset] = dim_out_overflow? (DT1_TYPE)0.0:
						act_tanh_node_nw((IDX_TYPE)0, node_idx_offset, dim_out_offset, act_in[node_idx_offset][dim_out]); }
				# endif
				# ifdef ACT_SIGM
				act_out[dim_out_offset] = dim_out_overflow? (DT1_TYPE)0.0:
						act_sigm_node_nw((IDX_TYPE)0, node_idx_offset, dim_out_offset, act_in[node_idx_offset][dim_out]); }
				# endif
			// Write the activation outputs into the FIFO.
			activation_out[node_idx_offset].write(act_out); }

		// Bypass channel.
		for (IDX_TYPE node_idx_offset=0; node_idx_offset<NODE_UNI_PARA; node_idx_offset++) {
			# pragma HLS UNROLL
			DT1_PARA cache;
			# pragma HLS AGGREGATE variable=cache
			node_emb_old_in[node_idx_offset].read(cache);
			node_emb_old_out[node_idx_offset].write(cache); }
	}
}

static void node_activation_1(
		IDX_TYPE num_of_nodes,
		hls::stream<DT0_VECTOR> activation_in[NODE_UNI_PARA],
		hls::stream<DT1_PARA> activation_out[NODE_UNI_PARA],
		hls::stream<DT1_PARA> node_emb_old_in[NODE_UNI_PARA],
		hls::stream<DT1_PARA> node_emb_old_out[NODE_UNI_PARA])
{
	# pragma HLS INLINE off

	// Temporary array of activation inputs.
	DT0_TYPE act_in[NODE_UNI_PARA][NODE_EMB_DIM];
	# pragma HLS ARRAY_PARTITION variable=act_in dim=1 type=complete
	# pragma HLS ARRAY_PARTITION variable=act_in dim=2 type=complete

	// Process NODE_UNI_PARA nodes each time.
	// Process DIM_PARA dimensions each time.
	IDX_TYPE iter_uni_num = ceildiv(num_of_nodes, (IDX_TYPE)NODE_UNI_PARA);
	IDX_TYPE iter_dim_num = ceildiv((IDX_TYPE)NODE_EMB_DIM, (IDX_TYPE)DIM_PARA);
	for (	IDX_TYPE iter_uni_cnt=(IDX_TYPE)0,
			node_idx_base=(IDX_TYPE)0;
			iter_uni_cnt<iter_uni_num;
			iter_uni_cnt+=(IDX_TYPE)1,
			node_idx_base+=(IDX_TYPE)NODE_UNI_PARA)
	# pragma HLS LOOP_TRIPCOUNT min=NODE_UNI_COUNT max=NODE_UNI_COUNT avg=NODE_UNI_COUNT
	for (	IDX_TYPE iter_dim_cnt=(IDX_TYPE)0,
			dim_out_base=(IDX_TYPE)0;
			iter_dim_cnt<iter_dim_num;
			iter_dim_cnt+=(IDX_TYPE)1,
			dim_out_base+=(IDX_TYPE)DIM_PARA)
	{
		# pragma HLS LOOP_TRIPCOUNT min=DIM_COUNT max=DIM_COUNT avg=DIM_COUNT
		# pragma HLS PIPELINE II=1 style=frp

		// Whether it is the first iteration.
		BOL_TYPE first = (dim_out_base == (IDX_TYPE)0);

		// Read new data from the last linear layer.
		if (first)
		for (IDX_TYPE node_idx_offset=0; node_idx_offset<NODE_UNI_PARA; node_idx_offset++) {
			# pragma HLS UNROLL
			DT0_VECTOR act_in_cache;
			# pragma HLS AGGREGATE variable=act_in_cache
			activation_in[node_idx_offset].read(act_in_cache);
			for (IDX_TYPE dim_in=0; dim_in<NODE_EMB_DIM; dim_in++)
				# pragma HLS UNROLL
				act_in[node_idx_offset][dim_in] = act_in_cache[dim_in]; }

		// Pass through the activation function.
		for (IDX_TYPE node_idx_offset=0; node_idx_offset<NODE_UNI_PARA; node_idx_offset++) {
			# pragma HLS UNROLL
			// Temporary array of activation outputs.
			DT1_PARA act_out;
			# pragma HLS AGGREGATE variable=act_out
			// Process DIM_PARA dimensions each time.
			for (IDX_TYPE dim_out_offset=0; dim_out_offset<DIM_PARA; dim_out_offset++) {
				# pragma HLS UNROLL
				// Get the index of the output dimension.
				IDX_TYPE dim_out = dim_out_base + dim_out_offset;
				BOL_TYPE dim_out_overflow = (dim_out >= (IDX_TYPE)NODE_EMB_DIM);
				// Pass through the activation function.
				# ifdef ACT_TANH
				act_out[dim_out_offset] = dim_out_overflow? (DT1_TYPE)0.0:
						act_tanh_node_nw((IDX_TYPE)1, node_idx_offset, dim_out_offset, act_in[node_idx_offset][dim_out]); }
				# endif
				# ifdef ACT_SIGM
				act_out[dim_out_offset] = dim_out_overflow? (DT1_TYPE)0.0:
						act_sigm_node_nw((IDX_TYPE)1, node_idx_offset, dim_out_offset, act_in[node_idx_offset][dim_out]); }
				# endif
			// Write the activation outputs into the FIFO.
			activation_out[node_idx_offset].write(act_out); }

		// Bypass channel.
		for (IDX_TYPE node_idx_offset=0; node_idx_offset<NODE_UNI_PARA; node_idx_offset++) {
			# pragma HLS UNROLL
			DT1_PARA cache;
			# pragma HLS AGGREGATE variable=cache
			node_emb_old_in[node_idx_offset].read(cache);
			node_emb_old_out[node_idx_offset].write(cache); }
	}
}

static void node_activation_2(
		IDX_TYPE num_of_nodes,
		hls::stream<DT0_VECTOR> activation_in[NODE_UNI_PARA],
		hls::stream<DT1_PARA> activation_out[NODE_UNI_PARA],
		hls::stream<DT1_PARA> node_emb_old_in[NODE_UNI_PARA],
		hls::stream<DT1_PARA> node_emb_old_out[NODE_UNI_PARA])
{
	# pragma HLS INLINE off

	// Temporary array of activation inputs.
	DT0_TYPE act_in[NODE_UNI_PARA][NODE_EMB_DIM];
	# pragma HLS ARRAY_PARTITION variable=act_in dim=1 type=complete
	# pragma HLS ARRAY_PARTITION variable=act_in dim=2 type=complete

	// Process NODE_UNI_PARA nodes each time.
	// Process DIM_PARA dimensions each time.
	IDX_TYPE iter_uni_num = ceildiv(num_of_nodes, (IDX_TYPE)NODE_UNI_PARA);
	IDX_TYPE iter_dim_num = ceildiv((IDX_TYPE)NODE_EMB_DIM, (IDX_TYPE)DIM_PARA);
	for (	IDX_TYPE iter_uni_cnt=(IDX_TYPE)0,
			node_idx_base=(IDX_TYPE)0;
			iter_uni_cnt<iter_uni_num;
			iter_uni_cnt+=(IDX_TYPE)1,
			node_idx_base+=(IDX_TYPE)NODE_UNI_PARA)
	# pragma HLS LOOP_TRIPCOUNT min=NODE_UNI_COUNT max=NODE_UNI_COUNT avg=NODE_UNI_COUNT
	for (	IDX_TYPE iter_dim_cnt=(IDX_TYPE)0,
			dim_out_base=(IDX_TYPE)0;
			iter_dim_cnt<iter_dim_num;
			iter_dim_cnt+=(IDX_TYPE)1,
			dim_out_base+=(IDX_TYPE)DIM_PARA)
	{
		# pragma HLS LOOP_TRIPCOUNT min=DIM_COUNT max=DIM_COUNT avg=DIM_COUNT
		# pragma HLS PIPELINE II=1 style=frp

		// Whether it is the first iteration.
		BOL_TYPE first = (dim_out_base == (IDX_TYPE)0);

		// Read new data from the last linear layer.
		if (first)
		for (IDX_TYPE node_idx_offset=0; node_idx_offset<NODE_UNI_PARA; node_idx_offset++) {
			# pragma HLS UNROLL
			DT0_VECTOR act_in_cache;
			# pragma HLS AGGREGATE variable=act_in_cache
			activation_in[node_idx_offset].read(act_in_cache);
			for (IDX_TYPE dim_in=0; dim_in<NODE_EMB_DIM; dim_in++)
				# pragma HLS UNROLL
				act_in[node_idx_offset][dim_in] = act_in_cache[dim_in]; }

		// Pass through the activation function.
		for (IDX_TYPE node_idx_offset=0; node_idx_offset<NODE_UNI_PARA; node_idx_offset++) {
			# pragma HLS UNROLL
			// Temporary array of activation outputs.
			DT1_PARA act_out;
			# pragma HLS AGGREGATE variable=act_out
			// Process DIM_PARA dimensions each time.
			for (IDX_TYPE dim_out_offset=0; dim_out_offset<DIM_PARA; dim_out_offset++) {
				# pragma HLS UNROLL
				// Get the index of the output dimension.
				IDX_TYPE dim_out = dim_out_base + dim_out_offset;
				BOL_TYPE dim_out_overflow = (dim_out >= (IDX_TYPE)NODE_EMB_DIM);
				// Pass through the activation function.
				# ifdef ACT_TANH
				act_out[dim_out_offset] = dim_out_overflow? (DT1_TYPE)0.0:
						act_tanh_node_nw((IDX_TYPE)2, node_idx_offset, dim_out_offset, act_in[node_idx_offset][dim_out]); }
				# endif
				# ifdef ACT_SIGM
				act_out[dim_out_offset] = dim_out_overflow? (DT1_TYPE)0.0:
						act_sigm_node_nw((IDX_TYPE)2, node_idx_offset, dim_out_offset, act_in[node_idx_offset][dim_out]); }
				# endif
			// Write the activation outputs into the FIFO.
			activation_out[node_idx_offset].write(act_out); }

		// Bypass channel.
		for (IDX_TYPE node_idx_offset=0; node_idx_offset<NODE_UNI_PARA; node_idx_offset++) {
			# pragma HLS UNROLL
			DT1_PARA cache;
			# pragma HLS AGGREGATE variable=cache
			node_emb_old_in[node_idx_offset].read(cache);
			node_emb_old_out[node_idx_offset].write(cache); }
	}
}

static void node_activation_3(
		BOL_TYPE last_layer,
		IDX_TYPE num_of_nodes,
		hls::stream<DT0_VECTOR> activation_in[NODE_UNI_PARA],
		hls::stream<DT1_PARA> node_emb_old_in[NODE_UNI_PARA],
		hls::stream<DT1_PARA> mlp_out_0[NODE_UNI_PARA],
		hls::stream<DT1_PARA> mlp_out_1[NODE_UNI_PARA])
{
	# pragma HLS INLINE off

	// Temporary array of activation inputs.
	DT0_TYPE act_in[NODE_UNI_PARA][NODE_EMB_DIM];
	# pragma HLS ARRAY_PARTITION variable=act_in dim=1 type=complete
	# pragma HLS ARRAY_PARTITION variable=act_in dim=2 type=complete

	// Process NODE_UNI_PARA nodes each time.
	// Process DIM_PARA dimensions each time.
	IDX_TYPE iter_uni_num = ceildiv(num_of_nodes, (IDX_TYPE)NODE_UNI_PARA);
	IDX_TYPE iter_dim_num = ceildiv((IDX_TYPE)NODE_EMB_DIM, (IDX_TYPE)DIM_PARA);
	for (	IDX_TYPE iter_uni_cnt=(IDX_TYPE)0,
			node_idx_base=(IDX_TYPE)0;
			iter_uni_cnt<iter_uni_num;
			iter_uni_cnt+=(IDX_TYPE)1,
			node_idx_base+=(IDX_TYPE)NODE_UNI_PARA)
	# pragma HLS LOOP_TRIPCOUNT min=NODE_UNI_COUNT max=NODE_UNI_COUNT avg=NODE_UNI_COUNT
	for (	IDX_TYPE iter_dim_cnt=(IDX_TYPE)0,
			dim_out_base=(IDX_TYPE)0;
			iter_dim_cnt<iter_dim_num;
			iter_dim_cnt+=(IDX_TYPE)1,
			dim_out_base+=(IDX_TYPE)DIM_PARA)
	{
		# pragma HLS LOOP_TRIPCOUNT min=DIM_COUNT max=DIM_COUNT avg=DIM_COUNT
		# pragma HLS PIPELINE II=1 style=frp

		// Whether it is the first iteration.
		BOL_TYPE first = (dim_out_base == (IDX_TYPE)0);

		// Read new data from the last linear layer.
		if (first)
		for (IDX_TYPE node_idx_offset=0; node_idx_offset<NODE_UNI_PARA; node_idx_offset++) {
			# pragma HLS UNROLL
			DT0_VECTOR act_in_cache;
			# pragma HLS AGGREGATE variable=act_in_cache
			activation_in[node_idx_offset].read(act_in_cache);
			for (IDX_TYPE dim_in=0; dim_in<NODE_EMB_DIM; dim_in++)
				# pragma HLS UNROLL
				act_in[node_idx_offset][dim_in] = act_in_cache[dim_in]; }

		// Output results of node network.
		for (IDX_TYPE node_idx_offset=0; node_idx_offset<NODE_UNI_PARA; node_idx_offset++) {
			# pragma HLS UNROLL
			// Temporary array of activation outputs.
			DT1_PARA act_out;
			# pragma HLS AGGREGATE variable=act_out
			// Process DIM_PARA dimensions each time.
			for (IDX_TYPE dim_out_offset=0; dim_out_offset<DIM_PARA; dim_out_offset++) {
				# pragma HLS UNROLL
				// Get the index of the output dimension.
				IDX_TYPE dim_out = dim_out_base + dim_out_offset;
				BOL_TYPE dim_out_overflow = (dim_out >= (IDX_TYPE)NODE_EMB_DIM);
				// Pass through the activation function.
				# ifdef ACT_TANH
				act_out[dim_out_offset] = dim_out_overflow? (DT1_TYPE)0.0:
						act_tanh_node_nw((IDX_TYPE)3, node_idx_offset, dim_out_offset, act_in[node_idx_offset][dim_out]); }
				# endif
				# ifdef ACT_SIGM
				act_out[dim_out_offset] = dim_out_overflow? (DT1_TYPE)0.0:
						act_sigm_node_nw((IDX_TYPE)3, node_idx_offset, dim_out_offset, act_in[node_idx_offset][dim_out]); }
				# endif

			// Sum of the old embedding and the output data of the activation function.
			DT1_PARA emb_old;
			DT1_PARA emb_new;
			# pragma HLS AGGREGATE variable=emb_old
			# pragma HLS AGGREGATE variable=emb_new
			node_emb_old_in[node_idx_offset].read(emb_old);
			for (IDX_TYPE dim_out_offset=0; dim_out_offset<DIM_PARA; dim_out_offset++)
				# pragma HLS UNROLL
				emb_new[dim_out_offset] = act_out[dim_out_offset] + emb_old[dim_out_offset];

			// Write it into FIFOs.
			if (!last_layer)
			mlp_out_1[node_idx_offset].write(emb_new);
			mlp_out_0[node_idx_offset].write(emb_new); }
	}
}

static void node_nw(
		BOL_TYPE last_layer,
		IDX_TYPE num_of_nodes,
		DT1_TYPE node_emb[NODE_NUM_MAX][NODE_EMB_DIM],
		DT0_TYPE msg_dst[2][EDGE_UNI_PARA][NODE_NUM_MAX][NODE_EMB_DIM],
		DT0_TYPE msg_src[2][EDGE_UNI_PARA][NODE_NUM_MAX][NODE_EMB_DIM],
		IDX_TYPE msg_dst_flag[EDGE_UNI_PARA][NODE_NUM_MAX],
		IDX_TYPE msg_src_flag[EDGE_UNI_PARA][NODE_NUM_MAX],
		hls::stream<DT1_PARA> mlp_out_0[NODE_UNI_PARA],
		hls::stream<DT1_PARA> mlp_out_1[NODE_UNI_PARA])
{
	# pragma HLS INLINE off
	# pragma HLS DATAFLOW

	// FIFO declarations.
	hls::stream<DT0_VECTOR> linear_out_0[NODE_UNI_PARA];
	hls::stream<DT0_VECTOR> linear_out_1[NODE_UNI_PARA];
	hls::stream<DT0_VECTOR> linear_out_2[NODE_UNI_PARA];
	hls::stream<DT0_VECTOR> linear_out_3[NODE_UNI_PARA];
	# pragma HLS STREAM variable=linear_out_0 depth=FIFO_DEPTH
	# pragma HLS STREAM variable=linear_out_1 depth=FIFO_DEPTH
	# pragma HLS STREAM variable=linear_out_2 depth=FIFO_DEPTH
	# pragma HLS STREAM variable=linear_out_3 depth=FIFO_DEPTH
	hls::stream<DT1_PARA> linear_emb_out_0[NODE_UNI_PARA];
	hls::stream<DT1_PARA> linear_emb_out_1[NODE_UNI_PARA];
	hls::stream<DT1_PARA> linear_emb_out_2[NODE_UNI_PARA];
	hls::stream<DT1_PARA> linear_emb_out_3[NODE_UNI_PARA];
	# pragma HLS STREAM variable=linear_emb_out_0 depth=(FIFO_DEPTH*DIM_COUNT)
	# pragma HLS STREAM variable=linear_emb_out_1 depth=(FIFO_DEPTH*DIM_COUNT)
	# pragma HLS STREAM variable=linear_emb_out_2 depth=(FIFO_DEPTH*DIM_COUNT)
	# pragma HLS STREAM variable=linear_emb_out_3 depth=(FIFO_DEPTH*DIM_COUNT)
	hls::stream<DT1_PARA> activation_out_0[NODE_UNI_PARA];
	hls::stream<DT1_PARA> activation_out_1[NODE_UNI_PARA];
	hls::stream<DT1_PARA> activation_out_2[NODE_UNI_PARA];
	# pragma HLS STREAM variable=activation_out_0 depth=FIFO_DEPTH
	# pragma HLS STREAM variable=activation_out_1 depth=FIFO_DEPTH
	# pragma HLS STREAM variable=activation_out_2 depth=FIFO_DEPTH
	hls::stream<DT1_PARA> activation_emb_out_0[NODE_UNI_PARA];
	hls::stream<DT1_PARA> activation_emb_out_1[NODE_UNI_PARA];
	hls::stream<DT1_PARA> activation_emb_out_2[NODE_UNI_PARA];
	# pragma HLS STREAM variable=activation_emb_out_0 depth=FIFO_DEPTH
	# pragma HLS STREAM variable=activation_emb_out_1 depth=FIFO_DEPTH
	# pragma HLS STREAM variable=activation_emb_out_2 depth=FIFO_DEPTH

	// Invoke functions.
	node_linear_0(
			num_of_nodes,
			msg_dst,
			msg_src,
			msg_dst_flag,
			msg_src_flag,
			node_emb,
			linear_out_0,
			linear_emb_out_0);
	node_activation_0(
			num_of_nodes,
			linear_out_0,
			activation_out_0,
			linear_emb_out_0,
			activation_emb_out_0);
	node_linear_1(
			num_of_nodes,
			activation_out_0,
			linear_out_1,
			activation_emb_out_0,
			linear_emb_out_1);
	node_activation_1(
			num_of_nodes,
			linear_out_1,
			activation_out_1,
			linear_emb_out_1,
			activation_emb_out_1);
	node_linear_2(
			num_of_nodes,
			activation_out_1,
			linear_out_2,
			activation_emb_out_1,
			linear_emb_out_2);
	node_activation_2(
			num_of_nodes,
			linear_out_2,
			activation_out_2,
			linear_emb_out_2,
			activation_emb_out_2);
	node_linear_3(
			num_of_nodes,
			activation_out_2,
			linear_out_3,
			activation_emb_out_2,
			linear_emb_out_3);
	node_activation_3(
			last_layer,
			num_of_nodes,
			linear_out_3,
			linear_emb_out_3,
			mlp_out_0,
			mlp_out_1);
}



// **************************** Adapter **************************** //

/*

	This part is to receive data from the node-processing network and
	send data to the edge-processing network. Once receive new node
	embeddings, upload them into FIFOs. Once the edge-processing net-
	work and FIFOs are ready, download them from FIFOs and send them
	into the edge-processing network.

	Version: 0.1
		Author: Hanqing Zhang
		Email : hanqing.zhang@zju.edu.cn
		Date  : July-2-2023
		Commit: Add function "refresh_array" and "receive_n_upload".

	Version: 0.2
		Author: Hanqing Zhang
		Email : hanqing.zhang@zju.edu.cn
		Date  : July-3-2023
		Commit: Merge function "refresh_array" into "receive_n_upload".

	Version: 0.3
		Author: Hanqing Zhang
		Email : hanqing.zhang@zju.edu.cn
		Date  : July-4-2023
		Commit: Add function "download_n_send" and "adapter".

	Version: 0.4
		Author: Hanqing Zhang
		Email : hanqing.zhang@zju.edu.cn
		Date  : July-19-2023
		Commit: Combine msg_out and msg_in together to reduce the area.

	Version: 0.5
		Author: Hanqing Zhang
		Email : hanqing.zhang@zju.edu.cn
		Date  : July-21-2023
		Commit: Fixed bugs of not being fully pipelined.

	Version: 0.6
		Author: Hanqing Zhang
		Email : hanqing.zhang@zju.edu.cn
		Date  : August-5-2023
		Commit: Implement FIFOs using shift registers to reduce the
				utilization of BRAMs.

	Version: 0.7
		Author: Hanqing Zhang
		Email : hanqing.zhang@zju.edu.cn
		Date  : August-16-2023
		Commit: Merge the functions "receive_n_upload" and
				"download_n_send" to reduce the usage of BRAMs.

	Version: 0.8
		Author: Hanqing Zhang
		Email : hanqing.zhang@zju.edu.cn
		Date  : October-1-2023
		Commit: Modify for dimension parallelism.

	Function:
		void adapter(
			IDX_TYPE num_of_nodes,
			IDX_TYPE num_of_edges,
			hls::stream<DT1_PARA> node_fifo[NODE_UNI_PARA],
			hls::stream<EDGE_TYPE> edge_fifo[EDGE_UNI_PARA]);
		Input arguments:
			num_of_nodes: the number of nodes.
			num_of_edges: the number of edges.
			node_fifo: new node embeddings from the node-processing network.
		Output arguments:
			edge_fifo: new node embeddings used to calculate edge labels.

*/

static void adapter(
		IDX_TYPE num_of_nodes,
		IDX_TYPE num_of_edges,
		hls::stream<DT1_PARA> node_fifo[NODE_UNI_PARA],
		hls::stream<EDGE_TYPE> edge_fifo[EDGE_UNI_PARA])
{
	# pragma HLS INLINE off

	// Temporary variables of edge indexes.
	IDX_TYPE edge_idx_base;
	IDX_TYPE edge_idx_base_trans = (IDX_TYPE)0;

	// Temporary arrays of node embeddings.
	DT1_TYPE node_emb_src[EDGE_UNI_PARA][NODE_NUM_MAX][NODE_EMB_DIM];
	DT1_TYPE node_emb_dst[EDGE_UNI_PARA][NODE_NUM_MAX][NODE_EMB_DIM];
	# pragma HLS ARRAY_PARTITION dim=1 type=complete variable=node_emb_src
	# pragma HLS ARRAY_PARTITION dim=1 type=complete variable=node_emb_dst
	# pragma HLS ARRAY_PARTITION dim=2 type=cyclic factor=NODE_UNI_PARA variable=node_emb_src
	# pragma HLS ARRAY_PARTITION dim=2 type=cyclic factor=NODE_UNI_PARA variable=node_emb_dst
	# pragma HLS ARRAY_PARTITION dim=3 type=cyclic factor=DIM_PARA variable=node_emb_src
	# pragma HLS ARRAY_PARTITION dim=3 type=cyclic factor=DIM_PARA variable=node_emb_dst

	// Phase 1: receive from node FIFOs, and send to edge FIFOs.
	int iter_num_0 = DIM_COUNT * ceildiv((int)num_of_nodes, NODE_UNI_PARA);
	for (int iter_cnt_0=0; iter_cnt_0<iter_num_0; iter_cnt_0++)
	{
		# pragma HLS LOOP_TRIPCOUNT min=(DIM_COUNT*NODE_UNI_COUNT/2) max=(DIM_COUNT*NODE_UNI_COUNT/2) avg=(DIM_COUNT*NODE_UNI_COUNT/2)
		# pragma HLS PIPELINE II=1 style=frp

		// Calculate node index and dimension.
		IDX_TYPE node_idx_base = (IDX_TYPE)(iter_cnt_0 / DIM_COUNT * NODE_UNI_PARA);
		IDX_TYPE dim_base = (IDX_TYPE)(iter_cnt_0 % DIM_COUNT * DIM_PARA);

		// Receive data from node FIFOs.
		for (IDX_TYPE node_idx_offset=0; node_idx_offset<NODE_UNI_PARA; node_idx_offset++) {
			# pragma HLS UNROLL
			// Get the node index.
			IDX_TYPE node_idx = node_idx_base + node_idx_offset;
			BOL_TYPE node_idx_overflow = (node_idx >= num_of_nodes);
			// Fetch the new embedding.
			DT1_PARA emb;
			# pragma HLS AGGREGATE variable=emb
			node_fifo[node_idx_offset].read(emb);
			// Save the new embedding in arrays.
			if (!node_idx_overflow)
			for (IDX_TYPE edge_idx_offset=0; edge_idx_offset<EDGE_UNI_PARA; edge_idx_offset++)
			# pragma HLS UNROLL
			for (IDX_TYPE dim_offset=0; dim_offset<DIM_PARA; dim_offset++) {
				# pragma HLS UNROLL
				IDX_TYPE dim = dim_base + dim_offset;
				BOL_TYPE dim_overflow = (dim >= (IDX_TYPE)NODE_EMB_DIM);
				if (!dim_overflow) {
					node_emb_src[edge_idx_offset][node_idx][dim] = emb[dim_offset];
					node_emb_dst[edge_idx_offset][node_idx][dim] = emb[dim_offset]; }}}

		// Whether the embedding can be output.
		if (dim_base == (IDX_TYPE)0) {
			edge_idx_base = g_determine[node_idx_base];
			if (edge_idx_base != (IDX_TYPE)-1)
			edge_idx_base_trans += (IDX_TYPE)EDGE_UNI_PARA; }

		// Send data to edge FIFOs.
		if (edge_idx_base != (IDX_TYPE)-1)
		for (IDX_TYPE edge_idx_offset=0; edge_idx_offset<EDGE_UNI_PARA; edge_idx_offset++) {
			# pragma HLS UNROLL
			// Get the edge index.
			IDX_TYPE edge_idx = edge_idx_base + edge_idx_offset;
			BOL_TYPE edge_idx_overflow = (edge_idx >= num_of_edges);
			// Create the edge-stream structure.
			EDGE_TYPE edge_cache;
			# pragma HLS AGGREGATE variable=edge_cache
			IDX_TYPE node_idx_src = edge_idx_overflow? (IDX_TYPE)-1:
					g_adj_list_adapter[edge_idx_offset][edge_idx][0];
			IDX_TYPE node_idx_dst = edge_idx_overflow? (IDX_TYPE)-1:
					g_adj_list_adapter[edge_idx_offset][edge_idx][1];
			for (IDX_TYPE dim_offset=0; dim_offset<DIM_PARA; dim_offset++) {
				# pragma HLS UNROLL
				IDX_TYPE dim = dim_base + dim_offset;
				BOL_TYPE dim_overflow = (dim >= (IDX_TYPE)NODE_EMB_DIM);
				edge_cache.node_emb_src[dim_offset] = (dim_overflow || edge_idx_overflow)?
						(DT1_TYPE)0.0: node_emb_src[edge_idx_offset][node_idx_src][dim];
				edge_cache.node_emb_dst[dim_offset] = (dim_overflow || edge_idx_overflow)?
						(DT1_TYPE)0.0: node_emb_dst[edge_idx_offset][node_idx_dst][dim]; }
			// Send the structure out.
			edge_fifo[edge_idx_offset].write(edge_cache); }
	}

	// Phase 2: send to edge FIFOs.
	int iter_num_1 = DIM_COUNT * ceildiv((int)(num_of_edges - edge_idx_base_trans), EDGE_UNI_PARA);
	for (int iter_cnt_1=0; iter_cnt_1<iter_num_1; iter_cnt_1++)
	{
		# pragma HLS LOOP_TRIPCOUNT min=(DIM_COUNT*EDGE_UNI_COUNT/2) max=(DIM_COUNT*EDGE_UNI_COUNT/2) avg=(DIM_COUNT*EDGE_UNI_COUNT/2)
		# pragma HLS PIPELINE II=1 style=frp

		// Calculate edge index and dimension.
		IDX_TYPE edge_idx_base = (IDX_TYPE)(iter_cnt_1 / DIM_COUNT * EDGE_UNI_PARA + (int)edge_idx_base_trans);
		IDX_TYPE dim_base = (IDX_TYPE)(iter_cnt_1 % DIM_COUNT * DIM_PARA);

		// Send data to edge FIFOs.
		for (IDX_TYPE edge_idx_offset=0; edge_idx_offset<EDGE_UNI_PARA; edge_idx_offset++) {
			# pragma HLS UNROLL
			// Get the edge index.
			IDX_TYPE edge_idx = edge_idx_base + edge_idx_offset;
			BOL_TYPE edge_idx_overflow = (edge_idx >= num_of_edges);
			// Create the edge-stream structure.
			EDGE_TYPE edge_cache;
			# pragma HLS AGGREGATE variable=edge_cache
			IDX_TYPE node_idx_src = edge_idx_overflow? (IDX_TYPE)-1:
					g_adj_list_adapter[edge_idx_offset][edge_idx][0];
			IDX_TYPE node_idx_dst = edge_idx_overflow? (IDX_TYPE)-1:
					g_adj_list_adapter[edge_idx_offset][edge_idx][1];
			for (IDX_TYPE dim_offset=0; dim_offset<DIM_PARA; dim_offset++) {
				IDX_TYPE dim = dim_base + dim_offset;
				BOL_TYPE dim_overflow = (dim >= (IDX_TYPE)NODE_EMB_DIM);
				edge_cache.node_emb_src[dim_offset] = (dim_overflow || edge_idx_overflow)?
						(DT1_TYPE)0.0: node_emb_src[edge_idx_offset][node_idx_src][dim];
				edge_cache.node_emb_dst[dim_offset] = (dim_overflow || edge_idx_overflow)?
						(DT1_TYPE)0.0: node_emb_dst[edge_idx_offset][node_idx_dst][dim]; }
			// Send the structure out.
			edge_fifo[edge_idx_offset].write(edge_cache); }
	}
}



// **************************** Edge-Processing Network **************************** //

/*

	This part is to implement edge-processing network. It contains
	4 linear layers and 4 activation layers. For processing a edge,
	the linear layer first fetches one input dimension, products it
	with corresponding weights, and sends them into all the output
	dimensions parallel. After iterating all the input dimension,
	the linear layer send all the output dimensions to the activation
	layer parallel. Once the activation layer receive the data. It
	will divide all the dimensions separately. Then, it'll make each
	data pass through the activation function, and send it into the
	next layer serially.

	Version: 0.1
		Author: Hanqing Zhang
		Email : hanqing.zhang@zju.edu.cn
		Date  : June-30-2023
		Commit: Add function "edge_nw".

	Version: 0.2
		Author: Hanqing Zhang
		Email : hanqing.zhang@zju.edu.cn
		Date  : July-2-2023
		Commit: Modify the output-data type.

	Version: 0.3
		Author: Hanqing Zhang
		Email : hanqing.zhang@zju.edu.cn
		Data  : July-17-2023
		Commit: Modify for the DATAFLOW pragma.

	Version: 0.4
		Author: Hanqing Zhang
		Email : hanqing.zhang@zju.edu.cn
		Data  : September-20-2023
		Commit: Modify for the dimension parallel.

	Version: 0.5
		Author: Hanqing Zhang
		Email : hanqing.zhang@zju.edu.cn
		Data  : October-10-2023
		Commit: Fix the bug that the last layer cannot be activated.

	Version: 0.6
		Author: Hanqing Zhang
		Email : hanqing.zhang@zju.edu.cn
		Data  : October-19-2023
		Commit: Modify the depth of the FIFOs.

	Function:
		void edge_nw(
			BOL_TYPE last_layer,
			IDX_TYPE num_of_edges,
			hls::stream<EDGE_TYPE> edge_stream[EDGE_UNI_PARA],
			hls::stream<DT1_TYPE> value_stream[EDGE_UNI_PARA],
			hls::stream<EDGE_TYPE> msg_stream[EDGE_UNI_PARA]);
		Input arguments:
			last_layer: whether it is the last layer.
			num_of_edges: the number of edges.
			edge_stream: new node embedding stream from the adapter.
		Output arguments:
			value_stream: attention value stream.
			msg_stream: new node embedding stream to the message scatter.

*/

static void edge_linear_0(
		BOL_TYPE last_layer,
		IDX_TYPE num_of_edges,
		hls::stream<DT0_VECTOR> linear_out[EDGE_UNI_PARA],
		hls::stream<EDGE_TYPE> edge_in[EDGE_UNI_PARA],
		hls::stream<EDGE_TYPE> edge_out[EDGE_UNI_PARA])
{
	# pragma HLS INLINE off
	# pragma HLS ARRAY_PARTITION variable=g_edge_nw_weight_0_0 dim=1 type=complete
	# pragma HLS ARRAY_PARTITION variable=g_edge_nw_weight_0_1 dim=1 type=complete
	# pragma HLS ARRAY_PARTITION variable=g_edge_nw_weight_0_0 dim=2 type=cyclic factor=DIM_PARA
	# pragma HLS ARRAY_PARTITION variable=g_edge_nw_weight_0_1 dim=2 type=cyclic factor=DIM_PARA
	# pragma HLS ARRAY_PARTITION variable=g_edge_nw_bias_0 dim=1 type=complete

	// Temporary array of linear outputs.
	DT0_VECTOR linear_out_cache[EDGE_UNI_PARA];
	# pragma HLS AGGREGATE variable=linear_out_cache
	# pragma HLS ARRAY_PARTITION variable=linear_out_cache dim=1 type=complete

	// Process EDGE_UNI_PARA edges each time.
	// Process DIM_PARA dimensions each time.
	IDX_TYPE iter_uni_num = ceildiv(num_of_edges, (IDX_TYPE)EDGE_UNI_PARA);
	IDX_TYPE iter_dim_num = ceildiv((IDX_TYPE)NODE_EMB_DIM, (IDX_TYPE)DIM_PARA);
	for (	IDX_TYPE iter_uni_cnt=(IDX_TYPE)0,
			edge_idx_base=(IDX_TYPE)0;
			iter_uni_cnt<iter_uni_num;
			iter_uni_cnt+=(IDX_TYPE)1,
			edge_idx_base+=(IDX_TYPE)EDGE_UNI_PARA)
	# pragma HLS LOOP_TRIPCOUNT min=EDGE_UNI_COUNT max=EDGE_UNI_COUNT avg=EDGE_UNI_COUNT
	for (	IDX_TYPE iter_dim_cnt=(IDX_TYPE)0,
			dim_in_base=(IDX_TYPE)0;
			iter_dim_cnt<iter_dim_num;
			iter_dim_cnt+=(IDX_TYPE)1,
			dim_in_base+=(IDX_TYPE)DIM_PARA)
	{
		# pragma HLS LOOP_TRIPCOUNT min=DIM_COUNT max=DIM_COUNT avg=DIM_COUNT
		# pragma HLS PIPELINE II=1 style=frp

		// Whether it is the first or the last iteration.
		BOL_TYPE first = (dim_in_base == (IDX_TYPE)0);
		BOL_TYPE last = (dim_in_base + (IDX_TYPE)DIM_PARA >= (IDX_TYPE)NODE_EMB_DIM);
		
		// Get biases of all the output dimensions.
		if (first)
		for (IDX_TYPE dim_out=0; dim_out<NODE_EMB_DIM; dim_out++) {
			# pragma HLS UNROLL
			DT1_TYPE bias = g_edge_nw_bias_0[dim_out];
			for (IDX_TYPE edge_idx_offset=0; edge_idx_offset<EDGE_UNI_PARA; edge_idx_offset++)
				# pragma HLS UNROLL
				linear_out_cache[edge_idx_offset][dim_out] = bias; }

		// Process EDGE_UNI_PARA edges each time.
		for (IDX_TYPE edge_idx_offset=0; edge_idx_offset<EDGE_UNI_PARA; edge_idx_offset++) {
			# pragma HLS UNROLL

			// Read the new node embeddings from FIFOs.
			EDGE_TYPE edge_cache;
			# pragma HLS AGGREGATE variable=edge_cache
			edge_in[edge_idx_offset].read(edge_cache);
			DT1_PARA emb_src = edge_cache.node_emb_src;
			DT1_PARA emb_dst = edge_cache.node_emb_dst;
			# pragma HLS AGGREGATE variable=emb_src
			# pragma HLS AGGREGATE variable=emb_dst

			// Bypass channel.
			if (!last_layer)
			edge_out[edge_idx_offset].write(edge_cache);

			// Get the weight and calculate product.
			for (IDX_TYPE dim_out=0; dim_out<NODE_EMB_DIM; dim_out++) {
				# pragma HLS UNROLL
				// Temporary variables of products.
				DT0_TYPE product_out = (DT0_TYPE)0.0;
				DT0_TYPE product_in = (DT0_TYPE)0.0;
				// Calculate products.
				for (IDX_TYPE dim_in_offset=0; dim_in_offset<DIM_PARA; dim_in_offset++) {
					# pragma HLS UNROLL
					// Get the input dimension.
					IDX_TYPE dim_in = dim_in_base + dim_in_offset;
					BOL_TYPE dim_in_overflow = (dim_in >= (IDX_TYPE)NODE_EMB_DIM);
					// Calculate products.
					DT1_TYPE product_out_cache = dim_in_overflow? (DT1_TYPE)0.0:
							(DT1_TYPE)(emb_src[dim_in_offset] * g_edge_nw_weight_0_0[dim_out][dim_in]);
					DT1_TYPE product_in_cache = dim_in_overflow? (DT1_TYPE)0.0:
							(DT1_TYPE)(emb_dst[dim_in_offset] * g_edge_nw_weight_0_1[dim_out][dim_in]);
					// Accumulate products.
					product_out += product_out_cache;
					product_in += product_in_cache; }
				// Get the old data and accumulate products.
				DT0_TYPE data_old = linear_out_cache[edge_idx_offset][dim_out];
				DT0_TYPE data_new = product_out + product_in + data_old;
				linear_out_cache[edge_idx_offset][dim_out] = data_new; }}

		// Send the linear outputs to the next layer.
		if (last)
		for (IDX_TYPE edge_idx_offset=0; edge_idx_offset<EDGE_UNI_PARA; edge_idx_offset++)
			# pragma HLS UNROLL
			linear_out[edge_idx_offset].write(linear_out_cache[edge_idx_offset]);
	}
}

static void edge_linear_1(
		BOL_TYPE last_layer,
		IDX_TYPE num_of_edges,
		hls::stream<DT1_PARA> linear_in[EDGE_UNI_PARA],
		hls::stream<DT0_VECTOR> linear_out[EDGE_UNI_PARA],
		hls::stream<EDGE_TYPE> edge_in[EDGE_UNI_PARA],
		hls::stream<EDGE_TYPE> edge_out[EDGE_UNI_PARA])
{
	# pragma HLS INLINE off
	# pragma HLS ARRAY_PARTITION variable=g_edge_nw_weight_1 dim=1 type=complete
	# pragma HLS ARRAY_PARTITION variable=g_edge_nw_weight_1 dim=2 type=cyclic factor=DIM_PARA
	# pragma HLS ARRAY_PARTITION variable=g_edge_nw_bias_1 dim=1 type=complete

	// Temporary array of linear outputs.
	DT0_VECTOR linear_out_cache[EDGE_UNI_PARA];
	# pragma HLS AGGREGATE variable=linear_out_cache
	# pragma HLS ARRAY_PARTITION variable=linear_out_cache dim=1 type=complete

	// Process EDGE_UNI_PARA edges each time.
	// Process DIM_PARA dimensions each time.
	IDX_TYPE iter_uni_num = ceildiv(num_of_edges, (IDX_TYPE)EDGE_UNI_PARA);
	IDX_TYPE iter_dim_num = ceildiv((IDX_TYPE)NODE_EMB_DIM, (IDX_TYPE)DIM_PARA);
	for (	IDX_TYPE iter_uni_cnt=(IDX_TYPE)0,
			edge_idx_base=(IDX_TYPE)0;
			iter_uni_cnt<iter_uni_num;
			iter_uni_cnt+=(IDX_TYPE)1,
			edge_idx_base+=(IDX_TYPE)EDGE_UNI_PARA)
	# pragma HLS LOOP_TRIPCOUNT min=EDGE_UNI_COUNT max=EDGE_UNI_COUNT avg=EDGE_UNI_COUNT
	for (	IDX_TYPE iter_dim_cnt=(IDX_TYPE)0,
			dim_in_base=(IDX_TYPE)0;
			iter_dim_cnt<iter_dim_num;
			iter_dim_cnt+=(IDX_TYPE)1,
			dim_in_base+=(IDX_TYPE)DIM_PARA)
	{
		# pragma HLS LOOP_TRIPCOUNT min=DIM_COUNT max=DIM_COUNT avg=DIM_COUNT
		# pragma HLS PIPELINE II=1 style=frp

		// Whether it is the first or the last iteration.
		BOL_TYPE first = (dim_in_base == (IDX_TYPE)0);
		BOL_TYPE last = (dim_in_base + (IDX_TYPE)DIM_PARA >= (IDX_TYPE)NODE_EMB_DIM);

		// Get biases of all the output dimensions.
		if (first)
		for (IDX_TYPE dim_out=0; dim_out<NODE_EMB_DIM; dim_out++) {
			# pragma HLS UNROLL
			DT1_TYPE bias = g_edge_nw_bias_1[dim_out];
			for (IDX_TYPE edge_idx_offset=0; edge_idx_offset<EDGE_UNI_PARA; edge_idx_offset++)
				# pragma HLS UNROLL
				linear_out_cache[edge_idx_offset][dim_out] = bias; }

		// Process EDGE_UNI_PARA edges each time.
		for (IDX_TYPE edge_idx_offset=0; edge_idx_offset<EDGE_UNI_PARA; edge_idx_offset++) {
			# pragma HLS UNROLL

			// Read data from FIFOs.
			DT1_PARA input;
			# pragma HLS AGGREGATE variable=input
			linear_in[edge_idx_offset].read(input);

			// Get the weight and calculate product.
			for (IDX_TYPE dim_out=0; dim_out<NODE_EMB_DIM; dim_out++) {
				# pragma HLS UNROLL
				// Calculate the sum of products.
				DT0_TYPE sum_of_product = (DT0_TYPE)0.0;
				for (IDX_TYPE dim_in_offset=0; dim_in_offset<DIM_PARA; dim_in_offset++) {
					# pragma HLS UNROLL
					IDX_TYPE dim_in = dim_in_base + dim_in_offset;
					BOL_TYPE dim_in_overflow = (dim_in >= (IDX_TYPE)NODE_EMB_DIM);
					DT1_TYPE product_cache = dim_in_overflow? (DT1_TYPE)0.0:
							(DT1_TYPE)(input[dim_in_offset] * g_edge_nw_weight_1[dim_out][dim_in]);
					sum_of_product += product_cache; }
				// Get the old data and accumulate products.
				DT0_TYPE data_old = linear_out_cache[edge_idx_offset][dim_out];
				DT0_TYPE data_new = sum_of_product + data_old;
				linear_out_cache[edge_idx_offset][dim_out] = data_new; }}

		// Bypass channel.
		if (!last_layer)
		for (IDX_TYPE edge_idx_offset=0; edge_idx_offset<EDGE_UNI_PARA; edge_idx_offset++) {
			# pragma HLS UNROLL
			EDGE_TYPE edge_cache;
			# pragma HLS AGGREGATE variable=edge_cache
			edge_in[edge_idx_offset].read(edge_cache);
			edge_out[edge_idx_offset].write(edge_cache); }

		// Send the linear outputs to the next layer.
		if (last)
		for (IDX_TYPE edge_idx_offset=0; edge_idx_offset<EDGE_UNI_PARA; edge_idx_offset++)
			# pragma HLS UNROLL
			linear_out[edge_idx_offset].write(linear_out_cache[edge_idx_offset]);

	}
}

static void edge_linear_2(
		BOL_TYPE last_layer,
		IDX_TYPE num_of_edges,
		hls::stream<DT1_PARA> linear_in[EDGE_UNI_PARA],
		hls::stream<DT0_VECTOR> linear_out[EDGE_UNI_PARA],
		hls::stream<EDGE_TYPE> edge_in[EDGE_UNI_PARA],
		hls::stream<EDGE_TYPE> edge_out[EDGE_UNI_PARA])
{
	# pragma HLS INLINE off
	# pragma HLS ARRAY_PARTITION variable=g_edge_nw_weight_2 dim=1 type=complete
	# pragma HLS ARRAY_PARTITION variable=g_edge_nw_weight_2 dim=2 type=cyclic factor=DIM_PARA
	# pragma HLS ARRAY_PARTITION variable=g_edge_nw_bias_2 dim=1 type=complete

	// Temporary array of linear outputs.
	DT0_VECTOR linear_out_cache[EDGE_UNI_PARA];
	# pragma HLS AGGREGATE variable=linear_out_cache
	# pragma HLS ARRAY_PARTITION variable=linear_out_cache dim=1 type=complete

	// Process EDGE_UNI_PARA edges each time.
	// Process DIM_PARA dimensions each time.
	IDX_TYPE iter_uni_num = ceildiv(num_of_edges, (IDX_TYPE)EDGE_UNI_PARA);
	IDX_TYPE iter_dim_num = ceildiv((IDX_TYPE)NODE_EMB_DIM, (IDX_TYPE)DIM_PARA);
	for (	IDX_TYPE iter_uni_cnt=(IDX_TYPE)0,
			edge_idx_base=(IDX_TYPE)0;
			iter_uni_cnt<iter_uni_num;
			iter_uni_cnt+=(IDX_TYPE)1,
			edge_idx_base+=(IDX_TYPE)EDGE_UNI_PARA)
	# pragma HLS LOOP_TRIPCOUNT min=EDGE_UNI_COUNT max=EDGE_UNI_COUNT avg=EDGE_UNI_COUNT
	for (	IDX_TYPE iter_dim_cnt=(IDX_TYPE)0,
			dim_in_base=(IDX_TYPE)0;
			iter_dim_cnt<iter_dim_num;
			iter_dim_cnt+=(IDX_TYPE)1,
			dim_in_base+=(IDX_TYPE)DIM_PARA)
	{
		# pragma HLS LOOP_TRIPCOUNT min=DIM_COUNT max=DIM_COUNT avg=DIM_COUNT
		# pragma HLS PIPELINE II=1 style=frp

		// Whether it is the first or the last iteration.
		BOL_TYPE first = (dim_in_base == (IDX_TYPE)0);
		BOL_TYPE last = (dim_in_base + (IDX_TYPE)DIM_PARA >= (IDX_TYPE)NODE_EMB_DIM);

		// Get biases of all the output dimensions.
		if (first)
		for (IDX_TYPE dim_out=0; dim_out<NODE_EMB_DIM; dim_out++) {
			# pragma HLS UNROLL
			DT1_TYPE bias = g_edge_nw_bias_2[dim_out];
			for (IDX_TYPE edge_idx_offset=0; edge_idx_offset<EDGE_UNI_PARA; edge_idx_offset++)
				# pragma HLS UNROLL
				linear_out_cache[edge_idx_offset][dim_out] = bias; }

		// Process EDGE_UNI_PARA edges each time.
		for (IDX_TYPE edge_idx_offset=0; edge_idx_offset<EDGE_UNI_PARA; edge_idx_offset++) {
			# pragma HLS UNROLL

			// Read data from FIFOs.
			DT1_PARA input;
			# pragma HLS AGGREGATE variable=input
			linear_in[edge_idx_offset].read(input);

			// Get the weight and calculate product.
			for (IDX_TYPE dim_out=0; dim_out<NODE_EMB_DIM; dim_out++) {
				# pragma HLS UNROLL
				// Calculate the sum of products.
				DT0_TYPE sum_of_product = (DT0_TYPE)0.0;
				for (IDX_TYPE dim_in_offset=0; dim_in_offset<DIM_PARA; dim_in_offset++) {
					# pragma HLS UNROLL
					IDX_TYPE dim_in = dim_in_base + dim_in_offset;
					BOL_TYPE dim_in_overflow = (dim_in >= (IDX_TYPE)NODE_EMB_DIM);
					DT1_TYPE product_cache = dim_in_overflow? (DT1_TYPE)0.0:
							(DT1_TYPE)(input[dim_in_offset] * g_edge_nw_weight_2[dim_out][dim_in]);
					sum_of_product += product_cache; }
				// Get the old data and accumulate products.
				DT0_TYPE data_old = linear_out_cache[edge_idx_offset][dim_out];
				DT0_TYPE data_new = sum_of_product + data_old;
				linear_out_cache[edge_idx_offset][dim_out] = data_new; }}

		// Bypass channel.
		if (!last_layer)
		for (IDX_TYPE edge_idx_offset=0; edge_idx_offset<EDGE_UNI_PARA; edge_idx_offset++) {
			# pragma HLS UNROLL
			EDGE_TYPE edge_cache;
			# pragma HLS AGGREGATE variable=edge_cache
			edge_in[edge_idx_offset].read(edge_cache);
			edge_out[edge_idx_offset].write(edge_cache); }

		// Send the linear outputs to the next layer.
		if (last)
		for (IDX_TYPE edge_idx_offset=0; edge_idx_offset<EDGE_UNI_PARA; edge_idx_offset++)
			# pragma HLS UNROLL
			linear_out[edge_idx_offset].write(linear_out_cache[edge_idx_offset]);
	}
}

static void edge_linear_3(
		BOL_TYPE last_layer,
		IDX_TYPE num_of_edges,
		hls::stream<DT1_PARA> linear_in[EDGE_UNI_PARA],
		hls::stream<DT0_TYPE> linear_out[EDGE_UNI_PARA],
		hls::stream<EDGE_TYPE> edge_in[EDGE_UNI_PARA],
		hls::stream<EDGE_TYPE> edge_out[EDGE_UNI_PARA])
{
	# pragma HLS INLINE off
	# pragma HLS ARRAY_PARTITION variable=g_edge_nw_weight_3 dim=1 type=cyclic factor=DIM_PARA

	// Temporary array of linear outputs.
	DT0_TYPE linear_out_cache[EDGE_UNI_PARA];
	# pragma HLS ARRAY_PARTITION variable=linear_out_cache dim=1 type=complete

	// Process EDGE_UNI_PARA edges each time.
	// Process DIM_PARA dimensions each time.
	IDX_TYPE iter_uni_num = ceildiv(num_of_edges, (IDX_TYPE)EDGE_UNI_PARA);
	IDX_TYPE iter_dim_num = ceildiv((IDX_TYPE)NODE_EMB_DIM, (IDX_TYPE)DIM_PARA);
	for (	IDX_TYPE iter_uni_cnt=(IDX_TYPE)0,
			edge_idx_base=(IDX_TYPE)0;
			iter_uni_cnt<iter_uni_num;
			iter_uni_cnt+=(IDX_TYPE)1,
			edge_idx_base+=(IDX_TYPE)EDGE_UNI_PARA)
	# pragma HLS LOOP_TRIPCOUNT min=EDGE_UNI_COUNT max=EDGE_UNI_COUNT avg=EDGE_UNI_COUNT
	for (	IDX_TYPE iter_dim_cnt=(IDX_TYPE)0,
			dim_in_base=(IDX_TYPE)0;
			iter_dim_cnt<iter_dim_num;
			iter_dim_cnt+=(IDX_TYPE)1,
			dim_in_base+=(IDX_TYPE)DIM_PARA)
	{
		# pragma HLS LOOP_TRIPCOUNT min=DIM_COUNT max=DIM_COUNT avg=DIM_COUNT
		# pragma HLS PIPELINE II=1 style=frp

		// Whether it is the first or the last iteration.
		BOL_TYPE first = (dim_in_base == (IDX_TYPE)0);
		BOL_TYPE last = (dim_in_base + (IDX_TYPE)DIM_PARA >= (IDX_TYPE)NODE_EMB_DIM);

		// Get biases of all the output dimensions.
		if (first) {
			DT1_TYPE bias = g_edge_nw_bias_3;
			for (IDX_TYPE edge_idx_offset=0; edge_idx_offset<EDGE_UNI_PARA; edge_idx_offset++)
				# pragma HLS UNROLL
				linear_out_cache[edge_idx_offset] = bias; }

		// Process EDGE_UNI_PARA edges each time.
		for (IDX_TYPE edge_idx_offset=0; edge_idx_offset<EDGE_UNI_PARA; edge_idx_offset++) {
			# pragma HLS UNROLL

			// Read data from FIFOs.
			DT1_PARA input;
			# pragma HLS AGGREGATE variable=input
			linear_in[edge_idx_offset].read(input);

			// Get the weight and calculate product.
			DT0_TYPE sum_of_product = (DT0_TYPE)0.0;
			for (IDX_TYPE dim_in_offset=0; dim_in_offset<DIM_PARA; dim_in_offset++) {
				# pragma HLS UNROLL
				IDX_TYPE dim_in = dim_in_base + dim_in_offset;
				BOL_TYPE dim_in_overflow = (dim_in >= (IDX_TYPE)NODE_EMB_DIM);
				DT1_TYPE product_cache = dim_in_overflow? (DT1_TYPE)0.0:
						(DT1_TYPE)(input[dim_in_offset] * g_edge_nw_weight_3[dim_in]);
				sum_of_product += product_cache; }
			// Get the old data and accumulate products.
			DT0_TYPE data_old = linear_out_cache[edge_idx_offset];
			DT0_TYPE data_new = sum_of_product + data_old;
			linear_out_cache[edge_idx_offset] = data_new; }

		// Bypass channel.
		if (!last_layer)
		for (IDX_TYPE edge_idx_offset=0; edge_idx_offset<EDGE_UNI_PARA; edge_idx_offset++) {
			# pragma HLS UNROLL
			EDGE_TYPE edge_cache;
			# pragma HLS AGGREGATE variable=edge_cache
			edge_in[edge_idx_offset].read(edge_cache);
			edge_out[edge_idx_offset].write(edge_cache); }

		// Send the linear outputs to the next layer.
		if (last)
		for (IDX_TYPE edge_idx_offset=0; edge_idx_offset<EDGE_UNI_PARA; edge_idx_offset++)
			# pragma HLS UNROLL
			linear_out[edge_idx_offset].write(linear_out_cache[edge_idx_offset]);
	}
}

static void edge_activation_0(
		BOL_TYPE last_layer,
		IDX_TYPE num_of_edges,
		hls::stream<DT0_VECTOR> activation_in[EDGE_UNI_PARA],
		hls::stream<DT1_PARA> activation_out[EDGE_UNI_PARA],
		hls::stream<EDGE_TYPE> edge_in[EDGE_UNI_PARA],
		hls::stream<EDGE_TYPE> edge_out[EDGE_UNI_PARA])
{
	# pragma HLS INLINE off

	// Temporary array of activation inputs.
	DT0_TYPE act_in[EDGE_UNI_PARA][NODE_EMB_DIM];
	# pragma HLS ARRAY_PARTITION variable=act_in dim=1 type=complete
	# pragma HLS ARRAY_PARTITION variable=act_in dim=2 type=complete

	// Process EDGE_UNI_PARA edges each time.
	// Process DIM_PARA dimensions each time.
	IDX_TYPE iter_uni_num = ceildiv(num_of_edges, (IDX_TYPE)EDGE_UNI_PARA);
	IDX_TYPE iter_dim_num = ceildiv((IDX_TYPE)NODE_EMB_DIM, (IDX_TYPE)DIM_PARA);
	for (	IDX_TYPE iter_uni_cnt=(IDX_TYPE)0,
			edge_idx_base=(IDX_TYPE)0;
			iter_uni_cnt<iter_uni_num;
			iter_uni_cnt+=(IDX_TYPE)1,
			edge_idx_base+=(IDX_TYPE)EDGE_UNI_PARA)
	# pragma HLS LOOP_TRIPCOUNT min=EDGE_UNI_COUNT max=EDGE_UNI_COUNT avg=EDGE_UNI_COUNT
	for (	IDX_TYPE iter_dim_cnt=(IDX_TYPE)0,
			dim_out_base=(IDX_TYPE)0;
			iter_dim_cnt<iter_dim_num;
			iter_dim_cnt+=(IDX_TYPE)1,
			dim_out_base+=(IDX_TYPE)DIM_PARA)
	{
		# pragma HLS LOOP_TRIPCOUNT min=DIM_COUNT max=DIM_COUNT avg=DIM_COUNT
		# pragma HLS PIPELINE II=1 style=frp

		// Whether it is the first iteration.
		BOL_TYPE first = (dim_out_base == (IDX_TYPE)0);

		// Get the activation inputs.
		if (first)
		for (IDX_TYPE edge_idx_offset=0; edge_idx_offset<EDGE_UNI_PARA; edge_idx_offset++) {
			# pragma HLS UNROLL
			DT0_VECTOR act_in_cache;
			# pragma HLS AGGREGATE variable=act_in_cache
			activation_in[edge_idx_offset].read(act_in_cache);
			for (IDX_TYPE dim_in=0; dim_in<NODE_EMB_DIM; dim_in++)
				# pragma HLS UNROLL
				act_in[edge_idx_offset][dim_in] = act_in_cache[dim_in]; }

		// Process EDGE_UNI_PARA edges each time.
		for (IDX_TYPE edge_idx_offset=0; edge_idx_offset<EDGE_UNI_PARA; edge_idx_offset++) {
			# pragma HLS UNROLL
			// Temporary array of activation outputs.
			DT1_PARA act_out;
			# pragma HLS AGGREGATE variable=act_out
			// Calculate activation outputs.
			for (IDX_TYPE dim_out_offset=0; dim_out_offset<DIM_PARA; dim_out_offset++) {
				# pragma HLS UNROLL
				// Get the output dimension.
				IDX_TYPE dim_out = dim_out_base + dim_out_offset;
				BOL_TYPE dim_out_overflow = (dim_out >= (IDX_TYPE)NODE_EMB_DIM);
				// Pass through activation functions.
				# ifdef ACT_TANH
				act_out[dim_out_offset] = dim_out_overflow? (DT1_TYPE)0.0:
						act_tanh_edge_nw((IDX_TYPE)0, edge_idx_offset, dim_out_offset, act_in[edge_idx_offset][dim_out]); }
				# endif
				# ifdef ACT_SIGM
				act_out[dim_out_offset] = dim_out_overflow? (DT1_TYPE)0.0:
						act_sigm_edge_nw((IDX_TYPE)0, edge_idx_offset, dim_out_offset, act_in[edge_idx_offset][dim_out]); }
				# endif
			// Send the activation outputs to the next layer.
			activation_out[edge_idx_offset].write(act_out); }

		// Bypass channel.
		if (!last_layer)
		for (IDX_TYPE edge_idx_offset=0; edge_idx_offset<EDGE_UNI_PARA; edge_idx_offset++) {
			# pragma HLS UNROLL
			EDGE_TYPE edge_cache;
			# pragma HLS AGGREGATE variable=edge_cache
			edge_in[edge_idx_offset].read(edge_cache);
			edge_out[edge_idx_offset].write(edge_cache); }
	}
}

static void edge_activation_1(
		BOL_TYPE last_layer,
		IDX_TYPE num_of_edges,
		hls::stream<DT0_VECTOR> activation_in[EDGE_UNI_PARA],
		hls::stream<DT1_PARA> activation_out[EDGE_UNI_PARA],
		hls::stream<EDGE_TYPE> edge_in[EDGE_UNI_PARA],
		hls::stream<EDGE_TYPE> edge_out[EDGE_UNI_PARA])
{
	# pragma HLS INLINE off

	// Temporary array of activation inputs.
	DT0_TYPE act_in[EDGE_UNI_PARA][NODE_EMB_DIM];
	# pragma HLS ARRAY_PARTITION variable=act_in dim=1 type=complete
	# pragma HLS ARRAY_PARTITION variable=act_in dim=2 type=complete

	// Process EDGE_UNI_PARA edges each time.
	// Process DIM_PARA dimensions each time.
	IDX_TYPE iter_uni_num = ceildiv(num_of_edges, (IDX_TYPE)EDGE_UNI_PARA);
	IDX_TYPE iter_dim_num = ceildiv((IDX_TYPE)NODE_EMB_DIM, (IDX_TYPE)DIM_PARA);
	for (	IDX_TYPE iter_uni_cnt=(IDX_TYPE)0,
			edge_idx_base=(IDX_TYPE)0;
			iter_uni_cnt<iter_uni_num;
			iter_uni_cnt+=(IDX_TYPE)1,
			edge_idx_base+=(IDX_TYPE)EDGE_UNI_PARA)
	# pragma HLS LOOP_TRIPCOUNT min=EDGE_UNI_COUNT max=EDGE_UNI_COUNT avg=EDGE_UNI_COUNT
	for (	IDX_TYPE iter_dim_cnt=(IDX_TYPE)0,
			dim_out_base=(IDX_TYPE)0;
			iter_dim_cnt<iter_dim_num;
			iter_dim_cnt+=(IDX_TYPE)1,
			dim_out_base+=(IDX_TYPE)DIM_PARA)
	{
		# pragma HLS LOOP_TRIPCOUNT min=DIM_COUNT max=DIM_COUNT avg=DIM_COUNT
		# pragma HLS PIPELINE II=1 style=frp

		// Whether it is the first iteration.
		BOL_TYPE first = (dim_out_base == (IDX_TYPE)0);

		// Get the activation inputs.
		if (first)
		for (IDX_TYPE edge_idx_offset=0; edge_idx_offset<EDGE_UNI_PARA; edge_idx_offset++) {
			# pragma HLS UNROLL
			DT0_VECTOR act_in_cache;
			# pragma HLS AGGREGATE variable=act_in_cache
			activation_in[edge_idx_offset].read(act_in_cache);
			for (IDX_TYPE dim_in=0; dim_in<NODE_EMB_DIM; dim_in++)
				# pragma HLS UNROLL
				act_in[edge_idx_offset][dim_in] = act_in_cache[dim_in]; }

		// Process EDGE_UNI_PARA edges each time.
		for (IDX_TYPE edge_idx_offset=0; edge_idx_offset<EDGE_UNI_PARA; edge_idx_offset++) {
			# pragma HLS UNROLL
			// Temporary array of activation outputs.
			DT1_PARA act_out;
			# pragma HLS AGGREGATE variable=act_out
			// Calculate activation outputs.
			for (IDX_TYPE dim_out_offset=0; dim_out_offset<DIM_PARA; dim_out_offset++) {
				# pragma HLS UNROLL
				// Get the output dimension.
				IDX_TYPE dim_out = dim_out_base + dim_out_offset;
				BOL_TYPE dim_out_overflow = (dim_out >= (IDX_TYPE)NODE_EMB_DIM);
				// Pass through activation functions.
				# ifdef ACT_TANH
				act_out[dim_out_offset] = dim_out_overflow? (DT1_TYPE)0.0:
						act_tanh_edge_nw((IDX_TYPE)1, edge_idx_offset, dim_out_offset, act_in[edge_idx_offset][dim_out]); }
				# endif
				# ifdef ACT_SIGM
				act_out[dim_out_offset] = dim_out_overflow? (DT1_TYPE)0.0:
						act_sigm_edge_nw((IDX_TYPE)1, edge_idx_offset, dim_out_offset, act_in[edge_idx_offset][dim_out]); }
				# endif
			// Send the activation outputs to the next layer.
			activation_out[edge_idx_offset].write(act_out); }

		// Bypass channel.
		if (!last_layer)
		for (IDX_TYPE edge_idx_offset=0; edge_idx_offset<EDGE_UNI_PARA; edge_idx_offset++) {
			# pragma HLS UNROLL
			EDGE_TYPE edge_cache;
			# pragma HLS AGGREGATE variable=edge_cache
			edge_in[edge_idx_offset].read(edge_cache);
			edge_out[edge_idx_offset].write(edge_cache); }
	}
}

static void edge_activation_2(
		BOL_TYPE last_layer,
		IDX_TYPE num_of_edges,
		hls::stream<DT0_VECTOR> activation_in[EDGE_UNI_PARA],
		hls::stream<DT1_PARA> activation_out[EDGE_UNI_PARA],
		hls::stream<EDGE_TYPE> edge_in[EDGE_UNI_PARA],
		hls::stream<EDGE_TYPE> edge_out[EDGE_UNI_PARA])
{
	# pragma HLS INLINE off

	// Temporary array of activation inputs.
	DT0_TYPE act_in[EDGE_UNI_PARA][NODE_EMB_DIM];
	# pragma HLS ARRAY_PARTITION variable=act_in dim=1 type=complete
	# pragma HLS ARRAY_PARTITION variable=act_in dim=2 type=complete

	// Process EDGE_UNI_PARA edges each time.
	// Process DIM_PARA dimensions each time.
	IDX_TYPE iter_uni_num = ceildiv(num_of_edges, (IDX_TYPE)EDGE_UNI_PARA);
	IDX_TYPE iter_dim_num = ceildiv((IDX_TYPE)NODE_EMB_DIM, (IDX_TYPE)DIM_PARA);
	for (	IDX_TYPE iter_uni_cnt=(IDX_TYPE)0,
			edge_idx_base=(IDX_TYPE)0;
			iter_uni_cnt<iter_uni_num;
			iter_uni_cnt+=(IDX_TYPE)1,
			edge_idx_base+=(IDX_TYPE)EDGE_UNI_PARA)
	# pragma HLS LOOP_TRIPCOUNT min=EDGE_UNI_COUNT max=EDGE_UNI_COUNT avg=EDGE_UNI_COUNT
	for (	IDX_TYPE iter_dim_cnt=(IDX_TYPE)0,
			dim_out_base=(IDX_TYPE)0;
			iter_dim_cnt<iter_dim_num;
			iter_dim_cnt+=(IDX_TYPE)1,
			dim_out_base+=(IDX_TYPE)DIM_PARA)
	{
		# pragma HLS LOOP_TRIPCOUNT min=DIM_COUNT max=DIM_COUNT avg=DIM_COUNT
		# pragma HLS PIPELINE II=1 style=frp

		// Whether it is the first iteration.
		BOL_TYPE first = (dim_out_base == (IDX_TYPE)0);

		// Get the activation inputs.
		if (first)
		for (IDX_TYPE edge_idx_offset=0; edge_idx_offset<EDGE_UNI_PARA; edge_idx_offset++) {
			# pragma HLS UNROLL
			DT0_VECTOR act_in_cache;
			# pragma HLS AGGREGATE variable=act_in_cache
			activation_in[edge_idx_offset].read(act_in_cache);
			for (IDX_TYPE dim_in=0; dim_in<NODE_EMB_DIM; dim_in++)
				# pragma HLS UNROLL
				act_in[edge_idx_offset][dim_in] = act_in_cache[dim_in]; }

		// Process EDGE_UNI_PARA edges each time.
		for (IDX_TYPE edge_idx_offset=0; edge_idx_offset<EDGE_UNI_PARA; edge_idx_offset++) {
			# pragma HLS UNROLL
			// Temporary array of activation outputs.
			DT1_PARA act_out;
			# pragma HLS AGGREGATE variable=act_out
			// Calculate activation outputs.
			for (IDX_TYPE dim_out_offset=0; dim_out_offset<DIM_PARA; dim_out_offset++) {
				# pragma HLS UNROLL
				// Get the output dimension.
				IDX_TYPE dim_out = dim_out_base + dim_out_offset;
				BOL_TYPE dim_out_overflow = (dim_out >= (IDX_TYPE)NODE_EMB_DIM);
				// Pass through activation functions.
				# ifdef ACT_TANH
				act_out[dim_out_offset] = dim_out_overflow? (DT1_TYPE)0.0:
						act_tanh_edge_nw((IDX_TYPE)2, edge_idx_offset, dim_out_offset, act_in[edge_idx_offset][dim_out]); }
				# endif
				# ifdef ACT_SIGM
				act_out[dim_out_offset] = dim_out_overflow? (DT1_TYPE)0.0:
						act_sigm_edge_nw((IDX_TYPE)2, edge_idx_offset, dim_out_offset, act_in[edge_idx_offset][dim_out]); }
				# endif
			// Send the activation outputs to the next layer.
			activation_out[edge_idx_offset].write(act_out); }

		// Bypass channel.
		if (!last_layer)
		for (IDX_TYPE edge_idx_offset=0; edge_idx_offset<EDGE_UNI_PARA; edge_idx_offset++) {
			# pragma HLS UNROLL
			EDGE_TYPE edge_cache;
			# pragma HLS AGGREGATE variable=edge_cache
			edge_in[edge_idx_offset].read(edge_cache);
			edge_out[edge_idx_offset].write(edge_cache); }
	}
}

static void edge_activation_3(
		BOL_TYPE last_layer,
		IDX_TYPE num_of_edges,
		hls::stream<DT0_TYPE> activation_in[EDGE_UNI_PARA],
		hls::stream<DT1_TYPE> activation_out[EDGE_UNI_PARA],
		hls::stream<EDGE_TYPE> edge_in[EDGE_UNI_PARA],
		hls::stream<EDGE_TYPE> edge_out[EDGE_UNI_PARA])
{
	# pragma HLS INLINE off

	// Process EDGE_UNI_PARA edges each time.
	// Process DIM_PARA dimensions each time.
	IDX_TYPE iter_uni_num = ceildiv(num_of_edges, (IDX_TYPE)EDGE_UNI_PARA);
	IDX_TYPE iter_dim_num = ceildiv((IDX_TYPE)NODE_EMB_DIM, (IDX_TYPE)DIM_PARA);
	for (	IDX_TYPE iter_uni_cnt=(IDX_TYPE)0,
			edge_idx_base=(IDX_TYPE)0;
			iter_uni_cnt<iter_uni_num;
			iter_uni_cnt+=(IDX_TYPE)1,
			edge_idx_base+=(IDX_TYPE)EDGE_UNI_PARA)
	# pragma HLS LOOP_TRIPCOUNT min=EDGE_UNI_COUNT max=EDGE_UNI_COUNT avg=EDGE_UNI_COUNT
	for (	IDX_TYPE iter_dim_cnt=(IDX_TYPE)0,
			dim_out_base=(IDX_TYPE)0;
			iter_dim_cnt<iter_dim_num;
			iter_dim_cnt+=(IDX_TYPE)1,
			dim_out_base+=(IDX_TYPE)DIM_PARA)
	{
		# pragma HLS LOOP_TRIPCOUNT min=DIM_COUNT max=DIM_COUNT avg=DIM_COUNT
		# pragma HLS PIPELINE II=1 style=frp

		// Whether it is the first iteration.
		BOL_TYPE first = (dim_out_base == (IDX_TYPE)0);

		// Pass through activation functions.
		if (first)
		for (IDX_TYPE edge_idx_offset=0; edge_idx_offset<EDGE_UNI_PARA; edge_idx_offset++) {
			# pragma HLS UNROLL
			DT0_TYPE act_in;
			activation_in[edge_idx_offset].read(act_in);
			DT1_TYPE act_out = act_sigm_result(edge_idx_offset, act_in);
			activation_out[edge_idx_offset].write(act_out); }

		// Bypass channel.
		if (!last_layer)
		for (IDX_TYPE edge_idx_offset=0; edge_idx_offset<EDGE_UNI_PARA; edge_idx_offset++) {
			# pragma HLS UNROLL
			EDGE_TYPE edge_cache;
			# pragma HLS AGGREGATE variable=edge_cache
			edge_in[edge_idx_offset].read(edge_cache);
			edge_out[edge_idx_offset].write(edge_cache); }
	}
}

static void edge_nw(
		BOL_TYPE last_layer,
		IDX_TYPE num_of_edges,
		hls::stream<EDGE_TYPE> edge_stream[EDGE_UNI_PARA],
		hls::stream<DT1_TYPE> value_stream[EDGE_UNI_PARA],
		hls::stream<EDGE_TYPE> msg_stream[EDGE_UNI_PARA])
{
	# pragma HLS INLINE off
	# pragma HLS DATAFLOW

	// FIFO declarations.
	hls::stream<DT0_VECTOR> linear_out_0[EDGE_UNI_PARA];
	hls::stream<DT0_VECTOR> linear_out_1[EDGE_UNI_PARA];
	hls::stream<DT0_VECTOR> linear_out_2[EDGE_UNI_PARA];
	hls::stream<DT0_TYPE> linear_out_3[EDGE_UNI_PARA];
	# pragma HLS STREAM variable=linear_out_0 depth=FIFO_DEPTH
	# pragma HLS STREAM variable=linear_out_1 depth=FIFO_DEPTH
	# pragma HLS STREAM variable=linear_out_2 depth=FIFO_DEPTH
	# pragma HLS STREAM variable=linear_out_3 depth=FIFO_DEPTH
	hls::stream<EDGE_TYPE> linear_edge_out_0[EDGE_UNI_PARA];
	hls::stream<EDGE_TYPE> linear_edge_out_1[EDGE_UNI_PARA];
	hls::stream<EDGE_TYPE> linear_edge_out_2[EDGE_UNI_PARA];
	hls::stream<EDGE_TYPE> linear_edge_out_3[EDGE_UNI_PARA];
	# pragma HLS STREAM variable=linear_edge_out_0 depth=(FIFO_DEPTH*DIM_COUNT)
	# pragma HLS STREAM variable=linear_edge_out_1 depth=(FIFO_DEPTH*DIM_COUNT)
	# pragma HLS STREAM variable=linear_edge_out_2 depth=(FIFO_DEPTH*DIM_COUNT)
	# pragma HLS STREAM variable=linear_edge_out_3 depth=(FIFO_DEPTH*DIM_COUNT)
	hls::stream<DT1_PARA> activation_out_0[EDGE_UNI_PARA];
	hls::stream<DT1_PARA> activation_out_1[EDGE_UNI_PARA];
	hls::stream<DT1_PARA> activation_out_2[EDGE_UNI_PARA];
	# pragma HLS STREAM variable=activation_out_0 depth=FIFO_DEPTH
	# pragma HLS STREAM variable=activation_out_1 depth=FIFO_DEPTH
	# pragma HLS STREAM variable=activation_out_2 depth=FIFO_DEPTH
	hls::stream<EDGE_TYPE> activation_edge_out_0[EDGE_UNI_PARA];
	hls::stream<EDGE_TYPE> activation_edge_out_1[EDGE_UNI_PARA];
	hls::stream<EDGE_TYPE> activation_edge_out_2[EDGE_UNI_PARA];
	# pragma HLS STREAM variable=activation_edge_out_0 depth=FIFO_DEPTH
	# pragma HLS STREAM variable=activation_edge_out_1 depth=FIFO_DEPTH
	# pragma HLS STREAM variable=activation_edge_out_2 depth=FIFO_DEPTH

	// Invoke functions.
	edge_linear_0(
			last_layer,
			num_of_edges,
			linear_out_0,
			edge_stream,
			linear_edge_out_0);
	edge_activation_0(
			last_layer,
			num_of_edges,
			linear_out_0,
			activation_out_0,
			linear_edge_out_0,
			activation_edge_out_0);
	edge_linear_1(
			last_layer,
			num_of_edges,
			activation_out_0,
			linear_out_1,
			activation_edge_out_0,
			linear_edge_out_1);
	edge_activation_1(
			last_layer,
			num_of_edges,
			linear_out_1,
			activation_out_1,
			linear_edge_out_1,
			activation_edge_out_1);
	edge_linear_2(
			last_layer,
			num_of_edges,
			activation_out_1,
			linear_out_2,
			activation_edge_out_1,
			linear_edge_out_2);
	edge_activation_2(
			last_layer,
			num_of_edges,
			linear_out_2,
			activation_out_2,
			linear_edge_out_2,
			activation_edge_out_2);
	edge_linear_3(
			last_layer,
			num_of_edges,
			activation_out_2,
			linear_out_3,
			activation_edge_out_2,
			linear_edge_out_3);
	edge_activation_3(
			last_layer,
			num_of_edges,
			linear_out_3,
			value_stream,
			linear_edge_out_3,
			msg_stream);
}



// **************************** Scatter **************************** //

/*

	This part is to implement data scattering. The message scattering
	function can receive edge labels and the corresponding new node
	embeddings, and then scatter messages to the message arrays. The
	node scattering function can receive new node embeddings and send
	them to the node embedding array.

	Version: 0.1
		Author: Hanqing Zhang
		Email : hanqing.zhang@zju.edu.cn
		Date  : July-8-2023
		Commit: Add function "msg_scatter".

	Version: 0.2
		Author: Hanqing Zhang
		Email : hanqing.zhang@zju.edu.cn
		Date  : July-19-2023
		Commit: Combine msg_out and msg_in together to reduce the area.

	Version: 0.3
		Author: Hanqing Zhang
		Email : hanqing.zhang@zju.edu.cn
		Date  : August-1-2023
		Commit: Add function "node_scatter".

	Version: 0.4
		Author: Hanqing Zhang
		Email : hanqing.zhang@zju.edu.cn
		Date  : September-1-2023
		Commit: Modify for FIFO balancing.

	Version: 0.5
		Author: Hanqing Zhang
		Email : hanqing.zhang@zju.edu.cn
		Date  : October-20-2023
		Commit: Fixe the bug of last layer.

	Function:
		void node_scatter(
			IDX_TYPE num_of_nodes,
			hls::stream<DT1_PARA> node_fifo[NODE_UNI_PARA],
			DT1_TYPE node_emb_new[NODE_NUM_MAX][NODE_EMB_DIM]);
		Input arguments:
			num_of_nodes: the number of nodes.
			node_fifo: data from node-processing network.
		Output arguments:
			node_emb_new: new node embedding array.

	Function:
		void msg_scatter(
			IDX_TYPE num_of_edges,
			hls::stream<DT1_TYPE> value_stream[EDGE_UNI_PARA],
			hls::stream<EDGE_TYPE> msg_stream[EDGE_UNI_PARA],
			DT0_TYPE msg_src_new[2][EDGE_UNI_PARA][NODE_NUM_MAX][NODE_EMB_DIM],
			DT0_TYPE msg_dst_new[2][EDGE_UNI_PARA][NODE_NUM_MAX][NODE_EMB_DIM],
			IDX_TYPE msg_src_new_flag[EDGE_UNI_PARA][NODE_NUM_MAX],
			IDX_TYPE msg_dst_new_flag[EDGE_UNI_PARA][NODE_NUM_MAX]);
		Input arguments:
			num_of_edges: the number of edges.
			value_stream: attention from edge-processing network.
			msg_stream: new node embedding from edge-processing network.
		Output arguments:
			msg_src_new: new source message array.
			msg_dst_new: new destination message array.
			msg_src_new_flag: flag of where the new source message is.
			msg_dst_new_flag: flag of where the new destination message is.

*/

static void node_scatter(
		IDX_TYPE num_of_nodes,
		hls::stream<DT1_PARA> node_fifo[NODE_UNI_PARA],
		DT1_TYPE node_emb_new[NODE_NUM_MAX][NODE_EMB_DIM])
{
	# pragma HLS INLINE off

	// Process NODE_UNI_PARA nodes each time.
	// Process DIM_PARA dimension each time.
	IDX_TYPE iter_uni_num = ceildiv(num_of_nodes, (IDX_TYPE)NODE_UNI_PARA);
	IDX_TYPE iter_dim_num = ceildiv((IDX_TYPE)NODE_EMB_DIM, (IDX_TYPE)DIM_PARA);
	for (	IDX_TYPE iter_uni_cnt=(IDX_TYPE)0,
			node_idx_base=(IDX_TYPE)0;
			iter_uni_cnt<iter_uni_num;
			iter_uni_cnt+=(IDX_TYPE)1,
			node_idx_base+=(IDX_TYPE)NODE_UNI_PARA)
	# pragma HLS LOOP_TRIPCOUNT min=NODE_UNI_COUNT max=NODE_UNI_COUNT avg=NODE_UNI_COUNT
	for (	IDX_TYPE iter_dim_cnt=(IDX_TYPE)0,
			dim_base=(IDX_TYPE)0;
			iter_dim_cnt<iter_dim_num;
			iter_dim_cnt+=(IDX_TYPE)1,
			dim_base+=(IDX_TYPE)DIM_PARA)
	{
		# pragma HLS LOOP_TRIPCOUNT min=DIM_COUNT max=DIM_COUNT avg=DIM_COUNT
		# pragma HLS PIPELINE II=1 style=frp

		// Process NODE_UNI_PARA nodes each time.
		for (IDX_TYPE node_idx_offset=0; node_idx_offset<NODE_UNI_PARA; node_idx_offset++) {
			# pragma HLS UNROLL
			// Get the node index.
			IDX_TYPE node_idx = node_idx_base + node_idx_offset;
			BOL_TYPE node_idx_overflow = (node_idx >= num_of_nodes);
			// Receive the new node embedding.
			DT1_PARA new_emb;
			# pragma HLS AGGREGATE variable=new_emb
			node_fifo[node_idx_offset].read(new_emb);
			// Update the node embedding.
			if (!node_idx_overflow)
			for (IDX_TYPE dim_offset=0; dim_offset<DIM_PARA; dim_offset++) {
				# pragma HLS UNROLL
				IDX_TYPE dim = dim_base + dim_offset;
				BOL_TYPE dim_overflow = (dim >= (IDX_TYPE)NODE_EMB_DIM);
				if (!dim_overflow)
				node_emb_new[node_idx][dim] = new_emb[dim_offset]; }}
	}
}

static void msg_scatter(
		IDX_TYPE num_of_edges,
		hls::stream<DT1_TYPE> value_stream[EDGE_UNI_PARA],
		hls::stream<EDGE_TYPE> msg_stream[EDGE_UNI_PARA],
		DT0_TYPE msg_src_new[2][EDGE_UNI_PARA][NODE_NUM_MAX][NODE_EMB_DIM],
		DT0_TYPE msg_dst_new[2][EDGE_UNI_PARA][NODE_NUM_MAX][NODE_EMB_DIM],
		IDX_TYPE msg_src_new_flag[EDGE_UNI_PARA][NODE_NUM_MAX],
		IDX_TYPE msg_dst_new_flag[EDGE_UNI_PARA][NODE_NUM_MAX])
{
	# pragma HLS INLINE off

	// Declarations of arrays.
	IDX_TYPE edge_idx[EDGE_UNI_PARA];
	BOL_TYPE edge_idx_overflow[EDGE_UNI_PARA];
	DT1_TYPE read_value_stream[EDGE_UNI_PARA];
	IDX_TYPE src_idx[EDGE_UNI_PARA];
	IDX_TYPE dst_idx[EDGE_UNI_PARA];
	IDX_TYPE src_flag[EDGE_UNI_PARA];
	IDX_TYPE dst_flag[EDGE_UNI_PARA];
	# pragma HLS ARRAY_PARTITION dim=1 type=complete variable=edge_idx
	# pragma HLS ARRAY_PARTITION dim=1 type=complete variable=edge_idx_overflow
	# pragma HLS ARRAY_PARTITION dim=1 type=complete variable=read_value_stream
	# pragma HLS ARRAY_PARTITION dim=1 type=complete variable=src_idx
	# pragma HLS ARRAY_PARTITION dim=1 type=complete variable=dst_idx
	# pragma HLS ARRAY_PARTITION dim=1 type=complete variable=src_flag
	# pragma HLS ARRAY_PARTITION dim=1 type=complete variable=dst_flag

	// Process EDGE_UNI_PARA edges each time.
	// Process DIM_PARA dimension each time.
	IDX_TYPE iter_uni_num = ceildiv(num_of_edges, (IDX_TYPE)EDGE_UNI_PARA);
	IDX_TYPE iter_dim_num = ceildiv((IDX_TYPE)NODE_EMB_DIM, (IDX_TYPE)DIM_PARA);
	for (	IDX_TYPE iter_uni_cnt=(IDX_TYPE)0,
			edge_idx_base=(IDX_TYPE)0;
			iter_uni_cnt<iter_uni_num;
			iter_uni_cnt+=(IDX_TYPE)1,
			edge_idx_base+=(IDX_TYPE)EDGE_UNI_PARA)
	# pragma HLS LOOP_TRIPCOUNT min=EDGE_UNI_COUNT max=EDGE_UNI_COUNT avg=EDGE_UNI_COUNT
	for (	IDX_TYPE iter_dim_cnt=(IDX_TYPE)0,
			dim_base=(IDX_TYPE)0;
			iter_dim_cnt<iter_dim_num;
			iter_dim_cnt+=(IDX_TYPE)1,
			dim_base+=(IDX_TYPE)DIM_PARA)
	{
		# pragma HLS LOOP_TRIPCOUNT min=DIM_COUNT max=DIM_COUNT avg=DIM_COUNT
		# pragma HLS PIPELINE II=1 style=frp
		# pragma HLS DEPENDENCE variable=msg_src_new_flag inter true distance=DIM_COUNT
		# pragma HLS DEPENDENCE variable=msg_dst_new_flag inter true distance=DIM_COUNT

		// Whether it is the first iteration.
		BOL_TYPE first = (dim_base == (IDX_TYPE)0);

		// If it is the first iteration.
		if (first) {
			// Get edge indexes.
			for (IDX_TYPE edge_idx_offset=0; edge_idx_offset<EDGE_UNI_PARA; edge_idx_offset++) {
				# pragma HLS UNROLL
				edge_idx[edge_idx_offset] = edge_idx_base + edge_idx_offset;
				edge_idx_overflow[edge_idx_offset] = (edge_idx[edge_idx_offset] >= num_of_edges); }
			// Fetch results from FIFOs.
			for (IDX_TYPE edge_idx_offset=0; edge_idx_offset<EDGE_UNI_PARA; edge_idx_offset++) {
				# pragma HLS UNROLL
				value_stream[edge_idx_offset].read(read_value_stream[edge_idx_offset]); }
			// Fetch indexes of message arrays.
			for (IDX_TYPE edge_idx_offset=0; edge_idx_offset<EDGE_UNI_PARA; edge_idx_offset++) {
				# pragma HLS UNROLL
				if (!edge_idx_overflow[edge_idx_offset]) {
				src_idx[edge_idx_offset] = g_adj_list_scatter[edge_idx_offset][edge_idx[edge_idx_offset]][0];
				dst_idx[edge_idx_offset] = g_adj_list_scatter[edge_idx_offset][edge_idx[edge_idx_offset]][1]; }}
			// Fetch flags of message arrays.
			for (IDX_TYPE edge_idx_offset=0; edge_idx_offset<EDGE_UNI_PARA; edge_idx_offset++) {
				# pragma HLS UNROLL
				if (!edge_idx_overflow[edge_idx_offset]) {
				src_flag[edge_idx_offset] = msg_src_new_flag[edge_idx_offset][src_idx[edge_idx_offset]];
				dst_flag[edge_idx_offset] = msg_dst_new_flag[edge_idx_offset][dst_idx[edge_idx_offset]]; }}
			// Update flags of message arrays.
			for (IDX_TYPE edge_idx_offset=0; edge_idx_offset<EDGE_UNI_PARA; edge_idx_offset++) {
				# pragma HLS UNROLL
				if (!edge_idx_overflow[edge_idx_offset]) {
				msg_src_new_flag[edge_idx_offset][src_idx[edge_idx_offset]] =
						(src_flag[edge_idx_offset] == (IDX_TYPE)-1)? (IDX_TYPE)0:
						((src_flag[edge_idx_offset] == (IDX_TYPE)0)? (IDX_TYPE)1: (IDX_TYPE)0);
				msg_dst_new_flag[edge_idx_offset][dst_idx[edge_idx_offset]] =
						(dst_flag[edge_idx_offset] == (IDX_TYPE)-1)? (IDX_TYPE)0:
						((dst_flag[edge_idx_offset] == (IDX_TYPE)0)? (IDX_TYPE)1: (IDX_TYPE)0); }}}

		// Scatter messages.
		for (IDX_TYPE edge_idx_offset=0; edge_idx_offset<EDGE_UNI_PARA; edge_idx_offset++) {
			# pragma HLS UNROLL

			// Fetch node embeddings from the FIFO.
			EDGE_TYPE edge_cache;
			# pragma HLS AGGREGATE variable=edge_cache
			msg_stream[edge_idx_offset].read(edge_cache);
			DT1_PARA src_emb = edge_cache.node_emb_dst;
			DT1_PARA dst_emb = edge_cache.node_emb_src;
			# pragma HLS AGGREGATE variable=src_emb
			# pragma HLS AGGREGATE variable=dst_emb
			DT1_TYPE value = read_value_stream[edge_idx_offset];

			// Scatter messages.
			if (!edge_idx_overflow[edge_idx_offset])
			for (IDX_TYPE dim_offset=0; dim_offset<DIM_PARA; dim_offset++) {
				# pragma HLS UNROLL
				// Get the output dimension.
				IDX_TYPE dim = dim_base + dim_offset;
				BOL_TYPE dim_overflow = (dim >= (IDX_TYPE)NODE_EMB_DIM);
				// Fetch old messages.
				DT0_TYPE old_msg_src =
						((src_flag[edge_idx_offset] == (IDX_TYPE)-1) || dim_overflow)? (DT0_TYPE)0.0:
						((src_flag[edge_idx_offset] == (IDX_TYPE)0)?
						msg_src_new[0][edge_idx_offset][src_idx[edge_idx_offset]][dim]:
						msg_src_new[1][edge_idx_offset][src_idx[edge_idx_offset]][dim]);
				DT0_TYPE old_msg_dst =
						((dst_flag[edge_idx_offset] == (IDX_TYPE)-1) || dim_overflow)? (DT0_TYPE)0.0:
						((dst_flag[edge_idx_offset] == (IDX_TYPE)0)?
						msg_dst_new[0][edge_idx_offset][dst_idx[edge_idx_offset]][dim]:
						msg_dst_new[1][edge_idx_offset][dst_idx[edge_idx_offset]][dim]);
				// Calculate new messages.
				DT0_TYPE new_msg_src = old_msg_src + src_emb[dim_offset] * value;
				DT0_TYPE new_msg_dst = old_msg_dst + dst_emb[dim_offset] * value;
				// Save new messages.
				if (!dim_overflow) {
					((src_flag[edge_idx_offset] == (IDX_TYPE)0)?
							msg_src_new[1][edge_idx_offset][src_idx[edge_idx_offset]][dim]:
							msg_src_new[0][edge_idx_offset][src_idx[edge_idx_offset]][dim])
									= new_msg_src;
					((dst_flag[edge_idx_offset] == (IDX_TYPE)0)?
							msg_dst_new[1][edge_idx_offset][dst_idx[edge_idx_offset]][dim]:
							msg_dst_new[0][edge_idx_offset][dst_idx[edge_idx_offset]][dim])
									= new_msg_dst; }}}
	}
}



// **************************** finalize.cpp **************************** //

/*

	This part is to implement the finalization function. It first
	receives data from edge-processing network, and then sends them
	into the output interfaces.

	Version: 0.1
		Author: Hanqing Zhang
		Email : hanqing.zhang@zju.edu.cn
		Date  : July-10-2023
		Commit: Add "finalize" function.

	Version: 0.2
		Author: Hanqing Zhang
		Email : hanqing.zhang@zju.edu.cn
		Date  : July-30-2023
		Commit: Modify the model because of data type changing.

	Function:
		void finalize(
			IDX_TYPE num_of_edges,
			hls::stream<DT1_TYPE> value_stream[EDGE_UNI_PARA],
			OUTPUT_RESULT_TYPE result[OUTPUT_RESULT_COUNT]);
		Input arguments:
			num_of_edges: the number of edges.
			value_stream: attention values from edge-processing network.
		Output arguments:
			result: the final results.

*/

static void finalize(
		IDX_TYPE num_of_edges,
		hls::stream<DT1_TYPE> value_stream[EDGE_UNI_PARA],
		OUTPUT_RESULT_TYPE result[OUTPUT_RESULT_COUNT])
{
	# pragma HLS INLINE off

	// Process EDGE_UNI_PARA edges each time.
	IDX_TYPE iter_uni_num = ceildiv(num_of_edges, (IDX_TYPE)EDGE_UNI_PARA);
	for (	IDX_TYPE iter_uni_cnt=(IDX_TYPE)0;
			iter_uni_cnt<iter_uni_num;
			iter_uni_cnt+=(IDX_TYPE)1)
	{
		# pragma HLS LOOP_TRIPCOUNT min=EDGE_UNI_COUNT max=EDGE_UNI_COUNT avg=EDGE_UNI_COUNT
		# pragma HLS PIPELINE II=1 style=frp

		// Fetch results from FIFOs.
		OUTPUT_RESULT_TYPE result_cache;
		# pragma HLS AGGREGATE variable=result_cache
		for (IDX_TYPE edge_idx_offset=0; edge_idx_offset<EDGE_UNI_PARA; edge_idx_offset++)
			# pragma HLS UNROLL
			value_stream[edge_idx_offset].read(result_cache[edge_idx_offset]);
		// Output results.
		result[iter_uni_cnt] = result_cache;
	}
}



// **************************** bypass.cpp **************************** //

/*

	This part is to implement bypass channels in the top function.
	The channel can send new node embeddings into the last block.

	Version: 0.1
		Author: Hanqing Zhang
		Email : hanqing.zhang@zju.edu.cn
		Date  : July-10-2023
		Commit: Add function "bypass".

	Function:
		void bypass(
			BOL_TYPE last_layer,
			IDX_TYPE num_of_nodes,
			hls::stream<DT1_PARA> node_fifo_input[NODE_UNI_PARA],
			hls::stream<DT1_PARA> node_fifo_output[NODE_UNI_PARA]);
		Input arguments:
			last_layer: flag of the last layer.
			num_of_nodes: the number of nodes.
			node_fifo_input: node embeddings from the previous block.
		Output arguments:
			node_fifo_output: node embeddings to the next block.

*/

static void bypass(
		BOL_TYPE last_layer,
		IDX_TYPE num_of_nodes,
		hls::stream<DT1_PARA> node_fifo_input[NODE_UNI_PARA],
		hls::stream<DT1_PARA> node_fifo_output[NODE_UNI_PARA])
{
	# pragma HLS INLINE off

	// If it isn't the last layer
	if (!last_layer)
	{
		// Get NODE_UNI_PARA nodes each time.
		// Get DIM_PARA dimensions each time.
		IDX_TYPE iter_uni_num = ceildiv(num_of_nodes, (IDX_TYPE)NODE_UNI_PARA);
		IDX_TYPE iter_dim_num = ceildiv((IDX_TYPE)NODE_EMB_DIM, (IDX_TYPE)DIM_PARA);
		for (	IDX_TYPE iter_uni_cnt=(IDX_TYPE)0,
				node_idx_base=(IDX_TYPE)0;
				iter_uni_cnt<iter_uni_num;
				iter_uni_cnt+=(IDX_TYPE)1,
				node_idx_base+=(IDX_TYPE)NODE_UNI_PARA)
		# pragma HLS LOOP_TRIPCOUNT min=NODE_UNI_COUNT max=NODE_UNI_COUNT avg=NODE_UNI_COUNT
		for (	IDX_TYPE iter_dim_cnt=(IDX_TYPE)0;
				iter_dim_cnt<iter_dim_num;
				iter_dim_cnt+=(IDX_TYPE)1)
		{
			# pragma HLS LOOP_TRIPCOUNT min=DIM_COUNT max=DIM_COUNT avg=DIM_COUNT
			# pragma HLS PIPELINE II=1 style=frp

			// Fetch node embeddings from the FIFO and then send them to the next block.
			for (IDX_TYPE node_idx_offset=0; node_idx_offset<NODE_UNI_PARA; node_idx_offset++) {
				# pragma HLS UNROLL
				DT1_PARA node_stream;
				# pragma HLS AGGREGATE variable=node_stream
				node_fifo_input[node_idx_offset].read(node_stream);
				node_fifo_output[node_idx_offset].write(node_stream); }
		}
	}
}



// **************************** Top Modules **************************** //

/*

	This part contains all the top modules, including the node-
	processing module, the adapter module, the edge-processing
	module, and the message-passing module. These four modules
	can be combined to form the layer-computation module. By
	invoking the layer-computation module, the kernel can send
	and receive data between ping-pong RAMs.

	Version: 0.1
		Author: Hanqing Zhang
		Email : hanqing.zhang@zju.edu.cn
		Date  : July-14-2023
		Commit: Add four sub-module functions.

	Version: 0.2
		Author: Hanqing Zhang
		Email : hanqing.zhang@zju.edu.cn
		Date  : July-20-2023
		Commit: Add top function.

	Version: 0.3
		Author: Hanqing Zhang
		Email : hanqing.zhang@zju.edu.cn
		Date  : August-3-2023
		Commit: Add bypass channel.

	Version: 0.4
		Author: Hanqing Zhang
		Email : hanqing.zhang@zju.edu.cn
		Date  : September-1-2023
		Commit: Change the design for multiple graphs calculating.

	Function:
		void compute_layer(
			// Flags.
			BOL_TYPE first_layer,
			BOL_TYPE last_layer,
			// Numbers.
			IDX_TYPE num_of_nodes,
			IDX_TYPE num_of_edges,
			// Input data.
			DT0_TYPE msg_dst_old[2][EDGE_UNI_PARA][NODE_NUM_MAX][NODE_EMB_DIM],
			DT0_TYPE msg_src_old[2][EDGE_UNI_PARA][NODE_NUM_MAX][NODE_EMB_DIM],
			IDX_TYPE msg_dst_old_flag[EDGE_UNI_PARA][NODE_NUM_MAX],
			IDX_TYPE msg_src_old_flag[EDGE_UNI_PARA][NODE_NUM_MAX],
			DT1_TYPE node_emb_old[NODE_NUM_MAX][NODE_EMB_DIM],
			// Output data.
			DT0_TYPE msg_dst_new[2][EDGE_UNI_PARA][NODE_NUM_MAX][NODE_EMB_DIM],
			DT0_TYPE msg_src_new[2][EDGE_UNI_PARA][NODE_NUM_MAX][NODE_EMB_DIM],
			IDX_TYPE msg_dst_new_flag[EDGE_UNI_PARA][NODE_NUM_MAX],
			IDX_TYPE msg_src_new_flag[EDGE_UNI_PARA][NODE_NUM_MAX],
			DT1_TYPE node_emb_new[NODE_NUM_MAX][NODE_EMB_DIM],
			// Results.
			OUTPUT_RESULT_TYPE result[OUTPUT_RESULT_COUNT])
		Input arguments:
			first_layer: flag of the first layer.
			last_layer: flag of the last layer.
			num_of_nodes: the number of nodes.
			num_of_edges: the number of edges.
			msg_dst_old: old source message array.
			msg_src_old: old destination message array.
			msg_dst_old_flag: flags of where the old source messages are.
			msg_src_old_flag: flags of where the old destination messages are.
			node_emb_old: old node embedding array.
		Output arguments:
			msg_dst_new: new source message array.
			msg_src_new: new destination message array.
			msg_dst_new_flag: flags of where the new source messages are.
			msg_src_new_flag: flags of where the new destination messages are.
			node_emb_new: new node embedding array.
			result: the final results to the m_axi interface.

*/

static void node_embedding(
		// Flags.
		BOL_TYPE first_layer,
		BOL_TYPE last_layer,
		// Numbers.
		IDX_TYPE num_of_nodes,
		// Input data.
		DT1_TYPE node_emb_old[NODE_NUM_MAX][NODE_EMB_DIM],
		DT0_TYPE msg_dst_old[2][EDGE_UNI_PARA][NODE_NUM_MAX][NODE_EMB_DIM],
		DT0_TYPE msg_src_old[2][EDGE_UNI_PARA][NODE_NUM_MAX][NODE_EMB_DIM],
		IDX_TYPE msg_dst_old_flag[EDGE_UNI_PARA][NODE_NUM_MAX],
		IDX_TYPE msg_src_old_flag[EDGE_UNI_PARA][NODE_NUM_MAX],
		// Output data.
		hls::stream<DT1_PARA> node_fifo[NODE_UNI_PARA],
		// Bypass data.
		hls::stream<DT1_PARA> node_fifo_bypass_out[NODE_UNI_PARA])
{
	# pragma HLS INLINE off

	// If it is the first layer.
	// Activate the input network.
	if (first_layer)
		input_nw(
				last_layer,
				num_of_nodes,
				node_fifo,
				node_fifo_bypass_out);
	// If not the first layer.
	// Activate the node-processing network.
	else
		node_nw(
				last_layer,
				num_of_nodes,
				node_emb_old,
				msg_dst_old,
				msg_src_old,
				msg_dst_old_flag,
				msg_src_old_flag,
				node_fifo,
				node_fifo_bypass_out);
}

static void adapter_top(
		// Flags.
		BOL_TYPE last_layer,
		// Numbers.
		IDX_TYPE num_of_nodes,
		IDX_TYPE num_of_edges,
		// Input data.
		hls::stream<DT1_PARA> node_fifo[NODE_UNI_PARA],
		// Output data.
		hls::stream<EDGE_TYPE> edge_fifo[EDGE_UNI_PARA],
		// Bypass data.
		hls::stream<DT1_PARA> node_fifo_bypass_in[NODE_UNI_PARA],
		hls::stream<DT1_PARA> node_fifo_bypass_out[NODE_UNI_PARA])
{
	# pragma HLS INLINE off

	// Activate the adapter.
	adapter(
			num_of_nodes,
			num_of_edges,
			node_fifo,
			edge_fifo);
	// Bypass channel.
	bypass(
			last_layer,
			num_of_nodes,
			node_fifo_bypass_in,
			node_fifo_bypass_out);
}

static void edge_embedding(
		// Flags.
		BOL_TYPE last_layer,
		// Numbers.
		IDX_TYPE num_of_nodes,
		IDX_TYPE num_of_edges,
		// Input data.
		hls::stream<EDGE_TYPE> edge_stream[EDGE_UNI_PARA],
		// Output data.
		hls::stream<DT1_TYPE> value_stream[EDGE_UNI_PARA],
		hls::stream<EDGE_TYPE> msg_stream[EDGE_UNI_PARA],
		// Bypass data.
		hls::stream<DT1_PARA> node_fifo_bypass_in[NODE_UNI_PARA],
		hls::stream<DT1_PARA> node_fifo_bypass_out[NODE_UNI_PARA])
{
	# pragma HLS INLINE off

	// Activate the edge-processing network.
	edge_nw(
			last_layer,
			num_of_edges,
			edge_stream,
			value_stream,
			msg_stream);
	// Bypass channel.
	bypass(
			last_layer,
			num_of_nodes,
			node_fifo_bypass_in,
			node_fifo_bypass_out);
}

static void message_passing(
		// Flags.
		BOL_TYPE last_layer,
		// Numbers.
		IDX_TYPE num_of_nodes,
		IDX_TYPE num_of_edges,
		// Input data.
		hls::stream<DT1_TYPE> value_stream[EDGE_UNI_PARA],
		hls::stream<EDGE_TYPE> msg_stream[EDGE_UNI_PARA],
		hls::stream<DT1_PARA> node_stream[NODE_UNI_PARA],
		// Output data.
		DT0_TYPE msg_src_new[2][EDGE_UNI_PARA][NODE_NUM_MAX][NODE_EMB_DIM],
		DT0_TYPE msg_dst_new[2][EDGE_UNI_PARA][NODE_NUM_MAX][NODE_EMB_DIM],
		IDX_TYPE msg_src_new_flag[EDGE_UNI_PARA][NODE_NUM_MAX],
		IDX_TYPE msg_dst_new_flag[EDGE_UNI_PARA][NODE_NUM_MAX],
		DT1_TYPE node_emb_new[NODE_NUM_MAX][NODE_EMB_DIM],
		// Results.
		OUTPUT_RESULT_TYPE result[OUTPUT_RESULT_COUNT])
{
	# pragma HLS INLINE off

	// If it is the last layer.
	// Activation finalization function.
	if (last_layer)
		finalize(
				num_of_edges,
				value_stream,
				result);
	// If it is not the last layer.
	// Activate node and message scatter.
	else {
		node_scatter(
				num_of_nodes,
				node_stream,
				node_emb_new);
		msg_scatter(
				num_of_edges,
				value_stream,
				msg_stream,
				msg_src_new,
				msg_dst_new,
				msg_src_new_flag,
				msg_dst_new_flag); }
}

static void compute_layer(
		// Flags.
		BOL_TYPE first_layer,
		BOL_TYPE last_layer,
		// Numbers.
		IDX_TYPE num_of_nodes,
		IDX_TYPE num_of_edges,
		// Input data.
		DT0_TYPE msg_dst_old[2][EDGE_UNI_PARA][NODE_NUM_MAX][NODE_EMB_DIM],
		DT0_TYPE msg_src_old[2][EDGE_UNI_PARA][NODE_NUM_MAX][NODE_EMB_DIM],
		IDX_TYPE msg_dst_old_flag[EDGE_UNI_PARA][NODE_NUM_MAX],
		IDX_TYPE msg_src_old_flag[EDGE_UNI_PARA][NODE_NUM_MAX],
		DT1_TYPE node_emb_old[NODE_NUM_MAX][NODE_EMB_DIM],
		// Output data.
		DT0_TYPE msg_dst_new[2][EDGE_UNI_PARA][NODE_NUM_MAX][NODE_EMB_DIM],
		DT0_TYPE msg_src_new[2][EDGE_UNI_PARA][NODE_NUM_MAX][NODE_EMB_DIM],
		IDX_TYPE msg_dst_new_flag[EDGE_UNI_PARA][NODE_NUM_MAX],
		IDX_TYPE msg_src_new_flag[EDGE_UNI_PARA][NODE_NUM_MAX],
		DT1_TYPE node_emb_new[NODE_NUM_MAX][NODE_EMB_DIM],
		// Results.
		OUTPUT_RESULT_TYPE result[OUTPUT_RESULT_COUNT])
{
	# pragma HLS INLINE off
	# pragma HLS DATAFLOW

	// FIFO declarations.
	hls::stream<DT1_PARA> node_fifo[NODE_UNI_PARA];
	hls::stream<EDGE_TYPE> edge_fifo[EDGE_UNI_PARA];
	hls::stream<DT1_TYPE> value_fifo[EDGE_UNI_PARA];
	hls::stream<EDGE_TYPE> msg_fifo[EDGE_UNI_PARA];
	hls::stream<DT1_PARA> node_fifo_bypass_0[NODE_UNI_PARA];
	hls::stream<DT1_PARA> node_fifo_bypass_1[NODE_UNI_PARA];
	hls::stream<DT1_PARA> node_fifo_bypass_2[NODE_UNI_PARA];
	# pragma HLS STREAM depth=FIFO_DEPTH variable=node_fifo
	# pragma HLS STREAM depth=FIFO_DEPTH variable=edge_fifo
	# pragma HLS STREAM depth=FIFO_DEPTH variable=value_fifo
	# pragma HLS STREAM depth=(FIFO_DEPTH*DIM_COUNT) variable=msg_fifo
	# pragma HLS STREAM depth=FIFO_DEPTH variable=node_fifo_bypass_0
	# pragma HLS STREAM depth=FIFO_DEPTH variable=node_fifo_bypass_1
	# pragma HLS STREAM depth=FIFO_DEPTH variable=node_fifo_bypass_2

	// Invoke functions.
	node_embedding(
			first_layer, 
			last_layer,
			num_of_nodes,
			node_emb_old, 
			msg_dst_old, 
			msg_src_old, 
			msg_dst_old_flag, 
			msg_src_old_flag,
			node_fifo,
			node_fifo_bypass_0);
	adapter_top(
			last_layer,
			num_of_nodes, 
			num_of_edges,
			node_fifo,
			edge_fifo,
			node_fifo_bypass_0,
			node_fifo_bypass_1);
	edge_embedding(
			last_layer,
			num_of_nodes, 
			num_of_edges,
			edge_fifo,
			value_fifo, msg_fifo,
			node_fifo_bypass_1,
			node_fifo_bypass_2);
	message_passing(
			last_layer,
			num_of_nodes, 
			num_of_edges,
			value_fifo, 
			msg_fifo,
			node_fifo_bypass_2,
			msg_src_new, 
			msg_dst_new, 
			msg_src_new_flag, 
			msg_dst_new_flag, 
			node_emb_new,
			result);
}

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
		OUTPUT_RESULT_TYPE result[GRAPH_NUM_MAX][OUTPUT_RESULT_COUNT])
{
	// Interface setting.
	# pragma HLS INTERFACE s_axilite port=return
	# pragma HLS INTERFACE m_axi port=num_of_nodes depth=(FIFO_DEPTH) offset=slave bundle=mem
	# pragma HLS INTERFACE m_axi port=num_of_edges depth=(FIFO_DEPTH) offset=slave bundle=mem
	# pragma HLS INTERFACE m_axi port=node_feature depth=(FIFO_DEPTH*INPUT_FEA_COUNT) offset=slave bundle=mem
	# pragma HLS INTERFACE m_axi port=adj_list depth=(FIFO_DEPTH*INPUT_ADJ_COUNT) offset=slave bundle=mem
	# pragma HLS INTERFACE m_axi port=result depth=(FIFO_DEPTH*OUTPUT_RESULT_COUNT) offset=slave bundle=mem

	// Global array partition.
	# pragma HLS ARRAY_PARTITION dim=1 type=cyclic factor=NODE_UNI_PARA variable=g_node_feature
	# pragma HLS ARRAY_PARTITION dim=2 type=complete variable=g_node_feature
	# pragma HLS ARRAY_PARTITION dim=1 type=complete variable=g_adj_list_adapter
	# pragma HLS ARRAY_PARTITION dim=3 type=complete variable=g_adj_list_adapter
	# pragma HLS ARRAY_PARTITION dim=1 type=complete variable=g_adj_list_scatter
	# pragma HLS ARRAY_PARTITION dim=3 type=complete variable=g_adj_list_scatter
	# pragma HLS ARRAY_PARTITION dim=1 type=cyclic factor=NODE_UNI_PARA variable=g_node_emb_ping
	# pragma HLS ARRAY_PARTITION dim=1 type=cyclic factor=NODE_UNI_PARA variable=g_node_emb_pong
	# pragma HLS ARRAY_PARTITION dim=2 type=cyclic factor=DIM_PARA variable=g_node_emb_ping
	# pragma HLS ARRAY_PARTITION dim=2 type=cyclic factor=DIM_PARA variable=g_node_emb_pong
	# pragma HLS ARRAY_PARTITION dim=1 type=complete variable=g_msg_dst_ping
	# pragma HLS ARRAY_PARTITION dim=1 type=complete variable=g_msg_dst_pong
	# pragma HLS ARRAY_PARTITION dim=1 type=complete variable=g_msg_src_ping
	# pragma HLS ARRAY_PARTITION dim=1 type=complete variable=g_msg_src_pong
	# pragma HLS ARRAY_PARTITION dim=2 type=complete variable=g_msg_dst_ping
	# pragma HLS ARRAY_PARTITION dim=2 type=complete variable=g_msg_dst_pong
	# pragma HLS ARRAY_PARTITION dim=2 type=complete variable=g_msg_src_ping
	# pragma HLS ARRAY_PARTITION dim=2 type=complete variable=g_msg_src_pong
	# pragma HLS ARRAY_PARTITION dim=3 type=cyclic factor=NODE_UNI_PARA variable=g_msg_dst_ping
	# pragma HLS ARRAY_PARTITION dim=3 type=cyclic factor=NODE_UNI_PARA variable=g_msg_dst_pong
	# pragma HLS ARRAY_PARTITION dim=3 type=cyclic factor=NODE_UNI_PARA variable=g_msg_src_ping
	# pragma HLS ARRAY_PARTITION dim=3 type=cyclic factor=NODE_UNI_PARA variable=g_msg_src_pong
	# pragma HLS ARRAY_PARTITION dim=4 type=cyclic factor=DIM_PARA variable=g_msg_dst_ping
	# pragma HLS ARRAY_PARTITION dim=4 type=cyclic factor=DIM_PARA variable=g_msg_dst_pong
	# pragma HLS ARRAY_PARTITION dim=4 type=cyclic factor=DIM_PARA variable=g_msg_src_ping
	# pragma HLS ARRAY_PARTITION dim=4 type=cyclic factor=DIM_PARA variable=g_msg_src_pong
	# pragma HLS ARRAY_PARTITION dim=1 type=complete variable=g_msg_dst_ping_flag
	# pragma HLS ARRAY_PARTITION dim=1 type=complete variable=g_msg_dst_pong_flag
	# pragma HLS ARRAY_PARTITION dim=1 type=complete variable=g_msg_src_ping_flag
	# pragma HLS ARRAY_PARTITION dim=1 type=complete variable=g_msg_src_pong_flag
	# pragma HLS ARRAY_PARTITION dim=2 type=cyclic factor=NODE_UNI_PARA variable=g_msg_dst_ping_flag
	# pragma HLS ARRAY_PARTITION dim=2 type=cyclic factor=NODE_UNI_PARA variable=g_msg_dst_pong_flag
	# pragma HLS ARRAY_PARTITION dim=2 type=cyclic factor=NODE_UNI_PARA variable=g_msg_src_ping_flag
	# pragma HLS ARRAY_PARTITION dim=2 type=cyclic factor=NODE_UNI_PARA variable=g_msg_src_pong_flag

	// Instance limit.
	# pragma HLS ALLOCATION function instances=load_graph limit=1
	# pragma HLS ALLOCATION function instances=compute_layer limit=1

	// Start process.
	for (IDX_TYPE graph_idx=(IDX_TYPE)0; graph_idx<num_of_graphs; graph_idx+=(IDX_TYPE)1)
	{
		# pragma HLS LOOP_TRIPCOUNT min=GRAPH_NUM_MAX max=GRAPH_NUM_MAX avg=GRAPH_NUM_MAX

		// Load one graph.
		IDX_TYPE num_of_nodes_cache = num_of_nodes[graph_idx];
		IDX_TYPE num_of_edges_cache = num_of_edges[graph_idx];
		num_of_nodes_cache = load_graph(num_of_nodes_cache, num_of_edges_cache, node_feature[graph_idx], adj_list[graph_idx]);

		// Calculate one graph.
		for (IDX_TYPE layer=(IDX_TYPE)0; layer<(IDX_TYPE)(NUM_OF_LAYERS+1); layer+=(IDX_TYPE)1)
		{
			# pragma HLS LOOP_TRIPCOUNT min=(NUM_OF_LAYERS+1) max=(NUM_OF_LAYERS+1) avg=(NUM_OF_LAYERS+1)

			// Flags.
			BOL_TYPE first_layer = (layer == (IDX_TYPE)0);
			BOL_TYPE last_layer = (layer == (IDX_TYPE)NUM_OF_LAYERS);
			BOL_TYPE even_or_odd = (layer % (IDX_TYPE)2 == (IDX_TYPE)0);

			// Invoke.
			if (even_or_odd)
			compute_layer(
					first_layer, last_layer,
					num_of_nodes_cache, num_of_edges_cache,
					g_msg_dst_ping, g_msg_src_ping, g_msg_dst_ping_flag, g_msg_src_ping_flag, g_node_emb_ping,
					g_msg_dst_pong, g_msg_src_pong, g_msg_dst_pong_flag, g_msg_src_pong_flag, g_node_emb_pong,
					result[graph_idx]);
			else
			compute_layer(
					first_layer, last_layer,
					num_of_nodes_cache, num_of_edges_cache,
					g_msg_dst_pong, g_msg_src_pong, g_msg_dst_pong_flag, g_msg_src_pong_flag, g_node_emb_pong,
					g_msg_dst_ping, g_msg_src_ping, g_msg_dst_ping_flag, g_msg_src_ping_flag, g_node_emb_ping,
					result[graph_idx]);
		}
	}
}
}



