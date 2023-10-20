
// Includes.
# include "config.h"
# include "fstream"
# include "vector"
# ifdef EMU
# include "xcl2.hpp"
# endif

// Global Variables.
# ifdef EMU
int TEST_TIMES = 1;
IDX_TYPE host_num_of_graphs = 10;
std::vector<IDX_TYPE, aligned_allocator<IDX_TYPE> > host_num_of_nodes(GRAPH_NUM_MAX);
std::vector<IDX_TYPE, aligned_allocator<IDX_TYPE> > host_num_of_edges(GRAPH_NUM_MAX);
std::vector<INPUT_FEA_TYPE, aligned_allocator<INPUT_FEA_TYPE> > host_node_feature(GRAPH_NUM_MAX*INPUT_FEA_COUNT);
std::vector<INPUT_ADJ_TYPE, aligned_allocator<INPUT_ADJ_TYPE> > host_adj_list(GRAPH_NUM_MAX*INPUT_ADJ_COUNT);
std::vector<OUTPUT_RESULT_TYPE, aligned_allocator<OUTPUT_RESULT_TYPE> > host_result_read(GRAPH_NUM_MAX*OUTPUT_RESULT_COUNT);
std::vector<OUTPUT_RESULT_TYPE, aligned_allocator<OUTPUT_RESULT_TYPE> > host_result_write(GRAPH_NUM_MAX*OUTPUT_RESULT_COUNT);
# endif
# ifdef SIM
IDX_TYPE tb_num_of_graphs = 10;
IDX_TYPE tb_num_of_nodes[GRAPH_NUM_MAX];
IDX_TYPE tb_num_of_edges[GRAPH_NUM_MAX];
INPUT_FEA_TYPE tb_node_feature[GRAPH_NUM_MAX][INPUT_FEA_COUNT];
INPUT_ADJ_TYPE tb_adj_list[GRAPH_NUM_MAX][INPUT_ADJ_COUNT];
OUTPUT_RESULT_TYPE tb_result_read[GRAPH_NUM_MAX][OUTPUT_RESULT_COUNT];
OUTPUT_RESULT_TYPE tb_result_write[GRAPH_NUM_MAX][OUTPUT_RESULT_COUNT];
# endif

static void load_graph();

# ifdef EMU
int main(int argc, char** argv)
{
    if (argc != 2)
	return EXIT_FAILURE;

	// **************** Declare data. **************** //
    cl_int err;
    cl::Context context;
    cl::Kernel kernel_compute_graph;
    cl::CommandQueue q;

    // **************** Write the kernel. **************** //
    // Get the binary file.
    std::string binaryFile = argv[1];
    auto fileBuf = xcl::read_binary_file(binaryFile);
    cl::Program::Binaries bins{{fileBuf.data(), fileBuf.size()}};
    // Get devices.
    auto devices = xcl::get_xil_devices();
    // Select a valid device.
    bool valid_device = false;
    for (unsigned int i=0; i<devices.size(); i++)
    {
    	// Select a device from the devices list.
        auto device = devices[i];
        OCL_CHECK(err, context = cl::Context(device,
        		nullptr, nullptr, nullptr, &err));
        OCL_CHECK(err, q = cl::CommandQueue(context,
        		device, CL_QUEUE_PROFILING_ENABLE, &err));
        // Write the binary file into the selected device.
        cl::Program program(context, {device}, bins, nullptr, &err);
        // Whether the device is valid.
        if (err == CL_SUCCESS)
        {
        	// Print INFO.
            std::cout << "Device[" << i << "]: program successful!" << std::endl;
            // Get the kernel.
            OCL_CHECK(err, kernel_compute_graph = cl::Kernel(program, "kernel_compute_graph", &err));
            valid_device = true;
            break;
        }
    }
    // Whether all the devices are invalid.
    if (!valid_device)
    {
        std::cout << "Failed to program any device found!" << std::endl;
        exit(EXIT_FAILURE);
    }

    // **************** Load data. **************** //
	// Load data from files.
	load_graph();
	// Creater buffers.
	OCL_CHECK(err, cl::Buffer buffer_num_of_nodes(context,
			CL_MEM_USE_HOST_PTR | CL_MEM_READ_ONLY,
			GRAPH_NUM_MAX * sizeof(IDX_TYPE),
			host_num_of_nodes.data(), &err));
	OCL_CHECK(err, cl::Buffer buffer_num_of_edges(context,
			CL_MEM_USE_HOST_PTR | CL_MEM_READ_ONLY,
			GRAPH_NUM_MAX * sizeof(IDX_TYPE),
			host_num_of_edges.data(), &err));
	OCL_CHECK(err, cl::Buffer buffer_node_feature(context,
			CL_MEM_USE_HOST_PTR | CL_MEM_READ_ONLY,
			GRAPH_NUM_MAX * INPUT_FEA_COUNT * sizeof(INPUT_FEA_TYPE),
			host_node_feature.data(), &err));
	OCL_CHECK(err, cl::Buffer buffer_adj_list(context,
			CL_MEM_USE_HOST_PTR | CL_MEM_READ_ONLY,
			GRAPH_NUM_MAX * INPUT_ADJ_COUNT * sizeof(INPUT_ADJ_TYPE),
			host_adj_list.data(), &err));
	OCL_CHECK(err, cl::Buffer buffer_result(context,
			CL_MEM_USE_HOST_PTR | CL_MEM_WRITE_ONLY,
			GRAPH_NUM_MAX * OUTPUT_RESULT_COUNT * sizeof(OUTPUT_RESULT_TYPE),
			host_result_write.data(), &err));

    // **************** Process data. **************** //
	// Set arguments of the kernel.
	int idx = 0;
	OCL_CHECK(err, err = kernel_compute_graph.setArg(idx++, host_num_of_graphs));
	OCL_CHECK(err, err = kernel_compute_graph.setArg(idx++, buffer_num_of_nodes));
	OCL_CHECK(err, err = kernel_compute_graph.setArg(idx++, buffer_num_of_edges));
	OCL_CHECK(err, err = kernel_compute_graph.setArg(idx++, buffer_node_feature));
	OCL_CHECK(err, err = kernel_compute_graph.setArg(idx++, buffer_adj_list));
	OCL_CHECK(err, err = kernel_compute_graph.setArg(idx++, buffer_result));
	// Launch the kernel.
	for (int test_idx=0; test_idx<TEST_TIMES; test_idx++)
	{
		OCL_CHECK(err, err = q.enqueueMigrateMemObjects({
			buffer_num_of_nodes,
			buffer_num_of_edges,
			buffer_node_feature,
			buffer_adj_list
			}, 0));
		OCL_CHECK(err, err = q.enqueueTask(
			kernel_compute_graph));
		OCL_CHECK(err, err = q.enqueueMigrateMemObjects({
			buffer_result
			}, CL_MIGRATE_MEM_OBJECT_HOST));
		OCL_CHECK(err, err = q.finish());
	}

    // **************** Test results. **************** //
    // Calculate the total number of edges.
	int total_edges = 0;
    for (int graph_idx=0; graph_idx<host_num_of_graphs; graph_idx++)
		total_edges += host_num_of_edges[graph_idx];
	// Calculate errors.
	double error_sum = 0.0;
	double error_2_sum = 0.0;
	double error_max = 0.0;
    for (int graph_idx=0; graph_idx<host_num_of_graphs; graph_idx++)
	for (int i=0; i<host_num_of_edges[graph_idx]; i++)
	{
		int idx_0 = graph_idx*OUTPUT_RESULT_COUNT + i/EDGE_UNI_PARA;
		int idx_1 = i%EDGE_UNI_PARA;
		double result_read = (double)host_result_read[idx_0][idx_1];
		double result_write = (double)host_result_write[idx_0][idx_1];
		double result_error = std::abs(result_read - result_write);
		error_sum += result_error;
		error_2_sum += result_error * result_error;
		error_max = (error_max > result_error)? error_max: result_error;
	}
	// Print results.
	std::cout.precision(4);
	std::cout << "There are "<< host_num_of_graphs << " graphs input. " << std::endl;
	std::cout << "Maximum error: " << error_max << std::endl;
	std::cout << "Mean error: " << error_sum / (double)total_edges << std::endl;
	std::cout << "Standard deviation: " << std::sqrt(error_2_sum / (double)total_edges) << std::endl;

    // **************** Finish processing. **************** //
	std::cout << "All done! " << std::endl;
    return EXIT_SUCCESS;
}
# endif

# ifdef SIM
int main(void)
{
	// Load data from files.
	load_graph();

	// Compute one graph.
	kernel_compute_graph(
			tb_num_of_graphs,
			tb_num_of_nodes,
			tb_num_of_edges,
			tb_node_feature,
			tb_adj_list,
			tb_result_write);

	// Test results.
	for (int graph_idx=0; graph_idx<tb_num_of_graphs; graph_idx++)
	{
		double error_sum = 0.0;
		double error_max = 0.0;
		for (int i=0; i<tb_num_of_edges[graph_idx]; i++)
		{
			double result_read = (double)tb_result_read[graph_idx][i/EDGE_UNI_PARA][i%EDGE_UNI_PARA];
			double result_write = (double)tb_result_write[graph_idx][i/EDGE_UNI_PARA][i%EDGE_UNI_PARA];
			double result_error = std::abs(result_read - result_write);
			error_sum += result_error;
			error_max = (error_max > result_error)? error_max: result_error;
		}
		std::cout << "Graph " << graph_idx << ": " << std::endl;
		std::cout << "    Mean error: " << error_sum / tb_num_of_edges[graph_idx] << std::endl;
		std::cout << "    Max error: " << error_max << std::endl;
		std::cout << " " << std::endl;
    }

	// Finish processing.
    return 0;
}
# endif

static void load_graph()
{
	// Data path.
	# ifdef EMU
	std::string path = "C:/Users/zju/Desktop/TrackGNN/generate_files/prj/input_files/";
	# endif
	# ifdef SIM
	std::string path = "C:/Users/zju/Desktop/TrackGNN/generate_files/prj/input_files/";
	# endif

	# ifdef EMU
	for (int graph_idx=0; graph_idx<host_num_of_graphs; graph_idx++)
	# endif
	# ifdef SIM
	for (int graph_idx=0; graph_idx<tb_num_of_graphs; graph_idx++)
	# endif	
	{
		// File names.
		std::string graph_idx_string = std::to_string(graph_idx);
		std::string node_feature_file = path + "node_feature_" + graph_idx_string + ".txt";
		std::string adj_list_file = path + "adj_list_" + graph_idx_string + ".txt";
		std::string result_file = path + "result_" + graph_idx_string + ".txt";

		// Temporary data.
		std::ifstream file_in;
		IDX_TYPE num;
		int uni_idx;
		int dim_idx;
		DT0_TYPE node_feature_cache;
		IDX_TYPE adj_list_cache;
		DT1_TYPE result_cache;

		# ifdef EMU
		// Get the number of nodes.
		file_in.open(node_feature_file);
		file_in >> num;
		host_num_of_nodes[graph_idx] = ceildiv(num, (IDX_TYPE)NODE_FEA_DIM);
		file_in.close();
		// Get the number of edges.
		file_in.open(adj_list_file);
		file_in >> num;
		host_num_of_edges[graph_idx] = ceildiv(num, (IDX_TYPE)2);
		file_in.close();
		// Get node features.
		file_in.open(node_feature_file);
		file_in >> num;
		for (int i=0; i<num; i++) {
			file_in >> node_feature_cache;
			uni_idx = graph_idx*INPUT_FEA_COUNT + i/NODE_FEA_DIM/INPUT_FEA_UNI_SIZE;
			dim_idx = (i/NODE_FEA_DIM)%INPUT_FEA_UNI_SIZE*INPUT_FEA_DIM_SIZE + i%NODE_FEA_DIM;
			host_node_feature[uni_idx][dim_idx] = node_feature_cache; }
		file_in.close();
		// Get the adjacent list.
		file_in.open(adj_list_file);
		file_in >> num;
		for (int i=0; i<num; i++) {
			file_in >> adj_list_cache;
			uni_idx = graph_idx*INPUT_ADJ_COUNT + i/2/INPUT_ADJ_UNI_SIZE;
			dim_idx = (i/2)%INPUT_ADJ_UNI_SIZE*INPUT_ADJ_DIM_SIZE + i%2;
			host_adj_list[uni_idx][dim_idx] = adj_list_cache; }
		file_in.close();
		// Get scalar results.
		file_in.open(result_file);
		file_in >> num;
		for (int i=0; i<num; i++) {
			file_in >> result_cache;
			int idx_0 = graph_idx*OUTPUT_RESULT_COUNT + i/EDGE_UNI_PARA;
			int idx_1 = i%EDGE_UNI_PARA;
			host_result_read[idx_0][idx_1] = result_cache; }
		file_in.close();
		# endif

		# ifdef SIM
		// Get the number of nodes.
		file_in.open(node_feature_file);
		file_in >> num;
		tb_num_of_nodes[graph_idx] = ceildiv(num, (IDX_TYPE)NODE_FEA_DIM);
		file_in.close();
		// Get the number of edges.
		file_in.open(adj_list_file);
		file_in >> num;
		tb_num_of_edges[graph_idx] = ceildiv(num, (IDX_TYPE)2);
		file_in.close();
		// Get node features.
		file_in.open(node_feature_file);
		file_in >> num;
		for (int i=0; i<tb_num_of_nodes[graph_idx]; i++)
		for (int j=0; j<NODE_FEA_DIM; j++) {
			file_in >> node_feature_cache;
			uni_idx = (i / INPUT_FEA_UNI_SIZE);
			dim_idx = (i % INPUT_FEA_UNI_SIZE) * INPUT_FEA_DIM_SIZE + j;
			tb_node_feature[graph_idx][uni_idx][dim_idx] = node_feature_cache; }
		file_in.close();
		// Get the adjacent list.
		file_in.open(adj_list_file);
		file_in >> num;
		for (int i=0; i<tb_num_of_edges[graph_idx]; i++)
		for (int j=0; j<2; j++) {
			file_in >> adj_list_cache;
			uni_idx = (i / INPUT_ADJ_UNI_SIZE);
			dim_idx = (i % INPUT_ADJ_UNI_SIZE) * INPUT_ADJ_DIM_SIZE + j;
			tb_adj_list[graph_idx][uni_idx][dim_idx] = adj_list_cache; }
		file_in.close();
		// Get scalar results.
		file_in.open(result_file);
		file_in >> num;
		for (int i=0; i<tb_num_of_edges[graph_idx]; i++) {
			file_in >> result_cache;
			tb_result_read[graph_idx][i/EDGE_UNI_PARA][i%EDGE_UNI_PARA] = result_cache; }
		file_in.close();
		# endif
	}
}

