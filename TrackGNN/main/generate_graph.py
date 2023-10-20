
# Import packages.
import numpy as np
import torch
from torch_scatter import scatter_add
import os
os.environ['KMP_DUPLICATE_LIB_OK']='True'

# Write data into files.
def write_data(
        data_type,
        file_path,
        file_name,
        file_data):
    
    # Open the file.
    file_name = open(file_path + file_name, 'w')

    # Write the length of the data.
    file_name.write(str(len(file_data)) + '\n')

    # Write data.
    for d in file_data:
        if data_type == 'float':
            file_name.write(str(float(d)) + '\n')
        elif data_type == 'integer':
            file_name.write(str(int(d)) + '\n')
        else:
            print('ERROR: Invalid data type.')
            exit()

# Main function.
def generate_graph(
        num_of_graphs,
        num_of_layers):
    
    # Print the information.
    print('INFO: Start generating input files.')

    # Get the data from the file.
    src_file_path = './src_files/wandb_src_files/'
    src_file_name = os.listdir(src_file_path)[0]
    file = torch.load(src_file_path + src_file_name, map_location='cuda:0')
    wandb_data = file['model']

    # Initialize the arrays.
    input_nw_weight = []
    input_nw_bias = []
    node_nw_weight = []
    node_nw_bias = []
    edge_nw_weight = []
    edge_nw_bias = []

    # Extract data.
    input_nw_weight.append(wandb_data['input_network.0.weight'].cpu().numpy().T)
    input_nw_bias.append(wandb_data['input_network.0.bias'].cpu().numpy())
    for i in range(4):
        node_nw_weight.append(wandb_data['node_network.network.' + str(2*i) + '.weight'].cpu().numpy().T)
        node_nw_bias.append(wandb_data['node_network.network.' + str(2*i) + '.bias'].cpu().numpy())
        edge_nw_weight.append(wandb_data['edge_network.network.' + str(2*i) + '.weight'].cpu().numpy().T)
        edge_nw_bias.append(wandb_data['edge_network.network.' + str(2*i) + '.bias'].cpu().numpy())
    
    # Get source data.
    src_file_path = './src_files/graph_src_files/'
    src_file_names = os.listdir(src_file_path)

    # Create destination path.
    dst_file_path = './prj/input_files/'
    if not os.path.exists(dst_file_path):
        os.makedirs(dst_file_path)

    # Initialize the max data.
    node_num_sum = 0
    edge_num_sum = 0
    node_num_max = 0
    edge_num_max = 0
    max_data_of_cpt = 0.0

    # Generate input files.
    for i in range(num_of_graphs):

        # Get a graph.
        src_file_name = src_file_names[i]

        # Get data.
        with np.load(src_file_path + src_file_name, allow_pickle=True) as f:
            # Node feature.
            scaled_hits = f['scaled_hits']
            layer_id = f['layer_id'].reshape(-1, 1)
            node_feature = np.concatenate([scaled_hits, layer_id], axis=1).astype(np.float32)
            # Adjacency list.
            adj_list = f['edge_index'].T
            src, dst = adj_list.T
            # Write data into files.
            node_feature_flatten = node_feature.flatten()
            adj_list_flatten = adj_list.flatten()
            write_data('float', dst_file_path, 'node_feature_' + str(i) + '.txt', node_feature_flatten)
            write_data('integer', dst_file_path, 'adj_list_' + str(i) + '.txt', adj_list_flatten)
            # Update max number of nodes and edges.
            node_num_sum += node_feature.shape[0]
            edge_num_sum += adj_list.shape[0]
            node_num_max = max(node_num_max, node_feature.shape[0])
            edge_num_max = max(edge_num_max, adj_list.shape[0])
            # print(node_feature.shape[0])
            # print(adj_list.shape[0])
            max_data_of_cpt = max(max_data_of_cpt, np.max(node_feature_flatten))
            label = np.logical_and(f['pid'][src] > 0, f['pid'][src] == f['pid'][dst])
        
        # Compute the result.
        x = np.matmul(node_feature, input_nw_weight[0]) + input_nw_bias[0]
        max_data_of_cpt = max(max_data_of_cpt, np.max(x))
        x = np.tanh(x)
        for n in range(num_of_layers + 1):
            x_src, x_dst = x[src], x[dst]
            y = np.concatenate([x_src, x_dst], axis=1)
            cache_0 = np.matmul(y, edge_nw_weight[0]) + edge_nw_bias[0]
            max_data_of_cpt = max(max_data_of_cpt, np.max(cache_0))
            cache_0 = np.tanh(cache_0)
            cache_0 = np.matmul(cache_0, edge_nw_weight[1]) + edge_nw_bias[1]
            max_data_of_cpt = max(max_data_of_cpt, np.max(cache_0))
            cache_0 = np.tanh(cache_0)
            cache_0 = np.matmul(cache_0, edge_nw_weight[2]) + edge_nw_bias[2]
            max_data_of_cpt = max(max_data_of_cpt, np.max(cache_0))
            cache_0 = np.tanh(cache_0)
            cache_0 = np.matmul(cache_0, edge_nw_weight[3]) + edge_nw_bias[3]
            max_data_of_cpt = max(max_data_of_cpt, np.max(cache_0))
            cache_0 = cache_0
            e = 1.0 / (1.0 + np.exp(-cache_0))
            if n == num_of_layers:
                break
            y_src, y_dst = e * x_src, e * x_dst
            mi = scatter_add(torch.from_numpy(y_src), torch.from_numpy(dst), dim=0, dim_size=x.shape[0]).numpy()
            mo = scatter_add(torch.from_numpy(y_dst), torch.from_numpy(src), dim=0, dim_size=x.shape[0]).numpy()
            max_data_of_cpt = max(max_data_of_cpt, np.max(mi))
            max_data_of_cpt = max(max_data_of_cpt, np.max(mo))
            z = np.concatenate([mi, mo, x], axis=1)
            cache_1 = np.matmul(z, node_nw_weight[0]) + node_nw_bias[0]
            max_data_of_cpt = max(max_data_of_cpt, np.max(cache_1))
            cache_1 = np.tanh(cache_1)
            cache_1 = np.matmul(cache_1, node_nw_weight[1]) + node_nw_bias[1]
            max_data_of_cpt = max(max_data_of_cpt, np.max(cache_1))
            cache_1 = np.tanh(cache_1)
            cache_1 = np.matmul(cache_1, node_nw_weight[2]) + node_nw_bias[2]
            max_data_of_cpt = max(max_data_of_cpt, np.max(cache_1))
            cache_1 = np.tanh(cache_1)
            cache_1 = np.matmul(cache_1, node_nw_weight[3]) + node_nw_bias[3]
            max_data_of_cpt = max(max_data_of_cpt, np.max(cache_1))
            cache_1 = np.tanh(cache_1)
            x = x + cache_1

        # Write the result into files.
        result_num = e.T[0]
        write_data('float', dst_file_path, 'result_' + str(i) + '.txt', result_num)

    # Print the information.
    print('INFO: Mean number of nodes: %.4f.' % (node_num_sum / num_of_graphs))
    print('INFO: Mean number of edges: %.4f.' % (edge_num_sum / num_of_graphs))
    print('INFO: Max number of nodes: %d.' % (node_num_max))
    print('INFO: Max number of edges: %d.' % (edge_num_max))
    print('INFO: Max data in processing: %.4f.' % (max_data_of_cpt))
    print('INFO: Finish generating input files.')

    # Return the result.
    return node_num_max, edge_num_max, max_data_of_cpt

