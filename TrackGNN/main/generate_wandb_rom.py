
# Import packages.
import torch
import os

# Write the data into a file.
def write_float_data(
        file_path,
        file_name,
        file_data):
    
    # Check the path.
    if not os.path.exists(file_path):
        os.makedirs(file_path)
    
    # Write the data into the file.
    file_name = open(file_path + file_name, 'w')
    for d in file_data[: -1]:
        file_name.write(str(float(d)) + ', ')
    file_name.write(str(float(file_data[-1])))

# Main function.
def generate_wandb_rom():
    
    # Print the information.
    print('INFO: Start generating weights and biases ROMs.')

    # Initialize the max data.
    max_data_of_wandb = 0.0

    # Get the data from the file.
    src_file_path = './src_files/wandb_src_files/'
    src_file_name = os.listdir(src_file_path)[0]
    file = torch.load(src_file_path + src_file_name, map_location='cuda:0')
    wandb_data = file['model']

    # Get dimensions.
    shape = list(wandb_data['input_network.0.weight'].shape)
    node_emb_dim = shape[0]
    node_fea_dim = shape[1]

    # Get the destination path.
    dst_file_path = './prj/data_files/'

    # Extract data of the input network.
    data = wandb_data['input_network.0.weight'].flatten()
    max_data_of_wandb = max(float(data.abs().max()), max_data_of_wandb)
    write_float_data(dst_file_path, 'input_nw_weight.txt', data)
    data = wandb_data['input_network.0.bias'].flatten()
    max_data_of_wandb = max(float(data.abs().max()), max_data_of_wandb)
    write_float_data(dst_file_path, 'input_nw_bias.txt', data)

    # Extract data of the node network.
    for i in range(4):
        if (i == 0):
            data = wandb_data['node_network.network.' + str(2*i) + '.weight']
            max_data_of_wandb = max(float(data.abs().max()), max_data_of_wandb)
            data_0 = data[:, 0*node_emb_dim: 1*node_emb_dim].flatten()
            data_1 = data[:, 1*node_emb_dim: 2*node_emb_dim].flatten()
            data_2 = data[:, 2*node_emb_dim: 3*node_emb_dim].flatten()
            write_float_data(dst_file_path, 'node_nw_weight_0_0.txt', data_0)
            write_float_data(dst_file_path, 'node_nw_weight_0_1.txt', data_1)
            write_float_data(dst_file_path, 'node_nw_weight_0_2.txt', data_2)
        else:
            data = wandb_data['node_network.network.' + str(2*i) + '.weight'].flatten()
            max_data_of_wandb = max(float(data.abs().max()), max_data_of_wandb)
            write_float_data(dst_file_path, 'node_nw_weight_' + str(i) + '.txt', data)
        data = wandb_data['node_network.network.' + str(2*i) + '.bias'].flatten()
        max_data_of_wandb = max(float(data.abs().max()), max_data_of_wandb)
        write_float_data(dst_file_path, 'node_nw_bias_' + str(i) + '.txt', data)

    # Extract data of the edge network.
    for i in range(4):
        if (i == 0):
            data = wandb_data['edge_network.network.' + str(2*i) + '.weight']
            max_data_of_wandb = max(float(data.abs().max()), max_data_of_wandb)
            data_0 = data[:, 0*node_emb_dim: 1*node_emb_dim].flatten()
            data_1 = data[:, 1*node_emb_dim: 2*node_emb_dim].flatten()
            write_float_data(dst_file_path, 'edge_nw_weight_0_0.txt', data_0)
            write_float_data(dst_file_path, 'edge_nw_weight_0_1.txt', data_1)
        else:
            data = wandb_data['edge_network.network.' + str(2*i) + '.weight'].flatten()
            max_data_of_wandb = max(float(data.abs().max()), max_data_of_wandb)
            write_float_data(dst_file_path, 'edge_nw_weight_' + str(i) + '.txt', data)
        data = wandb_data['edge_network.network.' + str(2*i) + '.bias'].flatten()
        max_data_of_wandb = max(float(data.abs().max()), max_data_of_wandb)
        write_float_data(dst_file_path, 'edge_nw_bias_' + str(i) + '.txt', data)

    # Print the information.
    print('INFO: The feature dimension is %d.' % (node_fea_dim))
    print('INFO: The embedding dimension is %d.' % (node_emb_dim))
    print('INFO: The max data of weights and biases is %.4f.' % (max_data_of_wandb))
    print('INFO: Finish generating weights and biases ROMs.')

    # Return the result.
    return node_fea_dim, node_emb_dim, max_data_of_wandb

