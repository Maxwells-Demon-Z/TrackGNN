
# Import packages.
import numpy as np
import os

# Global original data.
data_origin_sigm = []
data_origin_tanh = []

# Generate original data.
def generate_origin_data(frac):

    # Input accuracy and expansion number.
    accuracy = 2 ** (-frac)
    exp_num = 2 ** frac

    # Generate sigmoid data.
    x = 0.0
    while True:
        y = 1.0 / (1.0 + np.exp(-x))
        data_origin_sigm.append(y)
        x += accuracy
        if y > (1 - accuracy/2):
            for i in range(exp_num):
                y = 1.0 / (1.0 + np.exp(-x))
                data_origin_sigm.append(y)
                x += accuracy
            break

    # Generate tanh data.
    x = 0.0
    while True:
        y = (np.exp(x) - np.exp(-x)) / (np.exp(x) + np.exp(-x))
        data_origin_tanh.append(y)
        x += accuracy
        if y > (1 - accuracy/2):
            for i in range(exp_num):
                y = (np.exp(x) - np.exp(-x)) / (np.exp(x) + np.exp(-x))
                data_origin_tanh.append(y)
                x += accuracy
            break

# Calculate the difference between original data and restored data.
def calculate_diff(func, frac, tail):

    # Input accuracy and sample number.
    accuracy = 2 ** (-frac)
    num_sample = 2 ** tail

    # Select the original data.
    if func == 'SIGM':
        data_origin = data_origin_sigm
    elif func == 'TANH':
        data_origin = data_origin_tanh
    else:
        print('ERROR: Invalid function name!')
        exit(1)

    # Sample the original data.
    data_sample = []
    for i in range(len(data_origin)):
        if (i % num_sample) == 0:
            data_sample.append(data_origin[i])
            if (data_origin[i]) > (1 - accuracy/2):
                break
    
    # Restore the data.
    data_restore = []
    for i in range(len(data_sample)-1):
        y1 = data_sample[i+1]
        y0 = data_sample[i]
        for j in range(num_sample):
            x = j / num_sample
            data_restore.append(y1 * x + y0 * (1 - x))

    # Calculate the difference.
    data_origin = data_origin[: len(data_restore)]
    data_diff_max = np.max(np.abs(np.array(data_origin) - np.array(data_restore)))
    meet_acc = (data_diff_max < accuracy/2)
    
    # Return the result.
    return meet_acc, data_sample

# Calculate the tail and size of the ROM.
def calculate_tail(func, frac):
    tail = 1
    meet_acc_new, data_sample_new = False, []
    while True:
        meet_acc_old, data_sample_old = meet_acc_new, data_sample_new
        meet_acc_new, data_sample_new = calculate_diff(func, frac, tail)
        if meet_acc_new:
            tail += 1
        else:
            tail -= 1
            break
    return tail, len(data_sample_old), data_sample_old

# Main function.
def generate_activation_rom(
        frac,
        node_uni_para,
        edge_uni_para,
        dim_para):
    
    # Print the information.
    print('INFO: Start generating activation function ROMs.')
    
    # Number of MLP layers.
    node_mlp_layer = 4
    edge_mlp_layer = 4

    # Generate original data.
    generate_origin_data(frac)

    # Calculate the tail and size of the ROMs.
    sigm_tail, sigm_size, sigm_sample = calculate_tail('SIGM', frac)
    tanh_tail, tanh_size, tanh_sample = calculate_tail('TANH', frac)

    # Transform the data into string.
    sigm_data = ''
    tanh_data = ''
    for j in range(sigm_size):
        sigm_data += str(sigm_sample[j]) + ', '
    for j in range(tanh_size):
        tanh_data += str(tanh_sample[j]) + ', '

    # Write data into files.
    file_path = './prj/data_files/'
    if not os.path.exists(file_path):
        os.makedirs(file_path)
    with open(file_path + 'activation_sigm_result_rom.txt', 'w') as f:
        for i in range(edge_uni_para):
            f.write(sigm_data)
    with open(file_path + 'activation_sigm_node_nw_rom.txt', 'w') as f:
        for i in range(node_uni_para * node_mlp_layer * dim_para):
            f.write(sigm_data)
    with open(file_path + 'activation_sigm_edge_nw_rom.txt', 'w') as f:
        for i in range(edge_uni_para * edge_mlp_layer * dim_para):
            f.write(sigm_data)
    with open(file_path + 'activation_tanh_node_nw_rom.txt', 'w') as f:
        for i in range(node_uni_para * node_mlp_layer * dim_para):
            f.write(tanh_data)
    with open(file_path + 'activation_tanh_edge_nw_rom.txt', 'w') as f:
        for i in range(edge_uni_para * edge_mlp_layer * dim_para):
            f.write(tanh_data)

    # Print the information.
    print('INFO: Sigmoid tail = %d. ' % (sigm_tail))
    print('INFO: Sigmoid size = %d. ' % (sigm_size))
    print('INFO: Tanh tail = %d. ' % (tanh_tail))
    print('INFO: Tanh size = %d. ' % (tanh_size))
    print('INFO: Finish generating activation function ROMs.')

    # Return the result.
    return sigm_tail, sigm_size, tanh_tail, tanh_size

