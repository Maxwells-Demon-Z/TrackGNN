
# TrackGNN: A Highly Parallelized and FIFO-Balanced Accelerator for the Particle-Track-Classification GNN on FPGAs

The aim of this project is to develop an accelerator specifically tailored for processing data obtained from colliders, focusing on verifying the authenticity of tracks. The accelerator's design is rooted in and adapted from the principles outlined in the FlowGNN article at https://arxiv.org/abs/2204.13103. Our algorithm for Graph Neural Networks (GNN) draws inspiration from the sPhenix article and references open-source code available at https://arxiv.org/abs/1810.06111 and https://bitbucket.org/dtyu/trigger-detection-pipeline.  

Should you encounter any issues or have questions, please feel free to contact Hanqing Zhang via email at hanqing.zhang@zju.edu.cn.  

Using this project involves modifying parameters within the file located at './TrackGNN/main/generate.py' and executing it. This file will extract and generate necessary files from the folder './TrackGNN/src_files/'. All essential files required for building the accelerator will automatically populate in the './TrackGNN/prj/' folder shortly. Subsequently, the Vitis tool can be utilized for synthesis and implementation without requiring additional complex settings.  

The article was co-authored by Hanqing Zhang and Buqing Xu under the guidance of Dr. Cong Callie Hao. It remains unpublished and is currently under review by the International Conference on Field Programmable Logic and Applications (FPL).  

