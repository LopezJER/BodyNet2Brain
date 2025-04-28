# BodyNet2Brain
![Conceptual Framework](banner.png)
Mapping representations in the visual cortex with body-parsing neural networks
https://drive.google.com/drive/folders/1uMF26zAsWn0lrr61fyStb5PX_yv_S96f?usp=drive_link

The purpose of this study is to unravel representations in visual cortical regions with body-parsing neural networks (i.e., pose estimation and body segmentation models). Specifically, it aims to:

1. Perform Representational Similarity Analysis (RSA) between neural networks and visual cortical regions.
2. Conduct variance partitioning to determine the contributions of each neural network model to explaining the variance in the neural responses in visual cortical regions.

Hypotheses? Assumptions?


 # Project Structure: 

*create at the end*


Key stages of the project #

1. Extract lists of which out of ~4450 body-centric images from NSD dataset were shown to each subject.

Example?


2. Find corresponding fMRI activation patterns from visual cortex ROIs.

Example?


Model Feature Extraction3. 

Use YOLOv8-based models for pose and segmentation tasks.
or sapiens, mention other tasks. 

4. Extract layer-wise features for each image.

5. Representational Analysis

Run RSA to compare model features with neural RDMs.

Use variance partitioning to disentangle shared and unique variance6. .

Evaluation & Visualization7. 

Create bar plots of RSA values per region and model.

Plot variance Venn diagrams to show model overlaps.

Visualize searchlight RSA on cortical surface (glass brain).


Sample outputs:
Various graphs we got as results


# How To Use:
idk


# Code pipeline:

Trial list extraction:
Run findBodyImagesForSubjects2 for the excel that contains information about all subjects, trial numbers and images shown, alongside information about all images that included bodies (single_person_image_coco_ids.txt). the results will be all trialIDs for which single body images were shown to a subject (usually each image was shown three times).

ROI segmentation:
After downloading individually defined ROIs for each subject, segmenting using segment_nifti for body selective ROIs (EBA, FBA_1, FBA_2, MTL), segment_nifti_visual_ml for visual ROIs (V1v, V1d, V2v, V2d, V3v, V3d, hV4), or segment_nifti_face_ml for face selective ROIs (OFA, FFA_1, FFA_2).

Beta extraction:
Once extracting both ROIs and trials in which body images were shown for every subject, run batchProcessBetas to create excels of activations in each voxel of a given ROI for a chosen subject. 
Then, run formatExcel to transpose the data and dot multiply to correct for real beta values. 
Finally, batchExcelRDMs can be run to create RDMs in spreadsheets and visual graphs from the formatted excel spreadsheets.


RSA:


Variance Partitioning:

After splitting vectorized RDMs (vectorize_rdm.m) to train and val, for the neural data as well as all four deep neural models (Pose estimation, Body segmentation, Surface normal estimation, Depth estimation) for each subject, run sapiens_allmodels_batch_variance. The output will be an excel spreadsheet containing the values of metrics such as the full model R^2, unique variances of each model, shared variances between each combination of models, etc. 
The preliminary data can be visualized as a scatter plot for each ROI using sapiens_allmodels_graph. 
The results can also be averaged across subjects using updated_average_lateralized_var_stats and plotted using interactive_group_graphs. 
Permutation of the data for obtaining significance levels is run using group_permutation_excel. The results can be visualized with group_permutation_graphs.


Plotting:





 # References: 
Allen et al., Nature Neuroscience (2022) – Natural Scenes Dataset

Zhu et al., PNAS (2024) – Body pose representation

Kriegeskorte et al., Frontiers in Systems Neuro (2008) – RSA

Downing et al., Science (2001) – Body-selective cortical areas

Redmon et al., CVPR (2016) – YOLO model

Dwivedi et al., JOCN (2021) – Scene parsing + brain mapping

