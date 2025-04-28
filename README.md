# BodyNet2Brain
![Conceptual Framework](banner.png)
Mapping representations in the visual cortex with body-parsing neural networks

The purpose of this study is to unravel representations in visual cortical regions with body-parsing neural networks (i.e., pose estimation and body segmentation models). Specifically, it aims to:

1. Perform Representational Similarity Analysis (RSA) between neural networks and visual cortical regions.
2. Conduct variance partitioning to determine the contributions of each neural network model to explaining the variance in the neural responses in visual cortical regions.

# Analysis Results
To replicate our analysis results, you may simply run the corresponding notebooks under 'analyses'. As the notebook assumes you have certain files (e.g., preprocessed RDMs, image_indices), you can create a shortcut to our ['ML Project' folder assets](https://drive.google.com/drive/folders/1uMF26zAsWn0lrr61fyStb5PX_yv_S96f?usp=drive_link) containing thesee files in your Google Drive. If you'd rather get the RDMs from scratch, please follow the preprocessing pipeline described in the next section.

Our analysis notebooks include both RSA and variance partitioning, and running them should yield you the same values we arrived at as we fixed the random seed for the train/validation splits, which was the only source of variation in the main analysis.
 
# Preprocessing Pipeline:

Trial list extraction:
Run findBodyImagesForSubjects2 for the excel that contains information about all subjects, trial numbers and images shown, alongside information about all images that included bodies (single_person_image_coco_ids.txt). the results will be all trialIDs for which single body images were shown to a subject (usually each image was shown three times).

ROI segmentation:
After downloading individually defined ROIs for each subject, segmenting using segment_nifti for body selective ROIs (EBA, FBA_1, FBA_2, MTL), segment_nifti_visual_ml for visual ROIs (V1v, V1d, V2v, V2d, V3v, V3d, hV4), or segment_nifti_face_ml for face selective ROIs (OFA, FFA_1, FFA_2).

Beta extraction:
Once extracting both ROIs and trials in which body images were shown for every subject, run batchProcessBetas to create excels of activations in each voxel of a given ROI for a chosen subject. 
Then, run formatExcel to transpose the data and dot multiply to correct for real beta values. 
Finally, batchExcelRDMs can be run to create RDMs in spreadsheets and visual graphs from the formatted excel spreadsheets.
