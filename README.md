# BodyNet2Brain
![Conceptual Framework](banner.png)
Mapping representations in the visual cortex with body-parsing neural networks

The purpose of this study is to unravel representations in visual cortical regions with body-parsing neural networks (i.e., pose estimation and body segmentation models). Specifically, it aims to:

1. Perform Representational Similarity Analysis (RSA) between neural networks and visual cortical regions.
2. Conduct variance partitioning to determine the contributions of each neural network model to explaining the variance in the neural responses in visual cortical regions.

*draft*

Unraveling Representations in Body-selective Brain Regions with Body-aware Deep Neural Networks


BodyNet2Brain
This repository explores how the human visual cortex represents body information and compares it to representations from body-aware deep neural networks. Using Representational Similarity Analysis (RSA), the project investigates how different visual regions map onto body parsing tasks performed by computational models such as YOLO-based segmentation and pose networks.

Objectives & Hypothesis:
To link brain activity (from body-selective regions) with deep learning features extracted from body parsing models.

To isolate contributions of part-based and holistic representations in the human visual system.

Hypothesis: Visual areas show distinct representational alignments with segmentation and pose features, reflecting compositional processing of body structure.

Assumptions:
The fMRI data is filtered from the NSD (Natural Scenes Dataset) for ~4450 images containing single human bodies.

Model features are derived from YOLOv8, fine-tuned for body segmentation and pose estimation.

Neural activations and model features are compared using RSA and variance partitioning.

Project Structure:
bash
Copy
Edit
bodynet-brainRSA/
├── assets/                    # Visualization figures (RDMs, bar plots, venn, Grad-CAM)
├── data/                      # fMRI data and preprocessed model activations
├── src/
│   ├── models/                # YOLO-based models and wrappers
│   ├── rsa/                   # RSA and variance partitioning scripts
│   ├── preprocess/            # Scripts to extract and clean brain and model data
│   ├── analysis/              # Statistical analysis and plotting
├── results/                   # Output figures and tables
├── notebooks/                # Jupyter notebooks for exploration
├── README.md                 # This file
├── requirements.txt          # Core Python dependencies
├── LICENSE
Key Stages of the Project
Data Acquisition

Extract ~4450 body-centric images from NSD dataset.

Corresponding fMRI activation patterns from visual cortex ROIs.

Model Feature Extraction

Use YOLOv8-based models for pose and segmentation tasks.

Extract layer-wise features for each image.

Representational Analysis

Run RSA to compare model features with neural RDMs.

Use variance partitioning to disentangle shared and unique variance.

Evaluation & Visualization

Create bar plots of RSA values per region and model.

Plot variance Venn diagrams to show model overlaps.

Visualize searchlight RSA on cortical surface (glass brain).

Optional: regress out low-level features to isolate higher-level abstraction.

Sample Output:

RSA Correlation (bar plots)

Unique vs. shared variance (Venn diagrams)

Glass brain maps

Grad-CAM activations (optional)

Important Definitions & Key Parameters
RSA: Representational Similarity Analysis (compares structure of activation patterns).

Variance Partitioning: Disentangles what variance is uniquely explained by each model.

ROIs: V1, EBA, FBA, OPA, FFA, etc.

YOLOv8: State-of-the-art object detection model, fine-tuned for body understanding.

How To Use:
Clone the Repository

bash
Copy
Edit
git clone https://github.com/yourusername/bodynet-brainRSA.git  
cd bodynet-brainRSA  
Set Up Environment

bash
Copy
Edit
pip install virtualenv  
python3 -m venv venv  
source venv/bin/activate  
python -m pip install --upgrade pip  
Install Dependencies

bash
Copy
Edit
pip install -r requirements.txt  
Run RSA Pipeline

bash
Copy
Edit
python -m src.rsa.run_rsa --roi EBA --model segmentation  
Generate Visualizations

bash
Copy
Edit
python -m src.analysis.plot_rsa_results  
python -m src.analysis.plot_variance_partitioning  
Unit Tests:
Make sure to install dev tools:

bash
Copy
Edit
pip install pytest
Run:

bash
Copy
Edit
pytest tests/
References:
Allen et al., Nature Neuroscience (2022) – Natural Scenes Dataset

Zhu et al., PNAS (2024) – Body pose representation

Kriegeskorte et al., Frontiers in Systems Neuro (2008) – RSA

Downing et al., Science (2001) – Body-selective cortical areas

Redmon et al., CVPR (2016) – YOLO model

Dwivedi et al., JOCN (2021) – Scene parsing + brain mapping

