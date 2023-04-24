# SVM classification analysis for fMRI data
This MATLAB code aims to identify brain regions where different uncertainty conditions can be distinguished. It uses a n-fold partitioning scheme and a Support Vector Machine (SVM) classifier to classify the data.

## Requirements
CosmoMVPA toolbox (https://github.com/CoSMoMVPA/CoSMoMVPA)

libSVM for MATLAB (https://www.csie.ntu.edu.tw/~cjlin/libsvm/)

## Usage
1. Add paths to necessary toolboxes by modifying the following lines:

``` matlab
addpath(genpath('/path/to/CoSMoMVPA-master'));
addpath(genpath('/path/to/libsvm-3.25'));
```

2. Define a list of participants by modifying the following line:
``` matlab
subjects= [1 2 3];
```

3. Define the directory containing the data files, and the ROIs to be analyzed by modifying the following lines:
``` matlab
main_dir=('/path/to/data');
```

4. Save the results by modifying the following line to specify the path to the directory where the results should be saved:
``` matlab
name_to_save=sprintf('/path/to/data/MVPA_brain_location/%s',folder_name);
```

## Steps
1. For each subject and region of interest (ROI):
  - Load data files for each run of the experiment
  - Create datasets for each run
  - Stack all datasets together
  - Set sphere for searchlight
  - Define measure and its arguments; here cross-validation with SVM classifier to compute classification accuracies
  - Run searchlight
  - Save results
  - Create struct with all subjects
  

2. For all subjects and regions of interest:
  - Stack datasets together
  - Average samples
  - Save average results.
  
## Outputs
The code outputs a classification accuracy map for each ROI, which can be used to identify brain regions where uncertainty conditions can be distinguished.
