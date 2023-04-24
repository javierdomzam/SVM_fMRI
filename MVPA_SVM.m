%%% This analysis identified brain regions where uncertainty conditions can be
%%% distinguished using a n-fold partitioning scheme and a Support Vector Machine (SVM) classifier

%%% needs libSVM for matlab to work (https://www.csie.ntu.edu.tw/~cjlin/libsvm/).
%%% and CosmoMVPA toolbox

clear all % clear all variables from memory
clc % clear command window

% add paths to necessary toolboxes
addpath(genpath('/path/to/CoSMoMVPA-master'));
addpath(genpath('/path/to/libsvm-3.25'));

% Define participants list
subjects= [1 2 3];

% Define data filenames & load data from even and odd runs
main_dir=('/path/to/Attention_Decisions_fMRI');
rois={'Brain'};
nrois=numel(rois);

analysis_1=0;

% Loop through each region of interest and subject to perform the analysis

    for i_roi = 1:nrois
        for ss =1:length(subjects)
            
            %Participant code
            s=subjects(ss);
            if s>=10
                folder_name =sprintf('sub-%d',s);
            else
                folder_name =sprintf('sub-0%d',s);
            end
            disp (folder_name); % display the participant's folder name
            
            % brain mask in standard space (NMI) from third level FSL. It should be
            % the same for all participants (same number of voxels)
            mask_fn=fullfile('/usr/local/fsl/data/standard/MNI152_T1_1mm_brain_mask.nii.gz'); % whole brain mask
            
            % load data files for each run of the experiment
            data_path=fullfile(main_dir,folder_name);
            name_to_save=sprintf('/path/to/data/MVPA_brain_location/%s',folder_name);
            
            if ~exist(name_to_save, 'dir')
                mkdir(name_to_save) % create the directory if it does not already exist
            end
            
            for nc=1:4
                % cope files for each run
                cope_1=sprintf('%s/func/GLM1_run%d_output.feat/reg_standard/stats/cope5.nii.gz',data_path,nc); %cope 3 = HVR
                cope_2=sprintf('%s/func/GLM1_run%d_output.feat/reg_standard/stats/cope6.nii.gz',data_path,nc); %cope 4 = HVL
                
                % set condition labels for each run
                if nc==1 || nc==3
                    c=1;
                elseif nc==2 || nc==4
                    c=2;
                end
                
                % create datasets for each run
                ds_cope_1=cosmo_fmri_dataset(cope_1,'targets',1,'chunks',nc,'mask',mask_fn);
                ds_cope_2=cosmo_fmri_dataset(cope_2,'targets',2,'chunks',nc,'mask',mask_fn);
                
                if nc==1
                    ds_cope_run1=cosmo_stack({ds_cope_1, ds_cope_2});
                elseif nc==2
                    ds_cope_run2=cosmo_stack({ds_cope_1, ds_cope_2});
                elseif nc==3
                    ds_cope_run3=cosmo_stack({ds_cope_1, ds_cope_2});
                elseif nc==4
                    ds_cope_run4=cosmo_stack({ds_cope_1, ds_cope_2});
                end
            end
            
            % Stack all datasets together
            ds_all=cosmo_stack({ds_cope_run1, ds_cope_run2, ds_cope_run3, ds_cope_run4});
            
            % Set sphere for searchlight
            radius=3; % radius=3 is typical for fMRI datasets
            nbrhood=cosmo_spherical_neighborhood(ds_all,'radius',radius);
            
            % Define measure and its arguments; here crossvalidation with SVM
            % classifier to compute classification accuracies
            args=struct();
            args.classifier = @cosmo_classify_svm;
            args.partitions = cosmo_nfold_partitioner(ds_all); % Set partition scheme
            measure=@cosmo_crossvalidation_measure;
            
            % Run searchlight
            mvpa_result=cosmo_searchlight(ds_all,nbrhood,measure,args);
            
            % Save results
            file_name_save=sprintf('%s/%s_classificationAccuracy.nii',name_to_save,rois{i_roi});
            cosmo_map2fmri(mvpa_result,file_name_save);
            
            % Create struct with all subjects
            mvpa_first_level{ss}=cosmo_fmri_dataset(mvpa_result,'targets',1,'chunks',s,'mask',mask_fn);
            mvpa_first_level_2{ss}=cosmo_fmri_dataset(mvpa_result,'targets',1,'chunks',1,'mask',mask_fn);
            
            % Clear variables
            clear mvpa_result
            
            % Run across participants level with monte carlo perm
            % Prepare data:
            % - targets=1 for all participants
            % - chunks=1,2.. different value for each participant
            ds_stacked = cosmo_stack(mvpa_first_level);
            ds_stacked_2 = cosmo_stack(mvpa_first_level_2);
            
            % Average samples in ds_stacked_2
            ds_stacked_average=cosmo_average_samples(ds_stacked_2);
            
            % Save average results
            file_name_1=sprintf('/path/to/data/MVPA_brain_location/%s_classification.nii',rois{i_roi});
            cosmo_map2fmri(ds_stacked_average,file_name_1);
            
            % Define clustering neighborhood
            nh=cosmo_cluster_neighborhood(ds_stacked);
            
            % Set options for montecarlo based clustering statistic
            opt=struct();
            opt.cluster_stat='tfce'; % Threshold-Free Cluster Enhancement (Nichols & Smith, 2009, Neuroimage)
            opt.niter=3000%10000; %number of iterations
            opt.h0_mean=0.5; % If the data samples contains classification accuracies, then h0_mean=1/C, with C the number of classes (conditions)
            opt.seed=1; % For replication
            opt.progress=true; % show progress bar
            
            % Apply cluster-based correction
            z_mvpa_second_level=cosmo_montecarlo_cluster_stat(ds_stacked,nh,opt);
            
            % Show z-scores of each feature corrected for multiple comparisons;
            % abs(z_ds.samples)>1.96, survives correction for multiple comparisons
            cosmo_disp(z_mvpa_second_level.samples)
            
            % Save z-scores corrected for multiple comparisons
            file_name_2=sprintf('/path/to/dataMVPA_brain_location/%s_Zscore.nii',rois{i_roi});
            cosmo_map2fmri(z_mvpa_second_level,file_name_2);
            
        end
        
    end
