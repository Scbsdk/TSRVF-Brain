%% Addpath and read the fMRI4D data using SPM's functions
clear all;close all;

% Creating the output directory if not existing
if (~isfolder('output_dir'))
    mkdir output_dir
end
Fmri4d = spm_read_vols(spm_vol('\preprocessed_fmri_data\bold.nii'));

%% definition of path containing subjects' all tracking results from TractSeg
bundle_dir = 'TractSeg_tck_data';
bundles = dir(bundle_dir); 
bundles = bundles(3:end); 
num_bundles = length(bundles);

%% Calculate shape and function distance matrix
% Read the bundle information by mrtrix3'function and resample
Ds_all = zeros(2000,2000,num_bundles);
Df_all = zeros(2000,2000,num_bundles);
for count=1:num_bundles
    disp(bundles(count).name)
    fibers = read_mrtrix_tracks(fullfile(bundle_dir,bundles(count).name));
    M=fibers.data;
    % Resampled to 25 points per fiber
    N=25;
    for i=1:length(M)
    f=(M{i})';fsample(:,:,i)= ReSampleCurve(f,N);
    end
    % Transform the fiber to SRVF
    for i=1:size(fsample,3)
    q_fiber(:,:,i)=curve_to_qnos(fsample(:,:,i));
    end
    % Calculate the shape distance matrix
    Ds(:,:) = pairdistance(q_fiber,q_fiber);
    Ds_all(:,:,count) = Ds(:,:);
    % Extract the Bold signal from each point along the fiber
    for i=1:size(fsample,3)
        f=fsample(:,:,i);
        for j=1:size(f,2)
            W_all(:,j,i)=extract(f(:,j),Fmri4d);
        end
    end
    % PCA Dimensionality reduction on the Bold signal
    PCA_fmri=reshape(W_all,50000,1200);
    COEFF = pca(PCA_fmri');
    pW_all=reshape(COEFF(:,1:30),30,25,2000);
    clear W_all;clear PCA_fmri;
    % Transform the Bold signal to SRSF then Transform the SRSFs to TSRVF
    for i=1:size(pW_all,3)
        q(:,:,i)=f_to_srsf(pW_all(:,:,i),linspace(0,1,30));
    end
    h=qtr_to_h(q);
    % Calculate the function distance matrix
    Df(:,:) = pairdistance(h,h);
    Df_all(:,:,count) = Df(:,:);
end
    save E:\fibers_cluster\fiber_cluster_data\output_dir\distance.mat Ds_all Df_all 
    
  %% Estimate the number of clusters
  Alpha = 0.5;  % Choose the weight of shape and function
  load dis106016
   % Estimate the number of clusters
  for count=1:num_bundles
      D = Alpha*Ds_all(:,:,count)+(1-Alpha)*Df_all(:,:,count);
      EstimateK(count) = estimate(D);
  end
if (~isfolder('E:\fiber_data\output_dir'))
    mkdir E:\fiber_data\output_dir
end

%% Perform Spectral clustering
  Cluster_IDX = 0;
  for count=1:num_bundles
      fibers = read_mrtrix_tracks(fullfile(bundle_dir,bundles(count).name));
      M=fibers.data;
      D = Alpha*Ds_all(:,:,count)+(1-Alpha)*Df_all(:,:,count);
      Sim = exp(-D.^2);
      idx = spectralcluster(Sim,EstimateK(count),'Distance','precomputed','LaplacianNormalization','symmetric');
      IDX_fibers = tabulate(idx);
      BestK(count) = EstimateK(count)-sum(IDX_fibers(:,3)<5);
     
      % Write the clustering results back to the tck file for visualization
      counter1=1;
      for(i=1:EstimateK(count))
          if(~(IDX_fibers(i,3)<5))
              counter2=1;
              Cluster_fibers{counter1}=cell(1,IDX_fibers(i,2));
              for(j=1:2000)
                  if(idx(j)==i)
                    Cluster_fibers{counter1}{counter2}=M{j};counter2=counter2+1;
                  end
              end
              counter1=counter1+1;
          end
      end
      for(i=1:BestK(count))
          fibers.data=Cluster_fibers{i};
          cd E:\fiber_data\output_dir
          filename = ['Cluster' num2str(Cluster_IDX) '_' bundles(count).name];
          Cluster_IDX = Cluster_IDX+1;
          write_mrtrix_tracks (fibers, filename);
      end
  end
  
   
