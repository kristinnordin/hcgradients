%% Setup Paths
% These paths point to the directories where required data and surfaces are stored.
surf_path = '.\Congrads\SurfaceProjections\Dependencies\Surfaces\';
hcp_surf_path = '.\HCP_data\HCP_S1200_GroupAvg_v1';
margulies_path = '.\Margulies_etal_2016_gradients\Gradients_Margulies2016\';
atlas_path = '.\ATLASES\Schaefer\HCP\fslr32k\cifti\';

% Path to surface-projected pmap data from Workbench
% For Matlab integration: https://github.com/rudyvdbrink/Surface_projection
pmap_path = '.\Congrads\SurfaceProjections\Dynamic_v8\pmap_data\';

%% Load Surfaces
% Load the left and right hemisphere midthickness surfaces. These are standard geometry files in GIFTI format.
s_lh = convert_surface(fullfile(surf_path, 'S1200.L.midthickness_MSMAll.32k_fs_LR.surf.gii'));
s_rh = convert_surface(fullfile(surf_path, 'S1200.R.midthickness_MSMAll.32k_fs_LR.surf.gii'));

%% Load Medial Wall vertices to mask
% The medial wall mask indicates vertices that are not part of the cortical ribbon.
Medial_wall = cifti_read(fullfile(hcp_surf_path,'Human.MedialWall_Conte69.32k_fs_LR.dlabel.nii'));
Medial_wall_lh = cifti_struct_dense_extract_surface_data(Medial_wall, 'CORTEX_LEFT') == 0;
Medial_wall_rh = cifti_struct_dense_extract_surface_data(Medial_wall, 'CORTEX_RIGHT') == 0; 

%% Load Margulies Gradients
% Load precomputed functional gradients from Margulies et al. (2016).
Margulies_grads = cifti_read(fullfile(margulies_path,'hcp.gradients.dscalar.nii'));
Margulies_left = cifti_struct_dense_extract_surface_data(Margulies_grads, 'CORTEX_LEFT');
Margulies_right = cifti_struct_dense_extract_surface_data(Margulies_grads, 'CORTEX_RIGHT');

% Set medial wall vertices to NaN, as they are not cortical.
Margulies_left(Medial_wall_lh==0,:) = nan;
Margulies_right(Medial_wall_rh==0,:) = nan;

%% Load Hippocampus pmaps
% We have pmaps (cortical gradient maps) from hippocampus projections.
% We loop through each gradient (g = 1:3) and hemisphere (L,R) and load and mask them.
files = dir(fullfile(pmap_path,'*.func.gii'));
hemis = {'L','R'};
pmap = struct();

for g = 1:3
    for h = 1:2
        % Load and hippocampal pmaps for left and right surfaces
        pmap_lh_gii = gifti(fullfile(files(1).folder, sprintf('lh.hippocampus.pmaps000%i_smoothed_%s.func.gii',g-1,hemis{h})));
        pmap_rh_gii = gifti(fullfile(files(1).folder, sprintf('rh.hippocampus.pmaps000%i_smoothed_%s.func.gii',g-1,hemis{h})));
        
        % Store in a structured field naming convention for later retrieval
        pmap.(sprintf('%s_HIPP_G%i_LH', hemis{h}, g)) = pmap_lh_gii.cdata;
        pmap.(sprintf('%s_HIPP_G%i_RH', hemis{h}, g)) = pmap_rh_gii.cdata;
    end
end

%% Parcellating the data using Schaefer1000 parcels
ParcLUT = cifti_read(fullfile(atlas_path,'Schaefer2018_1000Parcels_7Networks_order.dlabel.nii'));

% Extract network names from the lookup table (LUT)
NetworkNames = string({ParcLUT.diminfo{1,2}.maps.table.name}');
% Remove the first entry (background) then extract network parts
NetworkNames = extractBetween(NetworkNames(2:end),14,'_');
[uniqueNetworkNames,~,NetworkNumbers] = unique(NetworkNames,'stable');

% Extract the parcel keys (unique identifiers) from LUT
NetworkKey = cell2mat({ParcLUT.diminfo{1,2}.maps.table.key}');
NetworkKey = NetworkKey(2:end); % remove the background key

% Create a label vector for all vertices (concatenate LH and RH)
labels = zeros(32492*2,1);
for i = 1:numel(NetworkKey)
    indx = ParcLUT.cdata == NetworkKey(i);
    labels(indx) = NetworkNumbers(i);
end

% Separate LH and RH labels if needed
labels_lh = labels(1:32492);
labels_rh = labels(32493:end);
parcs = ParcLUT.cdata;
parcs_no_zero = parcs(parcs~=0);

%% Spin Permutation
% Spin permutations to generate null distributions by rotating the data on a spherical representation.
n_perm = 1000; 
batch_size = 100; 

% Load spherical surfaces for performing spin permutations
% Using https://brainspace.readthedocs.io toolbox
lh_sphere = convert_surface(fullfile(hcp_surf_path,'S1200.L.sphere.32k_fs_LR.surf.gii'));
rh_sphere = convert_surface(fullfile(hcp_surf_path,'S1200.R.sphere.32k_fs_LR.surf.gii'));

% Combine spheres into a single structure
sphere = struct();
sphere.tri = [lh_sphere.tri; rh_sphere.tri];
sphere.coord = [lh_sphere.coord, rh_sphere.coord];

% Extract first 5 gradients from Margulies data for permutation
data_for_spins = [Margulies_left(:,1:5); Margulies_right(:,1:5)];

% Initialize cell array to store parcel-level spin results
Y_rand_parc = cell(1,5);

% Perform spin permutations in batches
for i = 1:(n_perm/batch_size)
    % spin_permutations is a function from BrainSpace library
    Y_rand = spin_permutations(data_for_spins, sphere, batch_size);

    for k = 1:5
        % For each gradient, compute the mean value within each parcel and store it
        Y_rand_parc{k} = [Y_rand_parc{k}, grpstats(squeeze(Y_rand{1}(parcs~=0,k,:)), parcs_no_zero, {'mean'})];
    end
end

% Save the permutation results for later use
save('Margulies_5grads_1000spins_Schaefer1000parcs.mat', 'Y_rand_parc');

%% Test Correlation between Hipp pmaps and Margulies Gradients
% Compute observed correlation and p-value by comparing to spin permutation null distribution.

% Preallocate a table to store results
spin_results = table('Size',[6*5,4],'VariableTypes',{'string', 'string','double','double'},...
    'VariableNames',{'Hipp_pmap', 'Margulies_grad','Spearmans_corr', 'P_spin'});

step = 1;
% Loop through each Margulies gradient, Hipp gradient, and hemisphere combination
for MargG = 1:5
    for pmapG = 1:3
        for hh = ['L','R']
            % Identify the correct pmaps from the structure (LH+RH combined)
            pmap_names = {sprintf('%s_HIPP_G%i_LH',hh,pmapG), sprintf('%s_HIPP_G%i_RH',hh,pmapG)};
            pmap_table = struct2table(pmap);
            colind = find(ismember(pmap_table.Properties.VariableNames, pmap_names));
            
            % Combine LH and RH data and remove parcels not in the atlas (parcs~=0)
            pmap_tmp = [pmap_table{:,colind(1)}; pmap_table{:,colind(2)}];
            pmap_tmp = pmap_tmp(parcs~=0);
            % Compute average pmap value per parcel
            pmap_tmp_parcs = grpstats(pmap_tmp,parcs_no_zero,{'mean'});

            % Similarly, extract the current Margulies gradient data and average per parcel
            Marg_tmp = [Margulies_left(:,MargG); Margulies_right(:,MargG)];
            Marg_tmp = Marg_tmp(parcs~=0);
            Marg_tmp_parcs = grpstats(Marg_tmp,parcs_no_zero,{'mean'});

            % Compute observed Spearman correlation
            obs_corr = corr(Marg_tmp_parcs, pmap_tmp_parcs,'rows','complete','type','Spearman');
            % Compute null correlation distribution from spin permutations
            null_corr = corr(Y_rand_parc{MargG}, pmap_tmp_parcs,'rows','complete','type','Spearman');
            
            % P-value is the proportion of null correlations greater than observed
            p_spin = mean(obs_corr > null_corr);
            if p_spin >= 0.975; p_spin = 1-p_spin; end

            % Store results
            spin_results(step,:) = {sprintf('%s_HIPP_G%i',hh,pmapG), sprintf('G%i',MargG), obs_corr, p_spin};
            step = step+1;
        end
    end
end

%% Plot Associations
% Plot the linear associations between Hipp pmaps (per parcel) and Margulies gradients (per parcel).

figure; % Create a figure to hold subplots
plot_step = 1;

for Pmap_Gnum = 1:3
    for HippHemi = ['L','R']
        for M_Gnum = 1:3
            
            % Extract Hipp pmaps from struct
            pmap_names = {sprintf('%s_HIPP_G%i_LH',HippHemi,Pmap_Gnum),...
                          sprintf('%s_HIPP_G%i_RH',HippHemi,Pmap_Gnum)};
            pmap_table = struct2table(pmap);
            colind = find(ismember(pmap_table.Properties.VariableNames,pmap_names));
            
            % Combine LH+RH and select only parcel vertices
            pmap_tmp = [pmap_table{:,colind(1)}; pmap_table{:,colind(2)}];
            pmap_tmp = pmap_tmp(parcs~=0);
            pmap_tmp_parcs = grpstats(pmap_tmp,parcs_no_zero,{'mean'});

            % Extract the selected Margulies gradient and average by parcel
            Marg_tmp = [Margulies_left(:,M_Gnum); Margulies_right(:,M_Gnum)];
            Marg_tmp = Marg_tmp(parcs~=0);
            Marg_tmp_parcs = grpstats(Marg_tmp, parcs_no_zero, {'mean'});
            
            % Fit a linear model to get slope, intercept, and confidence intervals
            mdl = fitlm(Marg_tmp_parcs, pmap_tmp_parcs);
            % Predict fitted values and confidence intervals at the given X-values
            [Ypred, Ypred_CI] = predict(mdl, Marg_tmp_parcs);
            
            % X, Y for plotting the data
            X = Marg_tmp_parcs;
            Y = pmap_tmp_parcs;
            
            % Confidence interval around the fitted line (upper and lower)
            Y_CI_upper = Ypred_CI(:,2) - Ypred; 
            Y_CI_lower = Ypred - Ypred_CI(:,1);
            
            % Create a subplot for each combination
            subplot(3,6,plot_step);
            % Plot the mean fit line with confidence intervals using confplot
            confplot(X, Ypred, Y_CI_lower, Y_CI_upper, 'Color',[1 0 0],'LineWidth',2)

            % Add correlation to the title
            current_corr = corr(Marg_tmp_parcs,pmap_tmp_parcs,'type','Spearman','rows','complete');
            title(sprintf('rho=%.2f',current_corr))
            xlabel(sprintf('Margulies G%i',M_Gnum));
            ylabel(sprintf('%s-Hipp G%i',HippHemi,Pmap_Gnum));
            plot_step = plot_step+1;
        end
    end
end
