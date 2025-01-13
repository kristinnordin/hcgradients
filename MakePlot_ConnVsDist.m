clear all;
fclose all;
close all;

%  extreme points of hippocampus noted after visual inspection
LStartPointmm = [-26 -6 -28];
LEndPointmm = [-26 -40 4];
RStartPointmm = [26 -6 -28];
REndPointmm = [26 -40 4];

Ldirection = LEndPointmm-LStartPointmm; % vector from anterior to posterior on left side
Ldirection = Ldirection./norm(Ldirection); %unity vector on left side
Rdirection = REndPointmm-RStartPointmm;
Rdirection = Rdirection./norm(Rdirection);

intervall = 0:2:46; % vector to divide axis in 2mm segments

stem = '/mnt/Data/Projects/DYNAMIC/congrads/congrads.v8/';

d = dir([stem '40*']);
fpc = 0;
for fp=1:numel(d) % loop over subject directories
    if ~d(fp).isdir
        continue
    end
    fpc=fpc+1;
    strSub{fpc} = d(fp).name;
    for hem=1:2 % first left, then right hemisphere
        if hem==1
            filename = [stem d(fp).name '/lh.hippocampus.norm/lh.hippocampus.cmaps.rearranged2.flipped.nii.gz'];
        else
            filename = [stem d(fp).name '/rh.hippocampus.norm/rh.hippocampus.cmaps.rearranged2.flipped.nii.gz'];
        end
        if ~exist(filename, 'file') % if file do not exist, exit with error message
            disp(sprintf('%s dont exist', filename));
            kdlsgdkslgjkl;
        end
        [head,Y] = dz_ReadNifti(filename); % load gzipped nifti file with inhouse function
        for grad=1:3 % loop over gradients
            
            Ytemp = Y(:,:,:,grad); % select gradient
            [x,y,z] = ndgrid(1:head.dim(1),1:head.dim(2),1:head.dim(3));
            xyz = [x(:) y(:) z(:) ones(numel(x), 1)];
            xyzmm = head.mat(1:3,:)*xyz'; % get voxels in mm
            clear x y z xyz;
            if hem==1
                ToUse = Ytemp(:)>0 & xyzmm(1,:)'<0; %select voxels with positive value and on left side
                SP = LStartPointmm;
                direction = Ldirection;
            else
                ToUse = Ytemp(:)>0 & xyzmm(1,:)'>0; %select voxels with positive value and on right side
                SP = RStartPointmm;
                direction = Rdirection;
            end
            Yr = Ytemp(ToUse); % select voxels with positive value
            xyzmmr = xyzmm(:,ToUse); % select voxels with positive value
            F = NaN(numel(Yr),1);
            for i=1:numel(Yr) % for each voxel
                coord = xyzmmr(:,i);
                P = coord(:)-SP(:); % coordinate relative to StartPoint (vector from anterior to "coord")
                F(i) = dot(direction,P); % vector projected on line from StartPoint to EndPoint
            end
            for j=1:numel(intervall)-1  % for each segment
                if j==1
                    iswithin = F<intervall(j+1); % select voxels within first segment
                elseif j<(numel(intervall)-1)
                    iswithin = F>=intervall(j) & F<intervall(j+1); % select voxels within segment
                else
                    iswithin = F>=intervall(j); % select voxels within last segment
                end
                Ym(fpc,j,hem,grad) = mean(Yr(iswithin)); % mean value of voxels within segment
                nof(fpc,j,hem,grad) = sum(iswithin); % number of voxels within segment
                Ydirmean(fpc,j,hem,grad) = mean(xyzmmr(2,iswithin)); % mean Y-position within segment
            end
        end % for grad
    end % for hem
end % for sub
%save ConnVsDist.mat intervall Ym nof strSub;
save ConnVsDist_20230919.mat intervall Ym nof strSub Ydirmean;

% intervall: distance from most anterior point
% Ym(subject, intervall, hemisphere, gradient)

% example for plotting left, first gradient:
% figure(35);plot(intervall(1:(end-1))+1,Ym(:,:,1,1)')



