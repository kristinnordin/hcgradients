clear all;
fclose all;
close all;


%  extreme points of hippocampus noted after visual inspection
LStartPointmm = [-26 -6 -28];
LEndPointmm = [-26 -40 4];
RStartPointmm = [26 -6 -28];
REndPointmm = [26 -40 4];


Ldirection = LEndPointmm-LStartPointmm; % vector from anterior to posterior on left side
Ldirection = Ldirection./norm(Ldirection); % unity vector on left side
Rdirection = REndPointmm-RStartPointmm; % vector from anterior to posterior on right side
Rdirection = Rdirection./norm(Rdirection); % unity vector on right side


for i=1:3
    figure(i);clf;hold on; % patch surface figure with colouring according to values on line between surface and middele axis
    figure(10+i);clf;hold on; % patch figure with squares around each voxel, coloured according to their values
    figure(20+i);clf;hold on; % path surface figure with colour according to values just beneath the surface
    figure(30+i);clf;hold on; % plot figure with all voxels coloured with their values
end
for hem = 1:2 % first left, then right hemisphere
    
    % load group average subject in MNI-space
    if hem==1
        stem = '/mnt/Data/Projects/DYNAMIC/congrads/congrads.v8.movie/group_all/lh.hippocampus.norm/';
        Vh = spm_vol([stem 'lh.hippocampus.cmaps.nii']);
        SP = LStartPointmm;
        EP = LEndPointmm;
        dirre = Ldirection; % unit vector along anterior-posterior-axis
    else
        stem = '/mnt/Data/Projects/DYNAMIC/congrads/congrads.v8.movie/group_all/rh.hippocampus.norm/';
        Vh = spm_vol([stem 'rh.hippocampus.cmaps.nii']);
        SP = RStartPointmm;
        EP = REndPointmm;
        dirre = Rdirection; % unit vector along anterior-posterior-axis
    end
    
    [Yh,XYZmm] = spm_read_vols(Vh);
    for gr=1:3 % loop over gradients
        Ys = Yh(:,:,:,gr);

        % find surface (faces and vertices) of HC, defined by having a value >0.01
        is = isosurface(reshape(XYZmm(1,:),size(Ys)),reshape(XYZmm(2,:),size(Ys)),reshape(XYZmm(3,:),size(Ys)),Ys,0.01, Ys*255);
        % select voxels within HC, defined by having a value >0.01
        within = Ys>0.01;
        XYZwithin = XYZmm(:,within);
        Ywithin = Ys(within);

        vJet = parula(1000); % set colour map
        clear colval* mFaceCols* dst*;
        for i=1:size(is.faces,1) % for each face (in the surface)
            verts = is.faces(i,:);
            coords(1,:) = is.vertices(verts(1),:);
            coords(2,:) = is.vertices(verts(2),:);
            coords(3,:) = is.vertices(verts(3),:);
            v12 = coords(1,:)-coords(2,:);
            v13 = coords(1,:)-coords(3,:);
            nrm = cross(v12,v13); % Normal vector to this face
            nrm = nrm(:); % make column-vector
            mp = mean(coords,1); % middle coordinate of face
            mp = mp(:); % make column-vector

            normpoint3 = mp+nrm/norm(nrm)*3; %3 mm into/outfrom HC-surface
            dst = sqrt((XYZwithin(1,:)-normpoint3(1)).^2+(XYZwithin(2,:)-normpoint3(2)).^2+(XYZwithin(3,:)-normpoint3(3)).^2); % distance from each voxel in HC to the norm
            mn = min(dst);
            if mn>2 % if distans to HC is larger than 2mm (=norm goes out from HC) then flip it
                nrm=-nrm;
            end
            
            % colour according to values just beneath the surface
            normvect = nrm/norm(nrm)*4; %4 mm into HC
            bollmitt = (mp+(mp+normvect))/2;
            for j=1:size(XYZwithin,2)
                dst(j) = sqrt((XYZwithin(1,j)-bollmitt(1)).^2+(XYZwithin(2,j)-bollmitt(2)).^2+(XYZwithin(3,j)-bollmitt(3)).^2); % distans from each point in HC to line between mp and mp+normvect
            end
            ToUse = dst<=4; % select voxels in a cylinder with 4mm radius between mp+normvect and mp
            values = Ywithin(ToUse);

            if any(values<=0) % value should be larger than 0, otherwise stop with error
                alskfjklasjf;
            end
            colvalNorm(i) = nanmedian(values); % take median (or mean) value, should be distributed between 0 and 1
            mFaceColsNorm(i,:) = vJet(ceil(colvalNorm(i)*1000),:); % Distribute colours of colormap
            
            % colour according to values on line between surface and middle axis
            cvec = mp-SP(:); % vector from anterior position to middle of face
            pdir = cvec/norm(cvec); % unit vector of cvec
            alfa = acos(dot(pdir, dirre)); % angle between cvec and dirre
            axlangd = cos(alfa)*norm(cvec); % length from anterior position to point on axis perpendicular to middle coord of face
            lpvec = axlangd*dirre; % vector from anterior position to perp-point
            pl = lpvec+SP; % perp-point
            pl = pl(:);
            for j=1:size(XYZwithin,2)
                dst(j) = norm(cross(pl-mp, XYZwithin(:,j)-mp))/norm(pl-mp); % distans from each voxel in HC to line between mp and pl
                dst2(j) = norm(XYZwithin(:,j)-mp);
            end
            ToUse = dst<=2 & dst2<=norm(pl-mp); % select voxels in a cylinder with 2mm radius between pl and mp
            if sum(ToUse)==0
                ToUse = dst<=2;
            end
            colval(i) = median(Ywithin(ToUse)); % take median (or mean) value, should be distributed between 0 and 1
            mFaceCols(i,:) = vJet(ceil(colval(i)*1000),:); % Distribute colours of colormap
        end
        is.facevertexcdata = mFaceCols; % colour according to values on line between surface and middele axis
        figure(gr);
        % display the figure
        hp = patch('Faces', is.faces, 'Vertices', is.vertices, 'EdgeColor', 'none', 'FaceColor', 'flat', 'FaceVertexCData', is.facevertexcdata, 'FaceAlpha', 0.5);
        axis equal;
        set(gca, 'color', [1 1 1]);
        set(gca, 'XLim', [-40 40], 'YLim', [-50 0], 'ZLim', [-30 10]);
        view(-170, 20);
        title(sprintf('gradient %d', gr'));
        xlabel('X');ylabel('Y');zlabel('Z');
        
        figure(20+gr);
        is.facevertexcdata = mFaceColsNorm; % colour according to values just beneath the surface
        % display the figure
        hp = patch('Faces', is.faces, 'Vertices', is.vertices, 'EdgeColor', 'none', 'FaceColor', 'flat', 'FaceVertexCData', is.facevertexcdata, 'FaceAlpha', 0.8);
        axis equal;
        set(gca, 'color', [1 1 1]);
        set(gca, 'XLim', [-40 40], 'YLim', [-50 0], 'ZLim', [-30 10]);
        view(-170, 20);
        title(sprintf('gradient %d', gr'));
        xlabel('X');ylabel('Y');zlabel('Z');
        
        % colour each cube
        figure(10+gr);
        cubeseq = [-1 1 1 -1; 1 1 -1 -1; 1 1 1 1];
        for i=1:numel(Ywithin)
            for j=1:3
                for k=-1:2:1
                    tcubeseq = cubeseq;
                    tcubeseq(3,:) = cubeseq(3,:)*k;
                    xseq = tcubeseq(j,:);
                    tmp = j+1;
                    if tmp>3
                        tmp=tmp-3;
                    end
                    yseq = tcubeseq(tmp,:);
                    tmp = j+2;
                    if tmp>3
                        tmp=tmp-3;
                    end
                    zseq = tcubeseq(tmp,:);
                    patch(XYZwithin(1,i)+xseq, XYZwithin(2,i)+yseq, XYZwithin(3,i)+zseq, vJet(ceil(Ywithin(i)*1000),:), 'EdgeColor', 'none', 'FaceColor', 'flat', 'FaceAlpha', 0.4, 'FaceVertexCData', repmat(vJet(ceil(Ywithin(i)*1000),:),[4 1]));
                end
            end
        end
        axis equal;
        set(gca, 'color', [1 1 1]);
        set(gca, 'XLim', [-40 40], 'YLim', [-50 0], 'ZLim', [-30 10]);
        view(-170, 20);
        title(sprintf('gradient %d', gr'));
        xlabel('X');ylabel('Y');zlabel('Z');
 
        figure(30+gr);
        for i=1:numel(Ywithin)
            plot3(XYZwithin(1,i), XYZwithin(2,i), XYZwithin(3,i), 'marker', '.', 'MarkerSize', 30, 'Color', vJet(ceil(Ywithin(i)*1000),:));
        end
        axis equal;
        set(gca, 'color', [1 1 1]);
        set(gca, 'XLim', [-40 40], 'YLim', [-50 0], 'ZLim', [-30 10]);
        view(-170, 20);
        title(sprintf('gradient %d', gr'));
        xlabel('X');ylabel('Y');zlabel('Z');
 
    end % for gr=1:3
end % for hem=1:2

