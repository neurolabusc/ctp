function scalar_atlas(pth)
%report intensity for each "scalar_" image in each "atlas_"
if ~exist('pth', 'var')
    pth = pwd;
end;
%load scalars:
fnms = dir(fullfile(pth, 'scalar_*.nii'));
if isempty(fnms)
    error('Unable to find "scalar_" files in %s', pth);
end
scalars = {};
for i = 1:numel(fnms)
    scalars{end+1} = fullfile(pth, fnms(i).name);
end
%load atlases
fnms = dir(fullfile(pth, 'atlas*.nii'));
if isempty(fnms)
    error('Unable to find "atlas_" files in %s', pth);
end
atlases = {};
for i = 1:numel(fnms)
    atlases{end+1} = fullfile(pth, fnms(i).name);
end
%analyze

for a = 1 : numel(atlases)
    rfnm = atlases{a};
    rhdr = spm_vol(rfnm); %load header 
    rimg = spm_read_vols(rhdr); %load image data
    nroi = max(rimg(:)); %number of regions
    for s = 1 : numel(scalars)
        fnm = scalars{s};
        hdr = spm_vol(fnm); %load header 
        img = spm_read_vols(hdr); %load image data
        fprintf("scalar:atlas = %s:%s\n", fnm, rfnm); 
        %https://github.com/neurolabusc/NiiStat/blob/master/nii_roi2stats.m
        mn = zeros(nroi,1);
        for r = 1:nroi
            idx = (rimg(:) == r) & (~isnan(img(:)) );
            mn(r) = mean(img(idx));
            fprintf('Region %d has %d voxels with a mean of %f\n',r,sum(idx), mn(r));
        end
    end
end