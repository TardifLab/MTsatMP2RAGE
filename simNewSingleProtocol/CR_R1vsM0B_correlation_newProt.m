%% Generate the M0B mapping to R1 from simulation results and acquired data

% We will assume that you have loaded your images for a single subject, and
% have the following [names as contained in square brackets]:
% - [T1map] (B1+ corrected)
% - [MTsat]  image (calculated using B1+ corrected T1 and S0, and corrected
%               for excitation flip angles in the equation)
% - [B1]+ map
% - brain [mask] - reasonable tight, remove skull and other tissue

%% Load images:

DATADIR = 'Image/Directory';
OutputDir = 'Directory\b1Correction\outputs'; % should be the same as in simSeq_M0b_R1obs_MNIprot.m
seqStr = 'GRE_MTsat';

% load in the fit results from simulations
fitValues= load(fullfile(OutputDir,'fitValues_',seqStr,'.mat'));
fitValues = fitValues.fitValues;


%% Load the mask - check orientation
figure; imshow3Dfullseg( T1map, [400 3000], mask)


%% You can choose whether you want to do additional processing on the B1 map
b1 = double(b1);
b1 = limitHandler(b1, 0.5,1.4);
b1 = CR_imgaussfilt3_withMask(b1, mask1, 5); %light smoothing to the map

figure; imshow3Dfull(b1, [0.6 1.2],jet) % check result. Values should be centered around 1, not 100.

%% Mask -> threshold based on T1 to remove CSF and some dura
mask(T1map > 2200) = 0;
mask(T1map < 650) = 0;
mask(isnan(T1map)) = 0;
mask = bwareaopen(mask, 10000,6);

figure; imshow3Dfullseg(T1map, [300 2500],mask) % confirm visually


%% With MTsat maps made, perform M0b mapping
% need to convert to 1/s from 1/ms
R1  = (1./T1map) *1000;

% initialize matrices
M0b  = zeros(size(R1));

%% SPEED IT UP BY DOING A FEW AXIAL SLICES
 
ax1 = 80;
ax2 = 110;
% Check to make sure you are going through right dimension
figure; imshow3Dfullseg(T1map(:,ax1:ax2,:), [300 2500],mask(:,ax1:ax2,:))

tic %  
for i = 1:size(R1_VFA,1) % coronal
    
    for j = ax1:ax2 % 1:size(R1_VFA,2) % for axial slices
        for k =  1:size(R1_VFA,3) % sagital slices  65
            
            if mask(i,j,k) > 0 %&& dual_s(i,j,k,3) > 0
                    
                 [M0b(i,j,k), ~,  ~]  = CR_fit_M0b_v2( b1(i,j,k), R1(i,j,k), MTsat(i,j,k), fitValues);

            end
        end
    end
    
    if rem( i, 20) == 0
        % save intermediate results, since it can take a while and your
        % computer might reset...
        hdr.file_name = fullfile(OutputDir,'M0b.mnc.gz'); niak_write_vol(hdr, M0b);
    end
    disp(i)
end
toc %% this took 30hours for 1mm isotropic full brain dataset. * was running fitting in another matlab
    % instance, so could be easily sped up running on its own and/or adding
    % the parfor loop. 


figure; imshow3Dfull(M0b, [0 0.15],jet) % If the image doesn't display well
% with this colour range, then the parameters are likely wrong somewhere.

% export
mkdir(fullfile(OutputDir,'processing'))
hdr.file_name = fullfile(OutputDir,'M0b.mnc.gz'); niak_write_vol(hdr, M0b);

%% With M0B maps made, correlate with R1 and update the fitValues file. 

% use this fake mask to get rid of dura. 
tempMask = mask;
tempMask = imerode(tempMask, strel('sphere',2));
figure; imshow3Dfullseg(M0b_VFA, [0 0.15],tempMask)

mkdir(fullfile(OutputDir,'figures'));

% Optimized Approach
fitValues  = CR_generate_R1vsM0B_correlation( R1, ...
    M0b, tempMask, fitValues, fullfile(OutputDir,'figures/R1vsM0b_csMP2RAGE.png'),...
    fullfile(OutputDir,'fitValues_MTsat.mat'));


%% Now use these results to B1 correct the data:
OutputDir = DATADIR;
corr_prot = MTsat_B1corr_factor_map(b1, R1, 1, fitValues);

% Part 2, apply correction map
MTsat_c   = (MTsat + MTsat.* corr_prot) .* mask1;
MTsat_c = double(limitHandler(MTsat_c,0, 0.1));

%% View results
figure; imshow3Dfull(MTsat_c , [0 0.046], jet); 


%% Save MTsat corrected
hdr.file_name = strcat(DATADIR,'/matlab/MTsat_corrected.mnc.gz'); niak_write_vol(hdr, MTsat_c);


% If you want to save the processed B1 map, you can run the following:
hdr.file_name = strcat(DATADIR,'/matlab/b1.mnc.gz'); niak_write_vol(hdr, b1);























