%%%%%% IN DEVELOPMENT - B1 FILES NOT YET ADDED %%%%%%%%%%%
%% Calculate MTsat maps from images acquired using the sample protocol provided by MNI
% 
% We assume the images have been loaded in with the following variable names:
%
% b1 - processed b1 map from fmap-b1_tfl
% MTw - MT-weighted image (anat-mpm_acq-mni_megre_MTon_1mm) - should include 6 echoes stacked in 4th dimension
% mp2r_uni -UNI image  from MP2RAGE (anat-t1w_acq-mp2rage_1mm_p3)
% mp2r_inv2 - Inversion 2 image (anat-t1w_acq-mp2rage_1mm_p3)


% Written by Christopher Rowley 2023
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Location of B1 correction files for MTsat:
fitValues_S_1 = load(fullfile(b1Dir,'fitValues_S_1.mat')); 

% Define the output image file extension
fileExt = '.nii'; % other options include '.nii.gz', '.mnc', '.mnc.gz'

% Specify if you want to view intermediate files:
viewDebugImages = 1;

% MTw image parameters:
low_flip_angle = 5;    % flip angle in degrees -> Customize
TR1 = 25;               % low flip angle repetition time of the GRE kernel in milliseconds -> Customize
    

% ** if you run into issues here, it might be how I rescaled the UNI
% image in the function. Please contact me with any issues so I can fix
% it for everyone :) 

% Customize this parameters:
MP2RAGE.B0 = 3;           % in Tesla
MP2RAGE.TR = 5;           % MP2RAGE TR in seconds     
MP2RAGE.TRFLASH = 7.7e-3; % TR of the GRE readout
MP2RAGE.TIs = [0.94 2.83];% inversion times - time between middle of refocusing pulse and excitatoin of the k-space center encoding
MP2RAGE.NZslices = [208/2 208/2];%  should be two values, number of excitations before k-space center, and number after. [Slices Per Slab * [PartialFourierInSlice-0.5  0.5] ]
MP2RAGE.NZslices = ceil(MP2RAGE.NZslices);
MP2RAGE.FlipDegrees = [4 5];% Flip angle of the two readouts in degrees


MP2RAGEimg.img = mp2r_uni; % load_untouch_nii(MP2RAGE.filenameUNI);
MP2RAGEINV2img.img = mp2r_inv2; % load_untouch_nii(MP2RAGE.filenameINV2);
B1.img = b1;
brain.img = mask1;

tic
[ T1map, ~, M0] = CR_T1B1correctpackageTFL_withM0( B1, MP2RAGEimg, MP2RAGEINV2img, MP2RAGE, brain, 0.96);
toc
    

T1_map = T1map.img;
M0 = M0.img;

T1_map = limitHandler(T1_map, 0 , 6000);
M0 = double(limitHandler(M0, 0, 20000));
    
hdr.file_name = strcat(DATADIR,'matlab/MP2RAGE_T1', fileExt); niak_write_vol(hdr, T1_map);
hdr.file_name = strcat(DATADIR,'matlab/MP2RAGE_M0', fileExt); niak_write_vol(hdr, M0);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Calculate MTsat
% Inital Parameters

if exist('mask', 'var')   
    MTsat  = calcMTsatThruLookupTablewithDummyV3( MTw, b1, T1_map, mask,  M0/2.5, 0, 1, TR1, low_flip_angle,0);
else
    MTsat  = calcMTsatThruLookupTablewithDummyV3( MTw, b1, T1corr, 1,  M0/2.5, 0, 1, TR1, low_flip_angle,0);
end

if viewDebugImages
    figure; imshow3Dfull(MTsat , [0 0.03], turbo) 
    title('M0 map pre- B1 correction')
end

%% Correct for B1 (code to do so is below)
R1 = (1./T1_map)*1000; % need to convert to 1/s from 1/ms

if exist('mask', 'var')   
    R1  = R1.*mask;
end

corr_MTsat_map = MTsat_B1corr_factor_map(b1, R1, 3.586, fitValues_S_1);
MTsat_corr = MTsat.*( 1 + corr_MTsat_map);

if exist('mask', 'var')   
    MTsat_corr  = MTsat_corr.*mask;
end

if viewDebugImages
    figure; imshow3Dfull(MTsat_corr , [0 0.03], turbo) 
    title('M0 map post- B1 correction')
end

hdr.file_name = strcat(DATADIR,'matlab/MP2RAGE_MTsat', fileExt); 
niak_write_vol(hdr, MTsat_corr);

