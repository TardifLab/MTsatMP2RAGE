%%%%%% IN DEVELOPMENT - B1 FILES NOT YET ADDED %%%%%%%%%%%
%% Calculate MTsat maps from images acquired using the sample protocol provided by MNI
% 
% We assume the images have been loaded in with the following variable names:
%
% b1 - processed b1 map from fmap-b1_tfl
% MTw - MT-weighted image (anat-mpm_acq-mni_megre_MTon_1mm) - should include 6 echoes stacked in 4th dimension
% mp2r_uni -UNI image  from MP2RAGE (anat-t1w_acq-mp2rage_1mm_p3)
% mp2r_inv1 - Inversion 1 Image
% mp2r_inv2 - Inversion 2 image (anat-t1w_acq-mp2rage_1mm_p3)

% Note: the dictionary mapping approach requires the INV1 image as well. 

% Written by Christopher Rowley 2023
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

OutputDir = 'Directory\b1Correction\outputs'; % should be the same as in simSeq_M0b_R1obs_MNIprot.m
DATADIR = 'Directory/to/Save/Output/Images'; % note that we then save images into a folder called 'matlab'

% Location of B1 correction files for MTsat:
fitValues_S_1 = load(fullfile( OutputDir, 'fitValues_MP2RAGE_MTsat.mat')); 
fitValues_S_1 = fitValues_S_1.fitValues;

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
MP2RAGE.TRFLASH = 6.8e-3; % TR of the GRE readout
MP2RAGE.TIs = [0.94 2.83];% inversion times - time between middle of refocusing pulse and excitatoin of the k-space center encoding
MP2RAGE.NZslices = [172/2 172/2];%  Listed in Sequence-Special Samples/TR
MP2RAGE.NZslices = ceil(MP2RAGE.NZslices);
MP2RAGE.FlipDegrees = [4 4];% Flip angle of the two readouts in degrees
    
% For dictionary mapping, correct INV1 and INV2 images:
INV1img.img = mp2r_inv1;
INV2img.img = mp2r_inv2;
UNIimg.img = mp2r_uni; 
interpolateT1 = true;
B1img.img = b1;

[INV1img, INV2img] = Correct_INV1INV2_withMP2RAGEuni(INV1img, INV2img, UNIimg, 0);

[T1, M0map.img, R1map.img] = MP2RAGE_dictionaryMatching(MP2RAGE, real(INV1img.img), INV2img.img, B1img.img, [0.002, 0.005], 1, B1img.img ~= 0);

T1_ms = limitHandler(T1*1000, 0 , 6000);
App_mp2 = double(limitHandler(M0map.img, 0, 20000));

hdr.file_name = strcat(DATADIR,'matlab/csMP2RAGE_T1.mnc.gz'); niak_write_vol(hdr, T1_map);
hdr.file_name = strcat(DATADIR,'matlab/csMP2RAGE_M0.mnc.gz'); niak_write_vol(hdr, App_mp);


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
    title('MTsat map pre- B1 correction')
end

%% Correct for B1 (code to do so is below)
R1 = (1./T1_map)*1000; % need to convert to 1/s from 1/ms

if exist('mask', 'var')   
    R1  = R1.*mask;
end

corr_MTsat_map = MTsat_B1corr_factor_map(b1, R1, 1, fitValues_S_1);
MTsat_corr = MTsat.*( 1 + corr_MTsat_map);

if exist('mask', 'var')   
    MTsat_corr  = MTsat_corr.*mask;
end

if viewDebugImages
    figure; imshow3Dfull(MTsat_corr , [0 0.03], turbo) 
    title('MTsat map post- B1 correction')
end

hdr.file_name = strcat(DATADIR,'matlab/csMP2RAGE_MTsat', fileExt); 
niak_write_vol(hdr, MTsat_corr);
