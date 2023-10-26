%%%%%% IN DEVELOPMENT - B1 FILES NOT YET ADDED %%%%%%%%%%%
%% Calculate MTsat maps from images acquired using the sample protocol provided by MNI
% 
% We assume the images have been loaded in with the following variable names:
%
% b1 - processed b1 map from fmap-b1_tfl
% MTw - MT-weighted image (anat-mtsat_acq-mni_gre_MTon_1mm)
% PDw - PD-weighted image (anat-mtsat_acq-mni_gre_MToff_1mm)
% T1w - T1-weighted image (anat-mtsat_acq-mni_gre_T1w_1mm)
% mask - optional brain mask to mask data.

% Written by Christopher Rowley 2023
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

OutputDir = 'Directory\b1Correction\outputs'; % should be the same as in simSeq_M0b_R1obs_MNIprot.m
DATADIR = 'Directory/to/Save/Output/Images'; % note that we then save images into a folder called 'matlab'

% Location of B1 correction files for MTsat:
fitValues_S_1 = load(fullfile(OutputDir,'fitValues_VFA_MTsat.mat')); 
fitValues_S_1 = fitValues_S_1.fitValues;

% Define the output image file extension
fileExt = '.nii'; % other options include '.nii.gz', '.mnc', '.mnc.gz'

% Specify if you want to view intermediate files:
viewDebugImages = 1;

%% Loading and prep done, now for calculations

%% VFA Helms/ hMRI approach.
% The model-based approach in the present paper uses a B1 map and
% spoiling correction here.

low_flip_angle = 5;    % flip angle in degrees -> Customize
high_flip_angle = 15;  % flip angle in degrees -> Customize
TR1 = 25;               % low flip angle repetition time of the GRE kernel in milliseconds -> Customize
TR2 = 15;               % high flip angle repetition time of the GRE kernel in milliseconds -> Customize
smallFlipApprox = true;

a1 = deg2rad(low_flip_angle); 
a2 = deg2rad(high_flip_angle); 

R1 = hmri_calc_R1(struct( 'data', PDw, 'fa', a1, 'TR', TR1, 'B1', b1),...
            struct( 'data', T1w, 'fa', a2, 'TR', TR2, 'B1', b1), smallFlipApprox);

if exist('mask', 'var')
    R1 = R1.*mask;            
    T1 = 1./R1.*mask; % convert to milliseconds
else
    T1 = 1./R1; % convert to milliseconds
end

T1 = limitHandler(T1, 0, 6000);

if viewDebugImages
    figure; imshow3Dfull(T1, [300 2500],jet)
    title('T1 map pre- spoiling correction')
end


%% Apply spoiling correction:
% I have copied the numbers from the JSON file from hMRI spoil correction:
% I have wrote to a mat file

% Note in Preibisch and Deichmann 2009, A and B are quadratic functions
% dependent on b1map.

% Get coefficients from json:
Acoef = [30.9471,-36.7804,23.9561];
Bcoef = [-0.0723,0.0474,0.9853];

% Make correction factors:
A = Acoef(1)*b1.^2 + Acoef(2)*b1 + Acoef(3);
B = Bcoef(1)*b1.^2 + Bcoef(2)*b1 + Bcoef(3);

% Apply to the T1 map - in milliseconds as simulations are done in ms
T1corr = A + B.*T1;
T1corr = limitHandler(T1corr, 0, 6000);

if exist('mask', 'var')   
    T1corr =T1corr.*mask; % convert to milliseconds
end

if viewDebugImages
    figure; imshow3Dfull(T1corr, [300 2500],jet)
    title('T1 map post- spoiling correction')
end

hdr.file_name = strcat(DATADIR,'matlab/VFA_T1_hMRIspoil_t2_80ms', fileExt); 
niak_write_vol(hdr,T1corr);    
    

%% Recalculate M0 using this and the flash equation

flip_a = (low_flip_angle*b1) * pi / 180; % correct for B1 and convert to radians
x = cos(flip_a) ;
y = exp(-TR1./T1corr);

% Solve for M0 using equation for flass image.
M0_1 = PDw.*(1-x.*y)./ ( (1-y) .*sin(flip_a));

% second image
flip_a = (high_flip_angle*b1) * pi / 180; % correct for B1 and convert to radians
x = cos(flip_a) ;
y = exp(-TR2./T1corr);

% Solve for M0 using equation for flass image.
M0_2 = T1w.*(1-x.*y)./ ( (1-y) .*sin(flip_a));

M0_spoil = 2.5* (M0_1 + M0_2)./2; % 2.5 is the gain factor. 
M0_spoil = limitHandler(M0_spoil, 0, 20000);

if exist('mask', 'var')   
    M0_spoil =M0_spoil.*mask; % convert to milliseconds
end

if viewDebugImages
    figure; imshow3Dfull(M0_spoil , [00 15000])  
    title('M0 map post- spoiling correction')
end

hdr.file_name = strcat(DATADIR,'matlab/VFA_M0_hMRIspoil_t2_80ms', fileExt); 
niak_write_vol(hdr, M0_spoil);   

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% T1 and M0 calculated. Now calculate MTsat
if exist('mask', 'var')   
    MTsat  = calcMTsatThruLookupTablewithDummyV3( MTw, b1, T1corr, mask,  M0_spoil, 0, 1, TR1, low_flip_angle,0);
else
    MTsat  = calcMTsatThruLookupTablewithDummyV3( MTw, b1, T1corr, 1,  M0_spoil, 0, 1, TR1, low_flip_angle,0);
end

if viewDebugImages
    figure; imshow3Dfull(MTsat , [0 0.03], turbo) 
    title('M0 map pre- B1 correction')
end

R1_spoil = (1./T1corr)*1000; % need to convert to 1/s from 1/ms

if exist('mask', 'var')   
    R1_spoil  = R1_spoil.*mask;
end

corr_MTsat_map = MTsat_B1corr_factor_map(b1, R1_spoil, 1,fitValues_S_1);

MTsat_corr = MTsat.*( 1 + corr_MTsat_map);

if exist('mask', 'var')   
    MTsat_corr  = MTsat_corr.*mask;
end

if viewDebugImages
    figure; imshow3Dfull(MTsat_corr , [0 0.03], turbo) 
    title('M0 map post- B1 correction')
end

hdr.file_name = strcat(DATADIR,'matlab/VFA_MTsat_spoil', fileExt); 
niak_write_vol(hdr, MTsat_corr);
