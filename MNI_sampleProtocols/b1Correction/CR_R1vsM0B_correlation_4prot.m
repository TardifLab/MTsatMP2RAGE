%% Generate the M0B mapping to R1 from simulation results and acquired data

addpath(genpath('Directory/niak-master' ))



%% Load images:

DATADIR = 'Image/Directory';
OutputDir = 'Directory\b1Correction\outputs'; % should be the same as in simSeq_M0b_R1obs_MNIprot.m

%image names:
% in the order of dual, hfa, neg, lfa, pos
mtw_fn = {'anat-mpm_acq-mni_megre_MToff_1mm.mnc.gz',...
    'anat-mpm_acq-mni_megre_MToff_1mm1.mnc.gz',...
    'anat-mpm_acq-mni_megre_MToff_1mm2.mnc.gz',...
    'anat-mpm_acq-mni_megre_MToff_1mm3.mnc.gz',...
    'anat-mpm_acq-mni_megre_MToff_1mm4.mnc.gz',...
    'anat-mpm_acq-mni_megre_MToff_1mm5.mnc.gz',...
    'anat-mpm_acq-mni_megre_MToff_1mm6.mnc.gz',...
    'anat-mpm_acq-mni_megre_MToff_1mm7.mnc.gz',...
    'anat-mpm_acq-mni_megre_MTon_1mm.mnc.gz',...
    'anat-mpm_acq-mni_megre_MTon_1mm1.mnc.gz',...
    'anat-mpm_acq-mni_megre_MTon_1mm2.mnc.gz',...
    'anat-mpm_acq-mni_megre_MTon_1mm3.mnc.gz',...
    'anat-mpm_acq-mni_megre_MTon_1mm4.mnc.gz',...
    'anat-mpm_acq-mni_megre_MTon_1mm5.mnc.gz',...
    'anat-mpm_acq-mni_megre_T1w_1mm.mnc.gz',...
    'anat-mpm_acq-mni_megre_T1w_1mm1.mnc.gz',...
    'anat-mpm_acq-mni_megre_T1w_1mm2.mnc.gz',...
    'anat-mpm_acq-mni_megre_T1w_1mm3.mnc.gz',...
    'anat-mpm_acq-mni_megre_T1w_1mm4.mnc.gz',...
    'anat-mpm_acq-mni_megre_T1w_1mm5.mnc.gz',...
    'anat-mtsat_acq-mni_gre_MToff_1mm.mnc',...
    'anat-mtsat_acq-mni_gre_MTon_1mm.mnc',...
    'anat-mtsat_acq-mni_gre_T1w_1mm.mnc',...
    'anat-t1w_acq-mp2rage_1mm_p3_INV1.mnc',...
    'anat-t1w_acq-mp2rage_1mm_p3_INV2.mnc',...
    'anat-t1w_acq-mp2rage_1mm_p3_UNI_Images.mnc',...
    'anat-t1w_acq-mp2rage_sparse_1mm_WIP925_INV1.mnc',...
    'anat-t1w_acq-mp2rage_sparse_1mm_WIP925_INV2.mnc',...
    'anat-t1w_acq-mp2rage_sparse_1mm_WIP925_UNI_Images.mnc'};
 
for i = 1:size(mtw_fn,2)
    fn = fullfile(DATADIR,mtw_fn{i});
    [hdr, img] = niak_read_vol(fn);
    comb_mtw(:,:,:,i) = img; %.img;
end


% load in the fit results from simulations
fitValues_csMP2 = load(fullfile(OutputDir,'fitValues_csMP2RAGE_MTsat.mat'));
fitValues_csMP2 = fitValues_csMP2.fitValues;

fitValues_MP2 = load(fullfile(OutputDir,'fitValues_MP2RAGE_MTsat.mat'));
fitValues_MP2 = fitValues_MP2.fitValues;

fitValues_meGRE= load(fullfile(OutputDir,'fitValues_meGRE_MTsat.mat'));
fitValues_meGRE = fitValues_meGRE.fitValues;

fitValues_VFA= load(fullfile(OutputDir,'fitValues_VFA_MTsat.mat'));
fitValues_VFA = fitValues_VFA.fitValues;


%% Load the mask
fn = fullfile(DATADIR,'itk_mask.nii.gz');
[~, mask] = niak_read_vol(fn);
mask1 = permute(mask,[2 3 1]); % conversion between minc and nii reorients it

figure; imshow3Dfullseg( comb_mtw(:,:,:,i), [400 5000], mask1)

%% Some B1 issues so lets try and load that
[~, b1] = niak_read_vol(fullfile(DATADIR,'B1map_rs.mnc')); 
b1 = double(b1);
b1 = limitHandler(b1, 0.5,1.4);
b1 = CR_imgaussfilt3_withMask(b1, mask1, 5); %light smoothing to the map

figure; imshow3Dfull(b1, [0.6 1.2],jet)

%% Run MP-PCA denoising
comb_mtw = double(comb_mtw);
all_PCAcorr = MPdenoising(comb_mtw);


%% separate the images then average the MTw ones
mgre_PDw = all_PCAcorr(:,:,:,1:8);
mgre_MTw  = all_PCAcorr(:,:,:,9:14);
mgre_T1w = all_PCAcorr(:,:,:,15:20);

gre_PDw = comb_mtw(:,:,:,21);
gre_MTw  = comb_mtw(:,:,:,22);
gre_T1w = comb_mtw(:,:,:,23);

mp2r_inv1  = all_PCAcorr(:,:,:,24);
mp2r_inv2 = all_PCAcorr(:,:,:,25);
mp2r_uni  = all_PCAcorr(:,:,:,26);

sp_mp2r_inv1  = all_PCAcorr(:,:,:,27);
sp_mp2r_inv2 = all_PCAcorr(:,:,:,28);
sp_mp2r_uni  = all_PCAcorr(:,:,:,29);

%%%%%%%%%%%%%%%%%%%%%%%%%
%% Calculate maps
%%%%%%%%%%%%%%%%%%%%%%%%%

%% MP2RAGE:
  
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
    
mkdir( fullfile(DATADIR, 'matlab'));
fileExt = '.mnc.gz'; % other options include '.nii.gz', '.mnc', '.mnc.gz'

hdr.file_name = strcat(DATADIR,'/matlab/MP2RAGE_T1', fileExt); niak_write_vol(hdr, T1_map);
hdr.file_name = strcat(DATADIR,'/matlab/MP2RAGE_M0', fileExt); niak_write_vol(hdr, M0);

% MTw image parameters:
low_flip_angle = 5;    % flip angle in degrees -> Customize
TR1 = 25;               % low flip angle repetition time of the GRE kernel in milliseconds -> Customize
mp2r_MTsat  = calcMTsatThruLookupTablewithDummyV3( gre_MTw, b1, T1_map, mask1,  M0/2.5, 0, 1, TR1, low_flip_angle,0);

figure; imshow3Dfull(T1_map, [300 2500], turbo)
figure; imshow3Dfull(M0, [0 15000], turbo)
figure; imshow3Dfull(mp2r_MTsat, [0 0.03])


%% csMP2RAGE:
  
MP2RAGE.B0 = 3;           % in Tesla
MP2RAGE.TR = 5;           % MP2RAGE TR in seconds     
MP2RAGE.TRFLASH = 6.8e-3; % TR of the GRE readout
MP2RAGE.TIs = [0.94 2.83];% inversion times - time between middle of refocusing pulse and excitatoin of the k-space center encoding
MP2RAGE.NZslices = [172/2 172/2];%  Listed in Sequence-Special Samples/TR
MP2RAGE.NZslices = ceil(MP2RAGE.NZslices);
MP2RAGE.FlipDegrees = [4 4];% Flip angle of the two readouts in degrees

MP2RAGEimg.img = sp_mp2r_uni; % load_untouch_nii(MP2RAGE.filenameUNI);
MP2RAGEINV2img.img = sp_mp2r_inv2; % load_untouch_nii(MP2RAGE.filenameINV2);
B1.img = b1;
brain.img = mask1;

tic
[ spT1map, spMP2RAGEcorrected, spAppmap2] = CR_T1B1correctpackageTFL_withM0( B1, MP2RAGEimg, MP2RAGEINV2img, MP2RAGE, brain, 0.96);
toc

spT1_map = spT1map.img;
spApp_mp2 = spAppmap2.img ;

spT1_map = limitHandler(spT1_map, 0, 6000);
spApp_mp2 = double(limitHandler(spApp_mp2, 0, 20000));


hdr.file_name = fullfile(DATADIR,'matlab/csMP2RAGE_T1.mnc.gz'); niak_write_vol(hdr, spT1_map);
hdr.file_name = fullfile(DATADIR,'matlab/csMP2RAGE_M0.mnc.gz'); niak_write_vol(hdr, spApp_mp2);

% MTw image parameters:
low_flip_angle = 5;    % flip angle in degrees -> Customize
TR1 = 25;               % low flip angle repetition time of the GRE kernel in milliseconds -> Customize
spmp2r_MTsat  = calcMTsatThruLookupTablewithDummyV3( gre_MTw, b1, spT1_map, mask1,  spApp_mp2/2.5, 0, 1, TR1, low_flip_angle,0);

figure; imshow3Dfull(spT1_map, [300 2500], turbo)
figure; imshow3Dfull(spApp_mp2, [0 15000], turbo)
figure; imshow3Dfull(spmp2r_MTsat, [0 0.03])

%% VFA

low_flip_angle = 5;    % flip angle in degrees -> Customize
high_flip_angle = 15;  % flip angle in degrees -> Customize
TR1 = 25;               % low flip angle repetition time of the GRE kernel in milliseconds -> Customize
TR2 = 15;               % high flip angle repetition time of the GRE kernel in milliseconds -> Customize
smallFlipApprox = true;

a1 = deg2rad(low_flip_angle); 
a2 = deg2rad(high_flip_angle); 

gre_R1 = hmri_calc_R1(struct( 'data', gre_PDw, 'fa', a1, 'TR', TR1, 'B1', b1),...
            struct( 'data', gre_T1w, 'fa', a2, 'TR', TR2, 'B1', b1), smallFlipApprox);

gre_T1 = 1./gre_R1.*mask1; % convert to milliseconds
gre_T1 = double(limitHandler(gre_T1, 0, 6000));

% Coeff from hMRI toolbox.
Acoef = [30.5365,-39.2652,24.6445];
Bcoef = [-0.0738,0.0576,0.98];

% Make correction factors:
A = Acoef(1)*b1.^2 + Acoef(2)*b1 + Acoef(3);
B = Bcoef(1)*b1.^2 + Bcoef(2)*b1 + Bcoef(3);

% Apply to the T1 map - in milliseconds as simulations are done in ms
gre_T1 = double(A + B.*gre_T1 .*mask1);

% Recalculate M0 using this and the flash equation
flip_a = (low_flip_angle*b1) * pi / 180; % correct for B1 and convert to radians
x = cos(flip_a) ;
y = exp(-TR1./gre_T1);

% Solve for M0 using equation for flass image.
M0_1 = gre_PDw.*(1-x.*y)./ ( (1-y) .*sin(flip_a));

% second image
flip_a = (high_flip_angle*b1) * pi / 180; % correct for B1 and convert to radians
x = cos(flip_a) ;
y = exp(-TR2./gre_T1);

% Solve for M0 using equation for flass image.
M0_2 = gre_T1w.*(1-x.*y)./ ( (1-y) .*sin(flip_a));
gre_M0 = (M0_1 + M0_2)./2; 
gre_M0 = double(limitHandler(gre_M0, 0, 20000));

hdr.file_name = fullfile(DATADIR,'matlab/gre_T1.mnc.gz'); niak_write_vol(hdr, gre_T1);
hdr.file_name = fullfile(DATADIR,'matlab/gre_M0.mnc.gz'); niak_write_vol(hdr, gre_M0);


gre_MTsat  = calcMTsatThruLookupTablewithDummyV3( gre_MTw, b1, gre_T1, mask1,  gre_M0, 0, 1, TR1, low_flip_angle,0);

figure; imshow3Dfull(gre_T1, [300 2500], turbo)
figure; imshow3Dfull(gre_M0, [0 15000], turbo)
figure; imshow3Dfull(gre_MTsat, [0 0.03])


%% ME-GRE
avgre_MTw = mean(mgre_MTw(:,:,:,1:6), 4);
avgre_PDw = mean(mgre_PDw(:,:,:,1:6), 4);
avgre_T1w = mean(mgre_T1w(:,:,:,1:6), 4);


low_flip_angle = 6;    % flip angle in degrees -> Customize
high_flip_angle = 20;  % flip angle in degrees -> Customize
TR1 = 23;               % low flip angle repetition time of the GRE kernel in milliseconds -> Customize
TR2 = 18;               % high flip angle repetition time of the GRE kernel in milliseconds -> Customize

smallFlipApprox = false;

a1 = deg2rad(low_flip_angle); 
a2 = deg2rad(high_flip_angle); 

megre_R1 = hmri_calc_R1(struct( 'data', avgre_PDw, 'fa', a1, 'TR', TR1, 'B1', b1),...
            struct( 'data', avgre_T1w, 'fa', a2, 'TR', TR2, 'B1', b1), smallFlipApprox);

megre_T1 = 1./megre_R1.*mask1; % convert to milliseconds
megre_T1 = double(limitHandler(megre_T1, 0, 6000));

% Coeff from hMRI toolbox.
Acoef = [62.5025,-72.3755,34.616];
Bcoef = [-0.1226,0.0852,0.9743000000000001];
% Make correction factors:
A = Acoef(1)*b1.^2 + Acoef(2)*b1 + Acoef(3);
B = Bcoef(1)*b1.^2 + Bcoef(2)*b1 + Bcoef(3);

% Apply to the T1 map - in milliseconds as simulations are done in ms
megre_T1 = double(A + B.*megre_T1 .*mask1);

% Recalculate M0 using this and the flash equation

flip_a = (low_flip_angle*b1) * pi / 180; % correct for B1 and convert to radians
x = cos(flip_a) ;
y = exp(-TR1./megre_T1);

% Solve for M0 using equation for flass image.
M0_1 = avgre_PDw.*(1-x.*y)./ ( (1-y) .*sin(flip_a));

% second image
flip_a = (high_flip_angle*b1) * pi / 180; % correct for B1 and convert to radians
x = cos(flip_a) ;
y = exp(-TR2./megre_T1);

% Solve for M0 using equation for flass image.
M0_2 = avgre_T1w.*(1-x.*y)./ ( (1-y) .*sin(flip_a));
megre_M0 = (M0_1 + M0_2)./2; 
megre_M0 = limitHandler(megre_M0, 0, 20000);

hdr.file_name = fullfile(DATADIR,'matlab/me-gre_T1.mnc.gz'); niak_write_vol(hdr, megre_T1);
hdr.file_name = fullfile(DATADIR,'matlab/me-gre_M0.mnc.gz'); niak_write_vol(hdr, megre_M0);

megre_MTsat  = calcMTsatThruLookupTablewithDummyV3( avgre_MTw, b1, megre_T1, mask1,  megre_M0, 0, 1, TR1, low_flip_angle,0);

figure; imshow3Dfull(megre_T1, [300 2500], turbo)
figure; imshow3Dfull(megre_M0, [0 15000], turbo)
figure; imshow3Dfull(megre_MTsat, [0 0.03])

%% List of T1 variables: 
% T1_map
% M0
% mp2r_MTsat
% 
% spT1_map
% spApp_mp2
% spmp2r_MTsat
% 
% gre_T1
% gre_M0
% gre_MTsat
% 
% megre_T1
% megre_M0
% megre_MTsat



%% Mask -> bet result touched up in itk, then threshold CSF and some dura

mask = mask1;
mask(spT1_map > 2200) = 0;
mask(spT1_map < 650) = 0;
mask(isnan(spT1_map)) = 0;
mask = bwareaopen(mask, 10000,6);
figure; imshow3Dfullseg(spT1_map, [300 2500],mask)


%% With MTsat maps made, perform M0b mapping

% need to convert to 1/s from 1/ms
R1_csMP2  = (1./spT1_map) *1000;
R1_MP2      = (1./T1_map) *1000;
R1_meGRE = (1./megre_T1) *1000;
R1_VFA      = (1./gre_T1) *1000;

% initialize matrices
M0b_csMP2  = zeros(size(R1_VFA));
M0b_MP2      = zeros(size(R1_VFA));
M0b_meGRE = zeros(size(R1_VFA));
M0b_VFA      = zeros(size(R1_VFA));



%% SPEED IT UP BY DOING A FEW AXIAL SLICES
 
ax1 = 80;
ax2 = 110;
% Check to make sure you are going through right dimension
figure; imshow3Dfullseg(spT1_map(:,ax1:ax2,:), [300 2500],mask(:,ax1:ax2,:))



tic %  
for i = 1:size(R1_VFA,1) % coronal
    
    for j = ax1:ax2 % 1:size(R1_VFA,2) % for axial slices
        for k =  1:size(R1_VFA,3) % sagital slices  65
            
            if mask(i,j,k) > 0 %&& dual_s(i,j,k,3) > 0
                    
                 [M0b_csMP2(i,j,k), ~,  ~]  = CR_fit_M0b_v2( b1(i,j,k), R1_csMP2(i,j,k), spmp2r_MTsat(i,j,k), fitValues_csMP2);
                 [M0b_MP2(i,j,k),  ~,  ~]     = CR_fit_M0b_v2( b1(i,j,k), R1_MP2(i,j,k),    mp2r_MTsat(i,j,k), fitValues_MP2);               
                 [M0b_meGRE(i,j,k),  ~,  ~] = CR_fit_M0b_v2( b1(i,j,k), R1_meGRE(i,j,k), megre_MTsat(i,j,k), fitValues_meGRE);
                 [M0b_VFA(i,j,k), ~, ~]        = CR_fit_M0b_v2( b1(i,j,k), R1_VFA(i,j,k),      gre_MTsat(i,j,k),fitValues_VFA);

            end
        end
    end
    
    if rem( i, 20) == 0
        % save intermediate results, since it can take a while and your
        % computer might reset...
        hdr.file_name = fullfile(OutputDir,'M0b_csMP2.mnc.gz'); niak_write_vol(hdr, M0b_csMP2);
        hdr.file_name = fullfile(OutputDir,'M0b_MP2.mnc.gz'); niak_write_vol(hdr,M0b_MP2);
        hdr.file_name = fullfile(OutputDir,'M0b_meGRE.mnc.gz'); niak_write_vol(hdr,M0b_meGRE);
        hdr.file_name = fullfile(OutputDir,'M0b_VFA.mnc.gz'); niak_write_vol(hdr,M0b_VFA);
    end
    disp(i)
end
toc %% this took 30hours for 1mm isotropic full brain dataset. * was running fitting in another matlab
    % instance, so could be easily sped up running on its own and/or adding
    % the parfor loop. 


figure; imshow3Dfull(M0b_csMP2, [0 0.15],jet)
figure; imshow3Dfull(M0b_MP2, [0 0.15],jet)
figure; imshow3Dfull(M0b_meGRE, [0 0.15],jet)
figure; imshow3Dfull(M0b_VFA, [0 0.15],jet)

% export
mkdir(fullfile(OutputDir,'processing'))
hdr.file_name = fullfile(OutputDir,'M0b_csMP2.mnc.gz'); niak_write_vol(hdr, M0b_csMP2);
hdr.file_name = fullfile(OutputDir,'M0b_MP2.mnc.gz'); niak_write_vol(hdr,M0b_MP2);
hdr.file_name = fullfile(OutputDir,'M0b_meGRE.mnc.gz'); niak_write_vol(hdr,M0b_meGRE);
hdr.file_name = fullfile(OutputDir,'M0b_VFA.mnc.gz'); niak_write_vol(hdr,M0b_VFA);




%% With M0B maps made, correlate with R1 and update the fitValues file. 

% use this fake mask to get rid of dura. 
tempMask = mask;
tempMask = imerode(tempMask, strel('sphere',2));
figure; imshow3Dfullseg(M0b_VFA, [0 0.15],tempMask)

mkdir(fullfile(OutputDir,'figures'));

% Optimized Approach
fitValues_csMP2  = CR_generate_R1vsM0B_correlation( R1_csMP2, ...
    M0b_csMP2, tempMask, fitValues_csMP2, fullfile(OutputDir,'figures/R1vsM0b_csMP2RAGE.png'),...
    fullfile(OutputDir,'fitValues_csMP2RAGE_MTsat.mat'));

fitValues_MP2 = CR_generate_R1vsM0B_correlation( R1_MP2, ...
    M0b_MP2, tempMask, fitValues_MP2, fullfile(OutputDir,'figures/R1vsM0b_MP2RAGE.png'),...
    fullfile(OutputDir,'fitValues_MP2RAGE_MTsat.mat'));

fitValues_meGRE = CR_generate_R1vsM0B_correlation( R1_meGRE, ...
    M0b_meGRE, tempMask, fitValues_meGRE, fullfile(OutputDir,'figures/R1vsM0b_meGRE.png'),...
    fullfile(OutputDir,'fitValues_meGRE_MTsat.mat'));

fitValues_VFA  = CR_generate_R1vsM0B_correlation( R1_VFA, ...
    M0b_VFA, tempMask, fitValues_VFA, fullfile(OutputDir,'figures/R1vsM0b_VFA.png'), ...
    fullfile(OutputDir,'fitValues_VFA_MTsat.mat'));


%% Now use these results to B1 correct the data:
OutputDir = DATADIR;

corr_prot1 = MTsat_B1corr_factor_map(b1, R1_csMP2, 1, fitValues_csMP2);
corr_prot2 = MTsat_B1corr_factor_map(b1, R1_MP2,     1, fitValues_MP2);
corr_prot3 = MTsat_B1corr_factor_map(b1, R1_meGRE, 1, fitValues_meGRE);
corr_prot4 = MTsat_B1corr_factor_map(b1, R1_VFA,      1, fitValues_VFA);


% Part 2, apply correction map
spmp2r_MTsat_c   = (spmp2r_MTsat + spmp2r_MTsat.* corr_prot1) .* mask1;
mp2r_MTsat_c   = (mp2r_MTsat + mp2r_MTsat.* corr_prot2) .* mask1;
megre_MTsat_c  = (megre_MTsat + megre_MTsat.* corr_prot3) .* mask1;
gre_MTsat_c       = (gre_MTsat     + gre_MTsat.* corr_prot4) .* mask1;



spmp2r_MTsat_c = double(limitHandler(spmp2r_MTsat_c,0, 0.1));
mp2r_MTsat_c = double(limitHandler(mp2r_MTsat_c,0, 0.1));
megre_MTsat_c = double(limitHandler(megre_MTsat_c,0, 0.1));
gre_MTsat_c = double(limitHandler(gre_MTsat_c,0, 0.1));


%% View results
figure; imshow3Dfull(spmp2r_MTsat_c , [0 0.046], jet); 
figure; imshow3Dfull(mp2r_MTsat_c , [0 0.046], jet); 
figure; imshow3Dfull(megre_MTsat_c , [0 0.03], jet)
figure; imshow3Dfull(gre_MTsat_c , [0 0.046], jet); 


%% Other things, save if you want

hdr.file_name = strcat(DATADIR,'/matlab/MTsat_csMP2RAGE.mnc.gz'); niak_write_vol(hdr, spmp2r_MTsat_c);
hdr.file_name = strcat(DATADIR,'/matlab/MTsat_MP2RAGE.mnc.gz'); niak_write_vol(hdr, mp2r_MTsat_c);
hdr.file_name = strcat(DATADIR,'/matlab/MTsat_me-gre.mnc.gz'); niak_write_vol(hdr, megre_MTsat_c);
hdr.file_name = strcat(DATADIR,'/matlab/MTsat_gre.mnc.gz'); niak_write_vol(hdr, gre_MTsat_c);

hdr.file_name = strcat(DATADIR,'/matlab/b1.mnc.gz'); niak_write_vol(hdr, b1);























