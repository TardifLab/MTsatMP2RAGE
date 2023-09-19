% This code calculates the different maps observed in the manuscript
% Currently, I have only provided data to do the MP2RAGE calcuation.
% If you require data to test these function, let me know and I will
% provide some.

% Search for the word 'Customize' for sections that I believe will require
% some tuning for your specific input images. I have included commented out
% image viewing code at any step where I suspect where problems could
% arise. On your first run through, it is best to just step through
% line/line to check

%% For all calculations the following images are required:
% - MP2RAGE (UNI and INV2 image)
% - B1 map

%% The following toolboxes/code will be needed:
% Note: I have provided sample input code, but you will need to download
% the code and modify 'media' in the below commands to the location of the
% folders

addpath( genpath( '/media/GitHub/TardifLab/OptimizeIHMTimaging/' )); % newer B1 correciton code for MTsat
addpath( genpath( '/media/GitHub/niak-master' )); % Image load code, or use what ever you want
addpath( genpath( '/media/GitHub/NeuroImagingMatlab')); % viewer code and other utilities
addpath( genpath('/media/Github/unring-master')) % For Gibbs unringing


%% Set up for image loading
subjectNames = {'hc01' , 'hc02' 'hc03' 'hc04' }; % this will loop through a few subjects
DD = 'folder/location'; % enter filepath to images UP TO but not including the subjectName above


% load in the fit results from simWith_M0b_R1obs.m
b1Dir = '/media/GitHub/MTsatMP2RAGE/b1' ;
fitValues_S_1 = load(fullfile(b1Dir,'fitValues_S_1.mat'));
fitValues_S_1 = fitValues_S_1.fitValues;


for z = 1:length(subjectNames)
    
    disp(subjectNames{z});

    DATADIR = fullfile( DD, subjectNames{z},'/images/');

    %image names:
    mtw_fn = { 'MT2k_reg.mnc' 'MP2RAGE_inv1_reg.mnc' 'MP2RAGE_inv2_reg.mnc' 'MP2RAGE_UNI_reg.mnc'}';                                            
  
    %% If running on sample data, just load all the input matrices and skip the image loading here
    % If you skip, jump to "START HERE" (CTRL+F search)
    
    % Customize this to use whatever you image loading code you want. I
    % find this to be robust across image types
    for i = 1:length(mtw_fn)
        fn = fullfile(DATADIR,mtw_fn{i});
        [hdr, img] = niak_read_vol(fn);
        comb_mtw(:,:,:,i) = img; %.img;
    end

    comb_mtw = double(comb_mtw);

    % Use this code to check that it loaded properly
    % figure; imshow3Dfull(comb_mtw(:,:,:,5) , [0 650], jet)
    % figure; imshow3Dfull(comb_mtw(:,:,:,6) , [0 650], jet)

    input_img = comb_mtw(:,:,:,1);
    
    
    %% Make a background mask 
    mask1 = zeros(size (input_img));
    mask1( (input_img > 100) )= 1;
    %figure; imshow3Dfullseg(input_img , [0 550], mask1)


    %% Some B1 issues so lets try and load that
    [hdr2, b1] = niak_read_vol(strcat(DATADIR,'../resampled_b1field.mnc'));
    b1 = double(b1);
    b1 = limitHandler(b1, 0.5, 1.4);
   
    % Customize this and use the viewing code to make sure you get a map
    % you want, and do not remove regions of brain. I had issues with the
    % sinuses blending into the brain, so this tries to mask it out and
    % smooth over the mask

    b1_t = permute(b1,[3 1 2]);
    b1_t = flip(b1_t,2);
    b1_t = flip(b1_t,3);
    b1_t = flip(b1_t,1);

    b1_t = CR_b1_erode_dilate(b1_t .*mask1, 7);
    mask2 = ones(size (b1_t));
    mask2(b1_t == 0 )= 0;
    b1 = CR_imgaussfilt3_withMask(b1_t, mask2, 2); %light smoothing to the map
    clear b1_t;
    
    figure; imshow3Dfullseg(b1, [0.6 1.2],mask1)

    %% Run MP-PCA denoising
    all_PCAcorr = MPdenoising(comb_mtw(:,:,:,1:end-2), [],[6,6,6]); % remove T1 maps

    mt2k = all_PCAcorr(:,:,:,4);
    lfa  = all_PCAcorr(:,:,:,5);
    hfa  = all_PCAcorr(:,:,:,6);
    
    mp2r_inv1 = all_PCAcorr(:,:,:,7);
    mp2r_inv2 = all_PCAcorr(:,:,:,8);
    mp2r_uni  = all_PCAcorr(:,:,:,9);
    
    sp_mp2r_inv1 = all_PCAcorr(:,:,:,10);
    sp_mp2r_inv2 = all_PCAcorr(:,:,:,11);
    sp_mp2r_uni  = all_PCAcorr(:,:,:,12);

    % figure; imshow3Dfullseg(mt2k , [0 550],mask1)
    %figure; imshow3Dfull(dual- input_img , [-50 50]) % see that there is a difference
    tic
    lfa_ur= unring3D(lfa,3);
    hfa_ur= unring3D(hfa,3);
    mt2k_ur= unring3D(mt2k,3);
    
    mp2r_inv1 = unring3D(mp2r_inv1 ,3);
    mp2r_inv2 = unring3D(mp2r_inv2 ,3);
    mp2r_uni = unring3D(mp2r_uni,3);
    
    sp_mp2r_inv1 = unring3D(sp_mp2r_inv1 ,3);
    sp_mp2r_inv2 = unring3D(sp_mp2r_inv2 ,3);
    sp_mp2r_uni = unring3D(sp_mp2r_uni,3);
    toc

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% Loading and prep done, now for calculations - START HERE
   
    
    %% MP2RAGE 

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
    
    
%   The parameters we used for sparse MP2RAGE    
%     MP2RAGE.TRFLASH = 6.4e-3; % TR of the GRE readout
%     MP2RAGE.TIs = [0.94 2.83];% inversion times - time between middle of refocusing pulse and excitatoin of the k-space center encoding
%     MP2RAGE.NZslices = [ceil(175/2) floor(175/2)];%  should be two values, number of excitations before k-space center, and number after. [Slices Per Slab * [PartialFourierInSlice-0.5  0.5] ]
%     MP2RAGE.FlipDegrees = [4 5];% Flip angle of the two readouts in degrees 
    
    MP2RAGEimg.img = mp2r_uni; % load_untouch_nii(MP2RAGE.filenameUNI);
    MP2RAGEINV2img.img = mp2r_inv2; % load_untouch_nii(MP2RAGE.filenameINV2);
    B1.img = b1;
    brain.img = mask1;
    
    tic
    [ T1map, ~, Appmap] = CR_T1B1correctpackageTFL_withM0( B1, MP2RAGEimg, MP2RAGEINV2img, MP2RAGE, brain, 0.96);
    toc
      
    
    T1_map = T1map.img;
    App_mp = Appmap.img;
    
    T1_map = limitHandler(T1_map, 0 , 6000);
    App_mp = double(limitHandler(App_mp, 0, 20000));
    
        
    % figure; imshow3Dfull(T1_map, [300 2500],jet)
    % figure; imshow3Dfull(App_mp2 , [00 15000])
    
    hdr.file_name = strcat(DATADIR,'matlab/MP2RAGE_T1.mnc.gz'); niak_write_vol(hdr, T1_map);
    hdr.file_name = strcat(DATADIR,'matlab/MP2RAGE_M0.mnc.gz'); niak_write_vol(hdr, App_mp);
  
    
    %% Calculate MTsat
    % Inital Parameters
    TR = 27; %milliseconds
    flipA = 6; % flip angle in degrees

    % calculate maps
    mp2_sat_2k     = calcMTsatThruLookupTablewithDummyV3( mt2k_ur, b1, T1_map, mask1, App_mp/2.5, 0, 1, TR, flipA, 0);

    % view to make sure all the parameters were set correctly
    figure; imshow3Dfull(mp2_sat_2k , [0 0.03], turbo)
  
    
    
       %% Correct for B1 (code to do so is below)
    
    R1_mp2_s = (1./T1_map)*1000; % need to convert to 1/s from 1/ms


    mp2_corr_MTsat_2k   = MTsat_B1corr_factor_map(b1, R1_mp2_s, 3.586,fitValues_S_1);
    mp2_sat_2k_c   = mp2_sat_2k.*(1 + mp2_corr_MTsat_2k) .* mask;



    %% Other things, save if you want

    mp2_sat_2k_c   = limitHandler(mp2_sat_2k_c, 0 , 0.3);
    figure; imshow3Dfull(mp2_sat_2k_c , [0 0.03], turbo)


    hdr.file_name = strcat(DATADIR,'matlab/MP2RAGE_MTsat_2k.mnc.gz'); niak_write_vol(hdr, mp2_sat_2k_c);

     
    % clear variables for next loop
    clearvars -except subjectNames z b1Dir
    z = z+1;
    
    
end
