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
% - VFA data (low flip angle and high flip angle images)
% - MP2RAGE (UNI and INV2 image)
% - sparse (cs)MP2RAGE (UNI and INV2 image)
% - B1 map

%% The following toolboxes/code will be needed:
% Note: I have provided sample input code, but you will need to download
% the code and modify 'media' in the below commands to the location of the
% folders

addpath( genpath( '/media/GitHub/TardifLab/OptimizeIHMTimaging/' )); % newer B1 correciton code for MTsat

addpath( genpath( '/media/GitHub/niak-master' )); % Image load code, or use what ever you want
addpath( genpath( '/media/GitHub/NeuroImagingMatlab'));
addpath( genpath( '/media/GitHub/spm12' )); % For hMRI toolbox function
addpath( genpath( '/media/GitHub/hMRI-toolbox'  )); % VFA spoiling correction and calculation
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
    mtw_fn = { 'MT2k_reg.mnc' 'lfa_reg.mnc' 'hfa_reg.mnc' ...                                % MTsat and VFA
          'MP2RAGE_inv1_reg.mnc' 'MP2RAGE_inv2_reg.mnc' 'MP2RAGE_UNI_reg.mnc'...           % Full MP2RAGE
          'spMP2RAGE_inv1_reg.mnc' 'spMP2RAGE_inv2_reg.mnc' 'spMP2RAGE_UNI_reg.mnc'...     % Sparse MP2RAGE
          'MP2RAGE_T1_reg.mnc', 'spMP2RAGE_T1_reg.mnc' }';                                            
  
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
    [hdr2, b1] = niak_read_vol(strcat(DATADIR,'../register/resampled_b1field.mnc'));
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
    %% Loading and prep done, now for calculations
    
    %% VFA Helms/ hMRI approach.
    % note the Helms MTsat approach requires no B1 map at this stage!
    % The model-based approach in the present paper uses a B1 map and
    % spoiling correction here.
    
    % No B1 correction
    low_flip_angle = 6;    % flip angle in degrees -> Customize
    high_flip_angle = 20;  % flip angle in degrees -> Customize
    TR1 = 27;               % low flip angle repetition time of the GRE kernel in milliseconds -> Customize
    TR2 = 15;               % high flip angle repetition time of the GRE kernel in milliseconds -> Customize

    
    [R1_h, App_h] = CR_fit_R1_Aapp_Helms2008( ...
    lfa, hfa, low_flip_angle, high_flip_angle, TR1, TR2, 1);

    R1_h = limitHandler(R1_h, 0 , 5);
    App_h = double(limitHandler(App_h, 0, 20000));

    MTsat_nob1 = CR_fitMTsat_Helms2008( ...
        R1_h, App_h, mt2k, low_flip_angle, TR1, 1);

    % With B1 correction
    MTsat_withB1 = CR_fitMTsat_Helms2008( ...
        R1_h, App_h, mt2k, low_flip_angle, TR1, b1);
    
    figure; imshow3Dfull(MTsat_nob1 , [0 0.03], turbo)
    figure; imshow3Dfull(MTsat_withB1 , [0 0.03], turbo)
    
    
    % Save outputs:
    T1_h = limitHandler( 1./R1_h, 0 , 6000);
    MTsat_nob1 = limitHandler( MTsat_nob1, 0 , 0.05);
    MTsat_withB1 = limitHandler( MTsat_withB1, 0 , 0.05);
    
    figure; imshow3Dfull(T1_h , [0 2500], turbo)
   
    hdr.file_name = strcat(DATADIR,'matlab/VFA_M0_noB1.mnc.gz'); niak_write_vol(hdr, mask1.*App_h);
    hdr.file_name = strcat(DATADIR,'matlab/VFA_T1_noB1.mnc.gz'); niak_write_vol(hdr, mask1.*T1_h);
    hdr.file_name = strcat(DATADIR,'matlab/VFA_MTsat_noB1.mnc.gz'); niak_write_vol(hdr, MTsat_nob1);
    hdr.file_name = strcat(DATADIR,'matlab/VFA_MTsat_fullHelms.mnc.gz'); niak_write_vol(hdr, MTsat_withB1);
    
    
    %% VFA model-based approach
        
        a1 = low_flip_angle*pi/180 .* b1; % note the inclusion of b1 here.
        a2 = high_flip_angle*pi/180 .* b1;
        
        % New code Aug 4, 2021 CR for two TR's
        R1 = 0.5 .* (hfa_ur.*a2./ TR2 - lfa_ur.*a1./TR1) ./ (lfa_ur./(a1) - hfa_ur./(a2));
        App = lfa_ur .* hfa_ur .* (TR1 .* a2./a1 - TR2.* a1./a2) ./ (hfa_ur.* TR1 .*a2 - lfa_ur.* TR2 .*a1);
    
        R1 = R1.*mask1;
        T1 = 1./R1.*mask1;
        App = App .* mask1;
        
        App = limitHandler(App, 0, 10000);
        T1 = limitHandler(T1, 0, 6000);

     %% Incomplete spoiling:
     
    param.outdir = b1Dir; % 
    
    T2 = [30, 50, 80, 110]; % test a few T2 values
    smallFlipApprox = false;
    
    param.prot_name  = strcat('vfa_t2_', num2str(T2(3)));    % Used in output
    fnJSON = fullfile(param.outdir,[strrep(param.prot_name,' ',''),'.json']);
    
    a1 = deg2rad(low_flip_angle); 
    a2 = deg2rad(high_flip_angle); 
    
    R1 = hmri_calc_R1(struct( 'data', lfa_ur, 'fa', a1, 'TR', TR1, 'B1', b1),...
               struct( 'data', hfa_ur, 'fa', a2, 'TR', TR2, 'B1', b1), smallFlipApprox);
    R1 = R1.*mask;
    
    
    App = hmri_calc_A(struct( 'data', lfa_ur, 'fa', a1, 'TR', TR1, 'B1', b1),...
               struct( 'data', hfa_ur, 'fa', a2, 'TR', TR2, 'B1', b1), smallFlipApprox);
    
    T1 = 1./R1.*mask ; % convert to milliseconds
    App = App .* mask;

    App = limitHandler(App, 0, 20000);
    T1 = limitHandler(T1, 0, 6000);

    % Make one to match the TFL (MP2RAGE) sequence
    App_tfl = App *2.5; % match the gain

    figure; imshow3Dfull(T1.* mask1, [300 2500],jet)


    hdr.file_name = strcat(DATADIR,'matlab/VFA_T1.mnc.gz'); niak_write_vol(hdr,T1 .* mask1);
    
    % Apply spoiling correction:
    param.prot_name  = strcat('vfa_t2_', num2str(T2(3)));    % Used in output - T2 = 80ms

    % I have wrote to a mat file
    coeff = load(fullfile( b1Dir, 'SpoilingFiles', [strrep(param.prot_name,' ',''),'_ABcoeff.mat']) );

    % Note in Preibisch and Deichmann 2009, A and B are quadratic functions
    % dependent on b1map.
    
    % Get coefficients from json:
    Acoef = coeff.polyCoeffA;
    Bcoef = coeff.polyCoeffB;
    
    % Make correction factors:
    A = Acoef(1)*b1.^2 + Acoef(2)*b1 + Acoef(3);
    B = Bcoef(1)*b1.^2 + Bcoef(2)*b1 + Bcoef(3);
    
    % Apply to the T1 map - in milliseconds as simulations are done in ms
    T1corr = A + B.*T1;
    T1corr = T1corr.*mask1;
    T1corr = limitHandler(T1corr, 0, 6000);

    figure; imshow3Dfull(T1corr.* mask1, [300 2500],jet)

    hdr.file_name = strcat(DATADIR,'matlab/VFA_T1_hMRIspoil_t2_', num2str(T2(3)),'.mnc.gz'); 
    niak_write_vol(hdr,T1corr);    
     
    
    %% Recalculate M0 using this and the flash equation
    % do both images then average: Rewrite variables here for fun
    low_flip_angle = 6;    % flip angle in degrees -> USER DEFINED
    high_flip_angle = 20;  % flip angle in degrees -> USER DEFINED
    TR1 = 27;               % low flip angle repetition time of the GRE kernel in milliseconds -> USER DEFINED
    TR2 = 15;               % high flip angle repetition time of the GRE kernel in milliseconds -> USER DEFINED

    flip_a = (low_flip_angle*b1) * pi / 180; % correct for B1 and convert to radians
    x = cos(flip_a) ;
    y = exp(-TR1./T1corr);

    % Solve for M0 using equation for flass image.
    M0_1 = lfa_ur.*(1-x.*y)./ ( (1-y) .*sin(flip_a));

    % second image
    flip_a = (high_flip_angle*b1) * pi / 180; % correct for B1 and convert to radians
    x = cos(flip_a) ;
    y = exp(-TR2./T1corr);

    % Solve for M0 using equation for flass image.
    M0_2 = hfa_ur.*(1-x.*y)./ ( (1-y) .*sin(flip_a));

    M0_spoil = 2.5* (M0_1 + M0_2)./2 .* mask1; % 2.5 is the gain factor. 
    M0_spoil = limitHandler(M0_spoil, 0, 20000);

    figure; imshow3Dfull(M0_spoil , [00 15000])   
     
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% MP2RAGE Section
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   
    % Start with no B1 map for comparison, then make maps with B1
    % correction:
    
    %% MP2RAGE approach using the T1 maps from scanner
    % this doesnt make sense, you have to run the scripts to get M0 anyway.
    % Just to please the crowd, ill run the MP2RAGE script with B1 = 1;

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
    B1.img = ones(size(b1));
    brain.img = mask1;
    
    tic
    [ ~, ~, m0_1] = CR_T1B1correctpackageTFL_withM0( B1, MP2RAGEimg, MP2RAGEINV2img, MP2RAGE, brain, 0.96);
    toc
      
    m0_1 = m0_1.img;
    m0_1 = double(limitHandler(m0_1, 0, 20000));

    % figure; imshow3Dfull(m0_1 , [00 15000])
    
    hdr.file_name = strcat(DATADIR,'matlab/MP2RAGE_M0_noB1.mnc.gz'); niak_write_vol(hdr, m0_1);


    %% T1 and M0 from the SPARSE MP2RAGE:


    MP2RAGE.TRFLASH = 6.4e-3; % TR of the GRE readout
    MP2RAGE.TIs = [0.94 2.83];% inversion times - time between middle of refocusing pulse and excitatoin of the k-space center encoding
    MP2RAGE.NZslices = [ceil(175/2) floor(175/2)];%  should be two values, number of excitations before k-space center, and number after. [Slices Per Slab * [PartialFourierInSlice-0.5  0.5] ]
    MP2RAGE.FlipDegrees = [4 5];% Flip angle of the two readouts in degrees


    MP2RAGEimg.img = sp_mp2r_uni; % load_untouch_nii(MP2RAGE.filenameUNI);
    MP2RAGEINV2img.img = sp_mp2r_inv2; % load_untouch_nii(MP2RAGE.filenameINV2);
    B1.img = ones(size(b1));
    brain.img = mask1;

    tic
    [ ~, ~, m0_2] = CR_T1B1correctpackageTFL_withM0( B1, MP2RAGEimg, MP2RAGEINV2img, MP2RAGE, brain, 0.96);
    toc

    m0_2 = m0_2.img;
    m0_2 = double(limitHandler(m0_2, 0, 20000));
    
    % figure; imshow3Dfull(m0_2 , [00 15000])
    hdr.file_name = strcat(DATADIR,'matlab/spMP2RAGE_M0_noB1.mnc.gz'); niak_write_vol(hdr, m0_2);
    

    mp2_sat_2k_raw   = calcMTsatThruLookupTablewithDummyV3( mt2k, ones(size(mt2k)), t1_1, mask1, m0_1/2.5, 0, 1, 27, 6, 0);
    spmp2_sat_2k_raw = calcMTsatThruLookupTablewithDummyV3( mt2k, ones(size(mt2k)), t1_2, mask1, m0_2/2.5,0, 1, 27, 6,0);

    figure; imshow3Dfull(mp2_sat_2k_raw , [0 0.03], turbo)
    figure; imshow3Dfull(spmp2_sat_2k_raw , [0 0.03], turbo)


    hdr.file_name = strcat(DATADIR,'matlab/MP2RAGE_MTsat_noB1.mnc.gz'); niak_write_vol(hdr, mp2_sat_2k_raw);
    hdr.file_name = strcat(DATADIR,'matlab/spMP2RAGE_MTsat_noB1.mnc.gz'); niak_write_vol(hdr, spmp2_sat_2k_raw);
    
    % Save the T1s as well just to be certain they are all being
    % consistently processed:
     hdr.file_name = strcat(DATADIR,'matlab/MP2RAGE_T1_noB1.mnc.gz'); niak_write_vol(hdr, t1_1.* mask1);
    hdr.file_name = strcat(DATADIR,'matlab/spMP2RAGE_T1_noB1.mnc.gz'); niak_write_vol(hdr, t1_2 .* mask1);

     %% Now from the MP2RAGE:
    
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
    [ T1map, MP2RAGEcorrected, Appmap2] = CR_T1B1correctpackageTFL_withM0( B1, MP2RAGEimg, MP2RAGEINV2img, MP2RAGE, brain, 0.96);
    toc
    
    T1_map = T1map.img;
    App_mp2 = Appmap2.img;
    
    T1_map = limitHandler(T1_map, 0 , 6000);
    App_mp2 = double(limitHandler(App_mp2, 0, 20000));
    
    % figure; imshow3Dfull(T1_map, [300 2500],jet)
    % figure; imshow3Dfull(t1_maps(:,:,:,1), [300 2500],jet) % close!
    % figure; imshow3Dfull(App_mp2 , [00 15000])
    
    
    hdr.file_name = strcat(DATADIR,'matlab/MP2RAGE_T1.mnc.gz'); niak_write_vol(hdr, T1_map);


    %% T1 and M0 from the SPARSE MP2RAGE:

    MP2RAGE.B0 = 3;           % in Tesla
    MP2RAGE.TR = 5;           % MP2RAGE TR in seconds
    MP2RAGE.TRFLASH = 6.4e-3; % TR of the GRE readout
    MP2RAGE.TIs = [0.94 2.83];% inversion times - time between middle of refocusing pulse and excitatoin of the k-space center encoding
    MP2RAGE.NZslices = [ceil(175/2) floor(175/2)];%  should be two values, number of excitations before k-space center, and number after. [Slices Per Slab * [PartialFourierInSlice-0.5  0.5] ]
    MP2RAGE.FlipDegrees = [4 5];% Flip angle of the two readouts in degrees


    MP2RAGEimg.img = sp_mp2r_uni; % load_untouch_nii(MP2RAGE.filenameUNI);
    MP2RAGEINV2img.img = sp_mp2r_inv2; % load_untouch_nii(MP2RAGE.filenameINV2);
    B1.img = b1;
    brain.img = mask1;

    tic
    [ spT1map, spMP2RAGEcorrected, spAppmap2] = CR_T1B1correctpackageTFL_withM0( B1, MP2RAGEimg, MP2RAGEINV2img, MP2RAGE, brain, 0.96);
    toc

    T1_ms = spT1map.img;
    spApp_mp2 = spAppmap2.img;

    T1_ms = limitHandler(T1_ms, 0 , 6000);
    spApp_mp2 = double(limitHandler(spApp_mp2, 0, 20000));
    
%     figure; imshow3Dfull(T1_ms, [300 2500],jet)
%     figure; imshow3Dfull(spApp_mp2 , [00 15000])

    %mkdir(strcat(DATADIR,'matlab') )
    hdr.file_name = strcat(DATADIR,'matlab/sparseMP2RAGE_T1.mnc.gz'); niak_write_vol(hdr, T1_ms);
    
    
    %% Confirm all the App/M0 images are scaled the same then save.
     % Scaling is important for the MTsat calculation
    figure; imshow3Dfull(App_tfl .* mask1 , [00 15000])
    figure; imshow3Dfull(App_mp2 , [00 15000])
    figure; imshow3Dfull(spApp_mp2 , [00 15000])
    figure; imshow3Dfull(M0_spoil , [00 15000])


    hdr.file_name = strcat(DATADIR,'matlab/sparseMP2RAGE_M0.mnc.gz'); niak_write_vol(hdr, spApp_mp2);
    hdr.file_name = strcat(DATADIR,'matlab/MP2RAGE_M0.mnc.gz'); niak_write_vol(hdr, App_mp2);
    hdr.file_name = strcat(DATADIR,'matlab/VFA_M0.mnc.gz'); niak_write_vol(hdr,App_tfl .* mask1);
    hdr.file_name = strcat(DATADIR,'matlab/VFA_M0_spoil.mnc.gz'); niak_write_vol(hdr,M0_spoil .* mask1);
    
    
    
    %% Calculate MTsat
    TR = 27; %milliseconds
    flipA = 6; % flip angle in degrees


    % calculate maps

    sat_2k             = calcMTsatThruLookupTablewithDummyV3( mt2k_ur, b1, T1,    mask1,  App, 0, 1, TR, flipA,0);
    sat_2k_spoil    = calcMTsatThruLookupTablewithDummyV3( mt2k_ur, b1, T1corr, mask1, M0_spoil/2.5, 0, 1, TR, flipA,0); % M0 was scaled to match TFL
    mp2_sat_2k     = calcMTsatThruLookupTablewithDummyV3( mt2k_ur, b1, T1_map, mask1, App_mp2/2.5, 0, 1, TR, flipA, 0);
    spmp2_sat_2k = calcMTsatThruLookupTablewithDummyV3( mt2k_ur, b1, T1_ms, mask1, spApp_mp2/2.5,0, 1, TR, flipA,0);

    % view to make sure all the parameters were set correctly
    figure; imshow3Dfull(sat_2k , [0 0.03], turbo)
    figure; imshow3Dfull(sat_2k_spoil , [0 0.03], turbo)
    figure; imshow3Dfull(mp2_sat_2k , [0 0.03], turbo)
    figure; imshow3Dfull(spmp2_sat_2k , [0 0.03], turbo)
    
    
       %% Correct for B1 (code to do so is below)
    
    R1_s = (1./T1)*1000; % need to convert to 1/s from 1/ms
    R1_spoil = (1./T1corr)*1000; % need to convert to 1/s from 1/ms
    R1_mp2_s = (1./T1_map)*1000;
    R1_spmp2_s = (1./T1_ms)*1000;

    corr_MTsat_2k       = MTsat_B1corr_factor_map(b1, R1_s,     3.586,fitValues_S_1);
    corr_MTsat_2k_spoil = MTsat_B1corr_factor_map(b1, R1_spoil, 3.586,fitValues_S_1);
    mp2_corr_MTsat_2k   = MTsat_B1corr_factor_map(b1, R1_mp2_s, 3.586,fitValues_S_1);
    spmp2_corr_MTsat_2k = MTsat_B1corr_factor_map(b1, R1_spmp2_s, 3.586,fitValues_S_1);

    sat_2k_c       = sat_2k.*(1 + corr_MTsat_2k) .* mask;
    sat_2k_spoil_c = sat_2k_spoil.*( 1 + corr_MTsat_2k_spoil) .* mask;
    mp2_sat_2k_c   = mp2_sat_2k.*(1 + mp2_corr_MTsat_2k) .* mask;
    spmp2_sat_2k_c = spmp2_sat_2k .*(1 + spmp2_corr_MTsat_2k) .* mask;


    %% Other things, save if you want
    sat_2k_c       = limitHandler(sat_2k_c, 0 , 0.3);
    sat_2k_spoil_c = limitHandler(sat_2k_spoil_c, 0 , 0.3);
    mp2_sat_2k_c   = limitHandler(mp2_sat_2k_c, 0 , 0.3);
    spmp2_sat_2k_c = limitHandler(spmp2_sat_2k_c, 0 , 0.3);

    figure; imshow3Dfull(sat_2k_c , [0 0.03], turbo)
    figure; imshow3Dfull(sat_2k_spoil_c , [0 0.03], turbo)
    figure; imshow3Dfull(mp2_sat_2k_c , [0 0.03], turbo)
    figure; imshow3Dfull(spmp2_sat_2k_c , [0 0.03], turbo)


    hdr.file_name = strcat(DATADIR,'matlab/VFA_MTsat_2k.mnc.gz'); niak_write_vol(hdr, sat_2k_c);
    hdr.file_name = strcat(DATADIR,'matlab/VFA_MTsat_2k_spoil.mnc.gz'); niak_write_vol(hdr, sat_2k_spoil_c);
    hdr.file_name = strcat(DATADIR,'matlab/MP2RAGE_MTsat_2k.mnc.gz'); niak_write_vol(hdr, mp2_sat_2k_c);
    hdr.file_name = strcat(DATADIR,'matlab/spMP2RAGE_MTsat_2k.mnc.gz'); niak_write_vol(hdr, spmp2_sat_2k_c);
      
    
    % clear variables for next loop
    clearvars -except subjectNames z b1Dir
    z = z+1;
    
    
    
end
