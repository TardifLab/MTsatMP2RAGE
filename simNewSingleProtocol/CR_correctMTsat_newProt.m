%% Correct MTsat maps from 3 protocols after running:
%   simSeq_M0B_R1obs_3prot.m  and...
%   CR_R1vsM0B_correlation.m

% We will assume that you have loaded your images for a single subject, and
% have the following [names as contained in square brackets]:
% - [T1map] (B1+ corrected)
% - [MTsat]  image (calculated using B1+ corrected T1 and S0, and corrected
%               for excitation flip angles in the equation)
% - [B1]+ map
% - brain [mask] - reasonable tight, remove skull and other tissue -
%                   optional

%% This requires the dependencies from: https://github.com/TardifLab/MTsatB1correction
OutputDir = 'Directory\b1Correction\outputs';

% load in the fit results from simulations
fitValues= load(fullfile(OutputDir,'fitValues_',seqStr,'.mat'));
fitValues = fitValues.fitValues;

%% You can choose whether you want to do additional processing on the B1 map
b1 = double(b1);
b1 = limitHandler(b1, 0.5,1.4);
b1 = CR_imgaussfilt3_withMask(b1, mask1, 5); %light smoothing to the map

figure; imshow3Dfull(b1, [0.6 1.2],jet) % check result. Values should be centered around 1, not 100.

%% MTsat correction
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


































            