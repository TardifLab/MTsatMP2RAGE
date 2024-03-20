%% 
% THe goal of this script is to estimate how the SNR would change in the
% VFA T1 experiment if the scan time was matched.

% Number of lines in VFA. 6min scan. TR = 27 ms. 6 mins = 360,000 ms
% 360,000/27 = 13 300 (there is a couple seconds of 'prescan' there too, so rounded down)

MP2RAGE_time = 427; % seconds
csMP2RAGE_time = 258; % seconds
VFA_time = 558; % seconds

numLines = 13300;


% %% The small angle approximation holds well for flip angles <25
% a1 = deg2rad(linspace(0, 35, 100));
% ttt = 2*tan(a1/2);
% plot(rad2deg(a1), rad2deg(ttt))
% hold on
% refline(1,0)

%% Need to come up with VFA TR's that would fit
%% MP2RAGE Match
tr1 = round(27*MP2RAGE_time/VFA_time);
tr2 = round(15*MP2RAGE_time/VFA_time);

tTotal = (tr1*numLines + tr2*numLines)/1000;

% update flip angles:
T1 = 1000; % milliseconds

% As is done by Helms et al, 2011 MRM Signal Bias in VFA.
ernst1 = acos(exp(-tr1/T1));
ernst2 = acos(exp(-tr2/T1));

a1= round(rad2deg(ernst1 * 0.4142));
a2= round(rad2deg(ernst2 * 2.4142));


%% csMP2RAGE Match
cstr1 = round(27*csMP2RAGE_time/VFA_time);
cstr2 = round(15*csMP2RAGE_time/VFA_time);

cstTotal = (cstr1*numLines + cstr2*numLines)/1000;

% As is done by Helms et al, 2011 MRM Signal Bias in VFA.
csernst1 = acos(exp(-cstr1/T1));
csernst2 = acos(exp(-cstr2/T1));

csa1= round(rad2deg(csernst1 * 0.4142));
csa2= round(rad2deg(csernst2 * 2.4142));

%% From this, we want to simulate the original , and the two modified protocols to see how the SNR


% see the function in Batch_sim/functions/CR_add_ihMTsat_SNR.m
Params.B0 = 3;
Params.MTC = 0; % Magnetization Transfer Contrast
Params.TissueType = 'WM';
Params.numExcitation = 1;
Params.echoSpacing = 7.7/1000;

Params = DefaultCortexTissueParams(Params);
Params = CalcImagingParams(Params);
Params.PerfectSpoiling = 1;
Params = CalcVariableImagingParams(Params);


%% Simulate VFA - acquired protocol

% two images with differing flip angles:
vfa_sig = zeros( 3, 2);

% Image 1 - PDw
Params.TR = 27/1000;
Params.flipAngle = 6;
[vfa_sig(1,1), ~, ~] = BlochSimFlashSequence_v2( Params );


% Image 2 - T1w
Params.TR = 15/1000;
Params.flipAngle = 20;

[vfa_sig(1,2), ~, ~] = BlochSimFlashSequence_v2( Params );


%% MP2RAGE Matched

% Image 1 - PDw
Params.TR = tr1/1000;
Params.flipAngle = a1;
[vfa_sig(2,1), ~, ~] = BlochSimFlashSequence_v2( Params );


% Image 2 - T1w
Params.TR = tr2/1000;
Params.flipAngle = a2;

[vfa_sig(2,2), ~, ~] = BlochSimFlashSequence_v2( Params );


%% MP2RAGE Matched

% Image 1 - PDw
Params.TR = cstr1/1000;
Params.flipAngle = csa1;
[vfa_sig(3,1), ~, ~] = BlochSimFlashSequence_v2( Params );


% Image 2 - T1w
Params.TR = cstr2/1000;
Params.flipAngle = csa2;

[vfa_sig(3,2), ~, ~] = BlochSimFlashSequence_v2( Params );


% Confirm consistent T1 (so calculations above worked)
vfa_T1 = zeros( 3, 1);
vfa_T1(1) = Helms_VFA_T1_2008( vfa_sig(1,1), vfa_sig(1,2), 6, 20, 27, 15, 1);
vfa_T1(2) = Helms_VFA_T1_2008( vfa_sig(2,1), vfa_sig(2,2), a1, a2, tr1, tr2, 1);
vfa_T1(3) = Helms_VFA_T1_2008( vfa_sig(3,1), vfa_sig(3,2), csa1, csa2, cstr1, cstr2, 1);

% Excellent!

%% Now calculate SNR:
SNR = 70;
noiseLvl = vfa_sig(1,1)/SNR; % This reference was used originally.

additive_noise = normrnd( 1,noiseLvl ,30000,6) - 1;

noisy_vfa_T1_1 = Helms_VFA_T1_2008( additive_noise(:,1) + vfa_sig(1,1), ...
    additive_noise(:,2) + vfa_sig(1,2), 6, 20, 27, 15, 1);

noisy_vfa_T1_2 = Helms_VFA_T1_2008( additive_noise(:,3) +vfa_sig(2,1), ...
    additive_noise(:,4) + vfa_sig(2,2), a1, a2, tr1, tr2, 1);

noisy_vfa_T1_3 = Helms_VFA_T1_2008( additive_noise(:,5) +vfa_sig(3,1), ...
    additive_noise(:,6) +vfa_sig(3,2), csa1, csa2, cstr1, cstr2, 1);


% Take the noise level propagated to 
noisy_vfa_T1_1 = std(noisy_vfa_T1_1);
noisy_vfa_T1_2 = std(noisy_vfa_T1_2);
noisy_vfa_T1_3 = std(noisy_vfa_T1_3);


%% Difference relative:
d1 = (noisy_vfa_T1_2-noisy_vfa_T1_1)/noisy_vfa_T1_1 *100;
d2 = (noisy_vfa_T1_3-noisy_vfa_T1_1)/noisy_vfa_T1_1 *100;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% this suggests that the noise level of the MP2RAGE time-matched VFA T1 map
% would be 13.5% more than what was acquired

% this suggests that the noise level of the cdMP2RAGE time-matched VFA T1 map
% would be 44.8% more than what was acquired

