% Look into the impact of MT effects on T1 mapping
% compare IR, MP2RAGE and VFA
addpath(genpath('E:\GitHub\NeuroImagingMatlab\QuantitativeFitting'))
addpath(genpath('E:\GitHub\NeuroImagingMatlab\colourmaps\OptimalMaps'))

% add necessary paths:
setupSimPaths_MTsatMP2RAGE;

clear all;
clc;

outputPath = 'E:\GitHub\MTsatMP2RAGE\MTeffect_in_T1mapping\savedOutputs\matrices';
outputPathFig = 'E:\GitHub\MTsatMP2RAGE\MTeffect_in_T1mapping\savedOutputs';
savePrefix ='R1B_4s'; % 'Sled2001', 'LorentzT270',

% Notes on savePrefix
% ModulateTP -> not realistic for all values to be maintained, so we use
% linear models to vary them.

%%

% Lets fix all the tissue parameters to be for the WM:
Params.B0 = 3;
Params.MTC = 0; % Magnetization Transfer Contrast
Params.TissueType = 'WM';
Params.numExcitation = 1;
Params.echoSpacing = 7.7/1000;

Params = DefaultCortexTissueParams(Params);
Params = CalcImagingParams(Params);
Params = CalcVariableImagingParams(Params);

Params.PerfectSpoiling = 1;
Params.R1b = 4; % The assumed value for Sled and Pike
Params.Ra = 0.4; % The assumed value for Sled and Pike
Params.R = []; % clear it to avoid errors

%% To test the MT effect, we will simulate for a range of bound pool parameters
% Compared to the main script, here we will vary the values based on Sled
% and Pike 2001. 


T1 = linspace(600, 2500, 50);
% Then vectorize for easy looping:
T1_m = T1(:)./1000;

% We actually input Ra_obs in, so invert:
Raobs_m = 1./T1_m;
%[M0b, kf, Ra, T2a, T2b, kr] = LinearModelVals_SledandPike2001( Raobs_m );
[M0b, kf, Ra, T2a, T2b, kr] = qMT_LinearModelVals( Raobs_m);
simNum = length(Raobs_m);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% We want to see how much the observed T1 changes as a function of the
% bound pool fraction used.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Simulate VFA
% two images with differing flip angles:
vfa_sig = zeros( simNum, 2);

% Image 1 - PDw
Params.TR = 27/1000;
Params.flipAngle = 6;

tic % About 18 seconds for 450 options

for i = 1 : simNum

    [vfa_sig(i,1), ~, ~] = BlochSimFlashSequence_v2( Params,...
        'M0b', M0b(i), 'T2a', T2a(i), 'T2b', T2b(i),'kf', kf(i), 'kr', kr(i));
end

% Image 2 - T1w
Params.TR = 15/1000;
Params.flipAngle = 20;

for i = 1: simNum
        [vfa_sig(i,2), ~, ~] = BlochSimFlashSequence_v2( Params,...
            'M0b', M0b(i),'T2a', T2a(i), 'T2b', T2b(i),'kf', kf(i), 'kr', kr(i) );
end
toc

vfa_T1 = Helms_VFA_T1_2008( vfa_sig(:,1), vfa_sig(:,2), 6, 20, 27, 15, 1);


save( fullfile( outputPath, [savePrefix,'_VFA_sim_rev2.mat']), "vfa_sig");
save( fullfile( outputPath, [savePrefix,'_vfa_T1_rev2.mat']), "vfa_T1");

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Simulate Inversion Recovery Signal.
clear Params
%% New IR Simulations
Params.B0 = 3;
Params.TissueType = 'WM';
Params = DefaultCortexTissueParams(Params);
Params.PerfectSpoiling = 1;
Params.R1b = 4; % The assumed value for Sled and Pike
Params.Ra = 0.4; % The assumed value for Sled and Pike

Params.InvPulseDur = 10.24/1000;
Params.ExcPulseDur = 3.072/1000;
Params.RefocusPulseDur = 3/1000;
Params.TI = 50/1000;
Params.TR = 2;
Params.TE = 8.5/1000;

Params.Inv.PulseOpt.beta = 672;  
Params.Inv.PulseOpt.n = 1; 
Params.Inv.PulseOpt.mu = 5;
Params.Inv.PulseOpt.A0 = 13.726;
Params.Inv.nSamples = 512;
Params.Inv.Trf = 10.24/1000;
Params.Inv.shape = 'hsn';


Params.Exc.alpha = 90; % degrees
Params.Exc.Trf = 3.072/1000;
Params.Exc.PulseOpt.TBW = 4;  
Params.Exc.nSamples = 512;
Params.Exc.shape = 'sinc';

Params.refocus.alpha = 180; % degrees
Params.refocus.Trf = 3.000/1000;
Params.refocus.nSamples = 512;
Params.refocus.shape = 'gaussian';
Params.refocus.PulseOpt = [];

% TI = [30; 250; 500; 750; 1000; 1500]/1000; % From first submission - matching NIST acquisition
TI = [10:5:60, 85, 125, 200, 350, 500, 750, 1000, 1200, 1400, 1600]'/1000;
numTI = length(TI);
IR_sig = zeros( simNum, numTI);

tic % Roughly 1 min for sim and fitting
for i = 1: simNum
    for j = 1:numTI
        
        [IR_sig(i,j), ~, ~] = BlochSim_SpinEcho_IR_Sequence( Params,...
            'M0b', M0b(i), 'TI', TI(j),...
            'T2a', T2a(i), 'T2b', T2b(i),'kf', kf(i), ...
            'kr', kr(i) );

    end
end
toc

tic
[IR_T1_mono, ~] = fit_T1_IR_data(IR_sig, TI,'complex');
toc

IR_T1_mono = IR_T1_mono *1000;

% biexponential - keep only some TIs

tic
[IR_T1_biL, IR_T1_biS, ~] = fit_biexp_T1_IR_data(IR_sig, TI);
toc

IR_T1_biL = IR_T1_biL *1000;
IR_T1_biS = IR_T1_biS *1000;


save( fullfile(outputPath,[savePrefix,'_IR_sim_rev2.mat']), "IR_sig");
save( fullfile(outputPath,[savePrefix,'_IR_T1_mono_rev2.mat']), "IR_T1_mono");
save( fullfile(outputPath,[savePrefix,'_IR_T1_biLong_rev2.mat']), "IR_T1_biL");
save( fullfile(outputPath,[savePrefix,'_IR_T1_biShort_rev2.mat']), "IR_T1_biS");



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Simulate MP2RAGE
clear Params

Params.B0 = 3;
Params.TissueType = 'WM';
Params = DefaultCortexTissueParams(Params);
Params.PerfectSpoiling = 1;
Params.WExcDur = 0.1/1000;
Params.PerfectSpoiling = 1;
Params.R1b = 4; % The assumed value for Sled and Pike
Params.Ra = 0.4; % The assumed value for Sled and Pike

% New Params
Params.Inv.PulseOpt.beta = 672;  
Params.Inv.PulseOpt.n = 1; 
Params.Inv.PulseOpt.mu = 5;
Params.Inv.PulseOpt.A0 = 17.5;
Params.Inv.nSamples = 512;
Params.Inv.Trf = 10.24/1000;
Params.Inv.shape = 'hsn';
Params.simAdiabatic = 1;
Params.Exc.shape = 'hard';
Params.Exc.Trf = 0.1/1000;
Params.Exc.nSamples = 20;

% Old
Params.TR = 5000/1000;
Params.flipAngle = [4,5];
Params.numExcitation = 175;
Params.echoSpacing = 7.7/1000;
Params.Readout = 'linear';
Params.TI = [940, 2830]./1000;
Params.DummyEcho = 0;
Params.PerfectSpoiling = 1;
Params.RFspoiling = 0;

MP2_sig = zeros( simNum, 2);

tic % roughly 3 mins
for i = 1: simNum
    [MP2_sig(i,1),MP2_sig(i,2), ~, ~] = BlochSim_MP2RAGESequence_6vec( Params,...
        'M0b', M0b(i), 'T2a', T2a(i), 'T2b', T2b(i),'kf', kf(i), 'kr', kr(i) );
end
toc

MP2RAGE.B0          = Params.B0;                  % In Tesla
MP2RAGE.TR          = Params.TR;                  % MP2RAGE TR in seconds
MP2RAGE.TRFLASH     = Params.echoSpacing;             % TR of the GRE readout
MP2RAGE.TIs         = Params.TI;   % Inversion times - time between middle of refocusing pulse and excitatoin of the k-space center encoding
MP2RAGE.NZslices    = Params.numExcitation;            % Slices Per Slab * [PartialFourierInSlice-0.5  0.5]
MP2RAGE.FlipDegrees = Params.flipAngle;              % Flip angle of the two readouts in degrees

[mp2rage_T1, PD, ~] = MP2RAGE_dictionaryMatching(MP2RAGE, MP2_sig(:,1), MP2_sig(:,2),...
    ones(size(MP2_sig(:,2))), [0.0005, 0.005], 0);
mp2rage_T1 = mp2rage_T1*1000;

save( fullfile(outputPath,[savePrefix,'_MP2_sim_rev2.mat']),"MP2_sig");
save( fullfile(outputPath,[savePrefix,'_mp2rage_T1_rev2.mat']),"mp2rage_T1");

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Make Plots!
% We are interested how the measured T1 varies as a function of the bound
% pool fraction. So plot x-axis as true T1, y-axis = measured T1, and then
% do separate lines for each bound pool fraction.


% load( fullfile( outputPath, [savePrefix,'_vfa_T1_rev2.mat']));
% load( fullfile(outputPath,[savePrefix,'_IR_T1_mono_rev2.mat']));
% load( fullfile(outputPath,[savePrefix,'_IR_T1_biLong_rev2.mat']));
% load( fullfile(outputPath,[savePrefix,'_IR_T1_biShort_rev2.mat']));
% load( fullfile(outputPath,[savePrefix,'_mp2rage_T1_rev2.mat']));


%% Set up values:
lLim = 600;
hLim = 1800;
ticV = lLim:400:hLim;

% x data = T1:
T1 = linspace(600, 2500, 50);

x = T1(:);
refy = x;
cm_data = 	[0, 0.4470, 0.9410; ...
    0.4660, 0.7740, 0.1880;...
    0.9290, 0.6940, 0.1250;...
    0.9350, 0.0780, 0.1840];



%% Make a second plot that shows the difference relative to T1 long. 

y = [(vfa_T1-IR_T1_biL)./IR_T1_biL,...
    (IR_T1_mono-IR_T1_biL)./IR_T1_biL,...,
    (mp2rage_T1-IR_T1_biL)./IR_T1_biL]*100;

FontSize = 20;
ylLim = -40;
yhLim = 40;
yticV = ylLim:20:yhLim;

figure; 
hold on;
for i = 1:3
    plot(x,y(:,i),'Color', cm_data(i,: ),'LineWidth',2)
end
xlabel( 'Input T_{1,obs} (ms)', 'FontSize', FontSize, 'FontWeight', 'bold' )
ylabel( '% Diff. to T_{1,Long}' , 'FontSize', FontSize, 'FontWeight', 'bold');
xlim([lLim, hLim]); ylim([ylLim, yhLim])
xticks(ticV); yticks(yticV)
%title('VFA', 'FontSize', FontSize, 'FontWeight', 'bold')
ax = gca;  ax.FontSize = FontSize; 
hold on
yline(0,':','Color',[0,0,0],'LineWidth',1.5)
legend('VFA','IR Mono','MP2RAGE', 'FontSize', FontSize-4, 'location', 'best')


saveas(gcf,fullfile(outputPathFig,[savePrefix,'_Comp2T1long.png']));












% %% VFA, IR and MP2RAGE
% 
% y = [vfa_T1, IR_T1_mono, IR_T1_biL, mp2rage_T1];
% 
% figure; 
% hold on;
% for i = 1:4
%     plot(x,y(:,i),'Color', cm_data(i,: ),'LineWidth',2)
% end
% xlabel( 'Input T_{1,obs} (ms)', 'FontSize', FontSize, 'FontWeight', 'bold' )
% ylabel( 'Simulated T_{1,obs} (ms)' , 'FontSize', FontSize, 'FontWeight', 'bold');
% xlim([lLim, hLim]); ylim([lLim, hLim])
% xticks(ticV); yticks(ticV)
% %title('VFA', 'FontSize', FontSize, 'FontWeight', 'bold')
% ax = gca;  ax.FontSize = FontSize; 
% hold on
% plot(x,refy,':','Color',[0,0,0],'LineWidth',1.5)
% 
% % Add lines for typical GM and WM T1s
% point = [850, 850];
% axLims = [lLim hLim lLim hLim];  %[x-min, x-max, y-min, y-max] axis limits
% plot([point(1), point(1)], [axLims(3), point(2)], 'k-')  %vertical line
% plot([axLims(1), point(1)], [point(2), point(2)], 'k-')  %horizontal line
% text(lLim + 20, point(1) + 70, strcat(num2str(point(1)), " ms"),'FontSize', FontSize-2)
% 
% point = [1400, 1400];
% axLims = [lLim hLim lLim hLim];  %[x-min, x-max, y-min, y-max] axis limits
% plot([point(1), point(1)], [axLims(3), point(2)], 'k-')  %vertical line
% plot([axLims(1), point(1)], [point(2), point(2)], 'k-')  %horizontal line
% text(lLim + 20, point(1) + 70, strcat(num2str(point(1)), " ms"),'FontSize', FontSize-2)
% 
% legend('VFA','IR Mono','IR Bi-Long','MP2RAGE', 'FontSize', FontSize, 'location', 'best')
% hold off
% 
% saveas(gcf,fullfile(outputPathFig,[savePrefix,'_Wang_T1plot_comparison.png']));
% 
% temp = mean(y(:,3)-y(:,1)); %Average difference was 70 ms between the two. 

















