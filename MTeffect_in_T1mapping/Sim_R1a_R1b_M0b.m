%% We don't know R1b, R1a or M0b. Simulate over them all!

setupSimPaths_MTsatMP2RAGE;

clear all;
clc;

outputPath = 'E:\GitHub\MTsatMP2RAGE\MTeffect_in_T1mapping\savedOutputs\matrices';
outputPathFig = 'E:\GitHub\MTsatMP2RAGE\MTeffect_in_T1mapping\savedOutputs';
savePrefix ='SimAll'; 

%% To test the MT effect, we will simulate for a range of bound pool parameters

M0b = 0:0.03:0.30;
R1a = linspace(0.2, 2, 10);
R1b = linspace(0.5, 4, 10);


[M0b_m, R1a_m, R1b_m] = ndgrid(M0b, R1a, R1b);

% Then vectorize for easy looping:
M0b_m = M0b_m(:); R1a_m = R1a_m(:); R1b_m = R1b_m(:);

simNum = length(R1a_m);


%% Simulate VFA
% Lets fix all the tissue parameters to be for the WM:
Params.B0 = 3;
Params.MTC = 0; % Magnetization Transfer Contrast
Params.TissueType = 'WM';
Params.numExcitation = 1;
Params.echoSpacing = 7.7/1000;
Params.WExcDur = 0.1/1000;
Params.kf = 4.45;
Params.kr = 28.34;
Params = DefaultCortexTissueParams(Params);
Params.PerfectSpoiling = 1;
Params = CalcImagingParams(Params);
Params = CalcVariableImagingParams(Params);

% two images with differing flip angles:
vfa_sig = zeros( simNum, 2);

% Image 1 - PDw
Params.TR = 27/1000;
Params.flipAngle = 6;

tic % About 18 seconds for 450 options

for i = 1 : simNum
    [vfa_sig(i,1), ~, ~] = BlochSimFlashSequence_v2( Params,'M0b', M0b_m(i), 'Ra',R1a_m(i),'R1b', R1b_m(i) );
end

% Image 2 - T1w
Params.TR = 15/1000;
Params.flipAngle = 20;

for i = 1: simNum
    [vfa_sig(i,2), ~, ~] = BlochSimFlashSequence_v2( Params,'M0b', M0b_m(i), 'Ra',R1a_m(i),'R1b', R1b_m(i) );
end
toc

vfa_R1 = CR_fit_R1_Aapp_Helms2008( vfa_sig(:,1), vfa_sig(:,2), 6, 20, 27, 15, 1);
vfa_T1 = 1./vfa_R1;

save( fullfile( outputPath, [savePrefix,'_VFA_sim.mat']), "vfa_sig");
save( fullfile( outputPath, [savePrefix,'_vfa_T1.mat']), "vfa_T1");


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear Params
%% New IR Simulations
Params.B0 = 3;
Params.MTC = 0; % Magnetization Transfer Contrast
Params.TissueType = 'WM';
Params.numExcitation = 1;
Params.echoSpacing = 7.7/1000;
Params.WExcDur = 0.1/1000;
Params.kf = 4.45;
Params.kr = 28.34;
Params = DefaultCortexTissueParams(Params);
Params = CalcImagingParams(Params);
Params = CalcVariableImagingParams(Params);
Params.PerfectSpoiling = 1;

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
                'M0b', M0b_m(i), 'TI', TI(j), 'Ra',R1a_m(i),'R1b', R1b_m(i) );

    end
    disp( i/length(IR_sig) *100)
end
toc

% Fit 2 different IR values
% Monoexponential 

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

save( fullfile(outputPath,[savePrefix,'_IR_sim_rev1.mat']), "IR_sig");
save( fullfile(outputPath,[savePrefix,'_IR_T1_mono_rev1.mat']), "IR_T1_mono");
save( fullfile(outputPath,[savePrefix,'_IR_T1_biLong_rev1.mat']), "IR_T1_biL");
save( fullfile(outputPath,[savePrefix,'_IR_T1_biShort_rev1.mat']), "IR_T1_biS");




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Simulate MP2RAGE
clear Params
Params.B0 = 3;
Params.MTC = 0; % Magnetization Transfer Contrast
Params.TissueType = 'WM';
Params.numExcitation = 1;
Params.echoSpacing = 7.7/1000;
Params.WExcDur = 0.1/1000;
Params.kf = 4.45;
Params.kr = 28.34;
Params = DefaultCortexTissueParams(Params);
Params = CalcImagingParams(Params);
Params = CalcVariableImagingParams(Params);
Params.PerfectSpoiling = 1;

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
                    'M0b', M0b_m(i), 'Ra',R1a_m(i),'R1b', R1b_m(i) );
    
    disp( i/length(IR_sig) *100)
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


%% Now plot each as some sort of surface. 


vfa_mat = reshape(vfa_T1, length(M0b), length(R1a), length(R1b) );
IRm_mat = reshape(IR_T1_mono, length(M0b), length(R1a), length(R1b) );
IRL_mat = reshape(IR_T1_biL, length(M0b), length(R1a), length(R1b) );
mp2_mat = reshape(mp2rage_T1, length(M0b), length(R1a), length(R1b) );


diffVFA = (vfa_mat-IRL_mat)./IRL_mat *100;
diffIRm = (IRm_mat-IRL_mat)./IRL_mat *100;
diffMP2 = (mp2_mat-IRL_mat)./IRL_mat *100;


% Maybe fix Ra? Lets see what a good value is:
fixR1obs = 1;
fixR = Params.kr;
fixR1b = 1;
fixM0b = 0.15;
fixR1a = fixR1obs - ((fixR * fixM0b * (fixR1b - fixR1obs)) / (fixR1b - fixR1obs + fixR));

% try 0.4 and 1?



[x,y] = ndgrid( M0b(2:end), R1b);

FSz = 18;
climLow = -80;
climHigh = 80;


r1a_slice = 2;
figure; tiledlayout(1,3, "TileSpacing","tight"); 
nexttile;
surfMap = squeeze(diffVFA(2:end, r1a_slice,:));
surf( x, y, surfMap,'FaceColor','interp',...
   'EdgeColor','none')
xlabel('M_{0B}', 'FontSize', FSz, 'FontWeight', 'bold')
ylabel('R_{1B} (s^{-1})', 'FontSize', FSz, 'FontWeight', 'bold')
zlabel( '% Diff to T_{1,Long}' , 'FontSize', FSz, 'FontWeight', 'bold');
ax = gca;  ax.FontSize = FSz; 
clim([ climLow, climHigh])
colormap(bipolar)
title('VFA')
axis tight
view(0,90)

nexttile;
surfMap = squeeze(diffMP2(2:end, r1a_slice,:));
surf( x, y, surfMap,'FaceColor','interp',...
   'EdgeColor','none')
xlabel('M_{0B}', 'FontSize', FSz, 'FontWeight', 'bold')
ylabel('R_{1B} (s^{-1})', 'FontSize', FSz, 'FontWeight', 'bold')
zlabel( '% Diff to T_{1,Long}' , 'FontSize', FSz, 'FontWeight', 'bold');
ax = gca;  ax.FontSize = FSz; 
clim([ climLow, climHigh])
colormap(bipolar)
title('MP2RAGE')
axis tight
view(0,90)


nexttile;
surfMap = squeeze(diffIRm(2:end, r1a_slice,:));
surf( x, y, surfMap,'FaceColor','interp',...
   'EdgeColor','none')
xlabel('M_{0B}', 'FontSize', FSz, 'FontWeight', 'bold')
ylabel('R_{1B} (s^{-1})', 'FontSize', FSz, 'FontWeight', 'bold')
zlabel( '% Diff to T_{1,Long}' , 'FontSize', FSz, 'FontWeight', 'bold');
ax = gca;  ax.FontSize = FSz; 
clim([ climLow, climHigh])
colormap(bipolar)
cbh = colorbar;
set(cbh,'YTick',linspace(climLow,climHigh, 5))
title('Mono IR')
axis tight
view(0,90)

set(gcf,'Position',[100 100 1100 400])
sgtitle( 'R_{1A} = 0.4 s^{-1}' , 'FontSize', FSz+4, 'FontWeight', 'bold')
saveas(gcf,fullfile(outputPathFig,[savePrefix,'_R1A_0p4.png']));

%% R1a = 5

r1a_slice = 5;
figure; tiledlayout(1,3, "TileSpacing","tight"); 
nexttile;
surfMap = squeeze(diffVFA(2:end, r1a_slice,:));
surf( x, y, surfMap,'FaceColor','interp',...
   'EdgeColor','none')
xlabel('M_{0B}', 'FontSize', FSz, 'FontWeight', 'bold')
ylabel('R_{1B} (s^{-1})', 'FontSize', FSz, 'FontWeight', 'bold')
zlabel( '% Diff to T_{1,Long}' , 'FontSize', FSz, 'FontWeight', 'bold');
ax = gca;  ax.FontSize = FSz; 
clim([ climLow, climHigh])
colormap(bipolar)
title('VFA')
axis tight
view(0,90)

nexttile;
surfMap = squeeze(diffMP2(2:end, r1a_slice,:));
surf( x, y, surfMap,'FaceColor','interp',...
   'EdgeColor','none')
xlabel('M_{0B}', 'FontSize', FSz, 'FontWeight', 'bold')
ylabel('R_{1B} (s^{-1})', 'FontSize', FSz, 'FontWeight', 'bold')
zlabel( '% Diff to T_{1,Long}' , 'FontSize', FSz, 'FontWeight', 'bold');
ax = gca;  ax.FontSize = FSz; 
clim([ climLow, climHigh])
colormap(bipolar)
title('MP2RAGE')
axis tight
view(0,90)


nexttile;
surfMap = squeeze(diffIRm(2:end, r1a_slice,:));
surf( x, y, surfMap,'FaceColor','interp',...
   'EdgeColor','none')
xlabel('M_{0B}', 'FontSize', FSz, 'FontWeight', 'bold')
ylabel('R_{1B} (s^{-1})', 'FontSize', FSz, 'FontWeight', 'bold')
zlabel( '% Diff to T_{1,Long}' , 'FontSize', FSz, 'FontWeight', 'bold');
ax = gca;  ax.FontSize = FSz; 
clim([ climLow, climHigh])
colormap(bipolar)
cbh = colorbar;
set(cbh,'YTick',linspace(climLow,climHigh, 5))
title('Mono IR')
axis tight
view(0,90)

set(gcf,'Position',[100 100 1100 400])
sgtitle( 'R_{1A} = 1.0 s^{-1}' , 'FontSize', FSz+4, 'FontWeight', 'bold')
saveas(gcf,fullfile(outputPathFig,[savePrefix,'_R1A_1.png']));












