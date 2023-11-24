% Look into the impact of MT effects on T1 mapping
% compare IR, MP2RAGE and VFA
addpath(genpath('E:\GitHub\NeuroImagingMatlab\QuantitativeFitting'))
addpath(genpath('E:\GitHub\NeuroImagingMatlab\colourmaps\OptimalMaps'))

% add necessary paths:
setupSimPaths;

clear all;
clc;

outputPath = 'E:\GitHub\Bloch_simulation\Investigations\MTeffect_in_T1mapping\savedOutputs\matrices';
outputPathFig = 'E:\GitHub\Bloch_simulation\Investigations\MTeffect_in_T1mapping\savedOutputs';
savePrefix ='FixedT1_r1'; % 'Sled2001', 'LorentzT270',

% Notes on savePrefix
% 'Sled2001' uses the relationship to calculate Ra from Raobs
% 'FixedT1' uses the assumption that the T1 of all pools is equal
% 'LorentzT270' changes lineshape to Lorenztian with T2b of 70us as in Van Gelderen et al 2016
% 'VanGeld2016' Builds on the above with  R1b = 4.0, R1a = 0.41 and
%   inversion staturation proportion of 87% 

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


if strcmp(savePrefix, "LorentzT270")
    % copying Van Gelderen et al 2016
    Params.lineshape = 'Lorentzian';
    Params.T2b = 70e-6;
elseif strcmp(savePrefix, "VanGeld2016")
    Params.lineshape = 'Lorentzian';
    Params.T2b = 70e-6;
    Params.Ra = 0.41;
    Params.R1b = 4;
    Params.Rrfb_inv = 665; % hack to get saturation level to 0.13 remaining
end

%% To test the MT effect, we will simulate for a range of bound pool parameters

M0b = 0:0.03:0.24;
T1 = linspace(600, 2500, 50);

[M0b_m, T1_m] = meshgrid(M0b, T1);

% Then vectorize for easy looping:
M0b_m = M0b_m(:); T1_m = T1_m(:)./1000;

% We actually input Ra_obs in, so invert:
Raobs_m = 1./T1_m;

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
    if strcmp(savePrefix, 'Sled2001') || strcmp(savePrefix, 'LorentzT270')...
            || strcmp(savePrefix, 'VanGeld2016')
        [vfa_sig(i,1), ~, ~] = BlochSimFlashSequence_v2( Params,...
            'Raobs', Raobs_m(i), 'M0b', M0b_m(i) );
    elseif strcmp( savePrefix, 'FixedT1')
        [vfa_sig(i,1), ~, ~] = BlochSimFlashSequence_v2( Params,...
            'Raobs', Raobs_m(i), 'M0b', M0b_m(i),'Ra', Raobs_m(i));
    end
end

% Image 2 - T1w
Params.TR = 15/1000;
Params.flipAngle = 20;

for i = 1: simNum
    if strcmp(savePrefix, 'Sled2001') || strcmp(savePrefix, 'LorentzT270')...
            || strcmp(savePrefix, 'VanGeld2016')
        [vfa_sig(i,2), ~, ~] = BlochSimFlashSequence_v2( Params,...
            'Raobs', Raobs_m(i), 'M0b', M0b_m(i) );
    elseif strcmp( savePrefix, 'FixedT1')
        [vfa_sig(i,2), ~, ~] = BlochSimFlashSequence_v2( Params,...
            'Raobs', Raobs_m(i), 'M0b', M0b_m(i),'Ra', Raobs_m(i));
    end
end
toc

vfa_T1 = Helms_VFA_T1_2008( vfa_sig(:,1), vfa_sig(:,2), 6, 20, 27, 15, 1);


save( fullfile( outputPath, [savePrefix,'_VFA_sim.mat']), "vfa_sig");
save( fullfile( outputPath, [savePrefix,'_vfa_T1.mat']), "vfa_T1");

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Simulate Inversion Recovery Signal.

TI = [30.0000; 250.0000; 500.0000; 750.0000; 1000.0000; 1500.0000]/1000;
numTI = length(TI);

% New parameters:
Params.InvPulseDur = 3/1000;
Params.Readout = 'linear';
Params.TR = 2000/1000;
Params.flipAngle = 90;
Params.numExcitation = 1;
Params.TE = 8.5/1000;
Params.CalcVector = 1;
Params.PerfectSpoiling = 1;
Params.DummyEcho = 0;
Params.echoSpacing = 0;

IR_sig = zeros( simNum, numTI);

tic % Roughly 1 min for sim and fitting
for i = 1: simNum
    for j = 1:numTI
        
        if strcmp(savePrefix, 'Sled2001') || strcmp(savePrefix, 'LorentzT270')...
            || strcmp(savePrefix, 'VanGeld2016')
            [IR_sig(i,j), ~, ~] = BlochSim_MPRAGESequence( Params,...
                'Raobs', Raobs_m(i), 'M0b', M0b_m(i), 'TI', TI(j) );
        elseif strcmp( savePrefix, 'FixedT1')
            [IR_sig(i,j), ~, ~] = BlochSim_MPRAGESequence( Params,...
                'Raobs', Raobs_m(i), 'M0b', M0b_m(i), 'TI', TI(j),...
                'Ra', Raobs_m(i) );
        end
    end
end
toc

tic
[IR_T1, resmap] = fit_T1_IR_data(IR_sig, TI,'complex');
toc

IR_T1 = IR_T1 *1000;


save( fullfile(outputPath,[savePrefix,'_IR_sim.mat']),"IR_sig");
save( fullfile(outputPath,[savePrefix,'_IR_T1.mat']),"IR_T1");




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Simulate MP2RAGE

% New Params
Params.InvPulseDur = 10.24/1000;
Params.invShape = 'hsn';
Params.PulseOpt = [];
Params.simAdiabatic = 1;

% Old
Params.TR = 5000/1000;
Params.flipAngle = [4,5];
Params.numExcitation = 175;
Params.echoSpacing = 7.7/1000;
Params.Readout = 'linear';
Params.TI = [940, 2830]./1000;
Params.DummyEcho = 0;
Params.PerfectSpoiling = 1;


MP2_sig = zeros( simNum, 2);

tic % roughly 3 mins
for i = 1: simNum

    if strcmp(savePrefix, 'Sled2001') || strcmp(savePrefix, 'LorentzT270')...
            || strcmp(savePrefix, 'VanGeld2016')
        [MP2_sig(i,1),MP2_sig(i,2), ~, ~] = BlochSim_MP2RAGESequence( Params,...
            'Raobs', Raobs_m(i), 'M0b', M0b_m(i) );
    elseif strcmp( savePrefix, 'FixedT1')
        [MP2_sig(i,1),MP2_sig(i,2), ~, ~] = BlochSim_MP2RAGESequence( Params,...
            'Raobs', Raobs_m(i), 'M0b', M0b_m(i),'Ra', Raobs_m(i));
    end

end
toc

MP2RAGE.B0          = Params.B0;                  % In Tesla
MP2RAGE.TR          = Params.TR;                  % MP2RAGE TR in seconds
MP2RAGE.TRFLASH     = Params.echoSpacing;             % TR of the GRE readout
MP2RAGE.TIs         = Params.TI;   % Inversion times - time between middle of refocusing pulse and excitatoin of the k-space center encoding
MP2RAGE.NZslices    = Params.numExcitation;            % Slices Per Slab * [PartialFourierInSlice-0.5  0.5]
MP2RAGE.FlipDegrees = Params.flipAngle;              % Flip angle of the two readouts in degrees

B1.img = ones(size(MP2_sig(:,2)));
brain.img = ones(size(MP2_sig(:,2)));
MP2RAGEINV2img.img = MP2_sig(:,2);
MP2RAGEimg.img = calculate_UNI_from_sims( MP2_sig(:,1), MP2_sig(:,2));
[ T1map, ~, ~] = CR_T1B1correctpackageTFL_withM0( B1, MP2RAGEimg, MP2RAGEINV2img, MP2RAGE, brain, 0.96);

mp2rage_T1 = T1map.img;

save( fullfile(outputPath,[savePrefix,'_MP2_sim.mat']),"MP2_sig");
save( fullfile(outputPath,[savePrefix,'_mp2rage_T1.mat']),"mp2rage_T1");

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Make Plots!
% We are interested how the measured T1 varies as a function of the bound
% pool fraction. So plot x-axis as true T1, y-axis = measured T1, and then
% do separate lines for each bound pool fraction.


%% Set up values:
FontSize = 20;
lLim = 600;
hLim = 1800;
ticV = lLim:400:hLim;

% x data = T1:
x = T1(:);
refy = x;
[cm_data]=inferno(12); % First one is black so lets avoid the darker shades so the reference line can be black. 


%% VFA
% y data = matrix with length(T1) rows, length(M0b) columns
y = reshape(vfa_T1, [length(T1), length(M0b)]);

figure; 
hold on;
for i = 1:length(M0b)
    plot(x,y(:,i),'Color', cm_data(end-length(M0b)+i-1,: ),'LineWidth',2)
end
xlabel( 'Input T_{1,obs} (ms)', 'FontSize', FontSize, 'FontWeight', 'bold' )
ylabel( 'Simulated T_{1,obs} (ms)' , 'FontSize', FontSize, 'FontWeight', 'bold');
xlim([lLim, hLim]); ylim([lLim, hLim])
xticks(ticV); yticks(ticV)
title('VFA', 'FontSize', FontSize, 'FontWeight', 'bold')
ax = gca;  ax.FontSize = FontSize; 
hold on
plot(x,refy,'--','Color',[0,0,0],'LineWidth',1.5)

% Add lines for typical GM and WM T1s
point = [850, 850];
axLims = [lLim hLim lLim hLim];  %[x-min, x-max, y-min, y-max] axis limits
plot([point(1), point(1)], [axLims(3), point(2)], 'k-')  %vertical line
plot([axLims(1), point(1)], [point(2), point(2)], 'k-')  %horizontal line
text(lLim + 20, point(1) + 70, strcat(num2str(point(1)), " ms"),'FontSize', FontSize-2)

point = [1400, 1400];
axLims = [lLim hLim lLim hLim];  %[x-min, x-max, y-min, y-max] axis limits
plot([point(1), point(1)], [axLims(3), point(2)], 'k-')  %vertical line
plot([axLims(1), point(1)], [point(2), point(2)], 'k-')  %horizontal line
text(lLim + 20, point(1) + 70, strcat(num2str(point(1)), " ms"),'FontSize', FontSize-2)
hold off

saveas(gcf,fullfile(outputPathFig,[savePrefix,'_T1plot_vfa.png']));


%% IR
% y data = matrix with length(T1) rows, length(M0b) columns
y = reshape(IR_T1, [length(T1), length(M0b)]);

figure; 
hold on;
for i = 1:length(M0b)
    plot(x,y(:,i),'Color', cm_data(end-length(M0b)+i-1,: ),'LineWidth',2)
end
xlabel( 'Input T_{1,obs} (ms)', 'FontSize', FontSize, 'FontWeight', 'bold' )
ylabel( 'Simulated T_{1,obs} (ms)' , 'FontSize', FontSize, 'FontWeight', 'bold');
xlim([lLim, hLim]); ylim([lLim, hLim])
xticks(ticV); yticks(ticV)
title('IR', 'FontSize', FontSize, 'FontWeight', 'bold')
ax = gca;  ax.FontSize = FontSize; 
hold on
plot(x,refy,'--','Color',[0,0,0],'LineWidth',1.5)

% Add lines for typical GM and WM T1s
point = [850, 850];
axLims = [lLim hLim lLim hLim];  %[x-min, x-max, y-min, y-max] axis limits
plot([point(1), point(1)], [axLims(3), point(2)], 'k-')  %vertical line
plot([axLims(1), point(1)], [point(2), point(2)], 'k-')  %horizontal line
text(lLim + 20, point(1) + 70, strcat(num2str(point(1)), " ms"),'FontSize', FontSize-2)

point = [1400, 1400];
axLims = [lLim hLim lLim hLim];  %[x-min, x-max, y-min, y-max] axis limits
plot([point(1), point(1)], [axLims(3), point(2)], 'k-')  %vertical line
plot([axLims(1), point(1)], [point(2), point(2)], 'k-')  %horizontal line
text(lLim + 20, point(1) + 70, strcat(num2str(point(1)), " ms"),'FontSize', FontSize-2)

hold off

saveas(gcf,fullfile(outputPathFig,[savePrefix,'_T1plot_IR.png']));

%% MP2RAGE
% y data = matrix with length(T1) rows, length(M0b) columns
y = reshape(mp2rage_T1, [length(T1), length(M0b)]);

figure; 
hold on;
for i = 1:length(M0b)
    plot(x,y(:,i),'Color', cm_data(end-length(M0b)+i-1,: ),'LineWidth',2)
end
xlabel( 'Input T_{1,obs} (ms)', 'FontSize', FontSize, 'FontWeight', 'bold' )
ylabel( 'Simulated T_{1,obs} (ms)' , 'FontSize', FontSize, 'FontWeight', 'bold');
xlim([lLim, hLim]); ylim([lLim, hLim])
xticks(ticV); yticks(ticV)
title('MP2RAGE', 'FontSize', FontSize, 'FontWeight', 'bold')
ax = gca;  ax.FontSize = FontSize; 
hold on
plot(x,refy,'--','Color',[0,0,0],'LineWidth',1.5)

% Add lines for typical GM and WM T1s
point = [850, 850];
axLims = [lLim hLim lLim hLim];  %[x-min, x-max, y-min, y-max] axis limits
plot([point(1), point(1)], [axLims(3), point(2)], 'k-')  %vertical line
plot([axLims(1), point(1)], [point(2), point(2)], 'k-')  %horizontal line
text(lLim + 20, point(1) + 70, strcat(num2str(point(1)), " ms"),'FontSize', FontSize-2)

point = [1400, 1400];
axLims = [lLim hLim lLim hLim];  %[x-min, x-max, y-min, y-max] axis limits
plot([point(1), point(1)], [axLims(3), point(2)], 'k-')  %vertical line
plot([axLims(1), point(1)], [point(2), point(2)], 'k-')  %horizontal line
text(lLim + 20, point(1) + 70, strcat(num2str(point(1)), " ms"),'FontSize', FontSize-2)

hold off

saveas(gcf,fullfile(outputPathFig,[savePrefix,'_T1plot_mp2rage.png']));


%% Difference From Gold Standard:
lLim2 = -30;
hLim2 = 30; %percents

% VFA - IR

y1 = reshape(vfa_T1, [length(T1), length(M0b)]);
y2 = reshape(IR_T1, [length(T1), length(M0b)]);

y = (y2-y1)./y2 *100; % percent difference

figure; 
hold on;
for i = 1:length(M0b)
    plot(x,y(:,i),'Color', cm_data(end-length(M0b)+i-1,: ),'LineWidth',2)
end
xlabel( 'Input T_{1,obs} (ms)', 'FontSize', FontSize, 'FontWeight', 'bold' )
ylabel( 'Sim. T_{1,obs} Diff (%)' , 'FontSize', FontSize, 'FontWeight', 'bold');
xlim([lLim, hLim]); ylim([lLim2, hLim2])
xticks(ticV)
title('IR - VFA', 'FontSize', FontSize, 'FontWeight', 'bold')
ax = gca;  ax.FontSize = FontSize; 
hold on
plot(x,refy,'--','Color',[0,0,0],'LineWidth',1.5)

% Add lines for typical GM and WM T1s
xline(850)
text(850 + 20, hLim2 - 5, strcat(num2str(850), " ms"),'FontSize', FontSize-2)

xline(1400)
text(1400 + 20, hLim2 - 5, strcat(num2str(point(1)), " ms"),'FontSize', FontSize-2)

hold off

saveas(gcf,fullfile(outputPathFig,[savePrefix,'_T1plot_diff_IR_vs_vfa.png']));


% MP2RAGE - IR

y1 = reshape(mp2rage_T1, [length(T1), length(M0b)]);
y2 = reshape(IR_T1, [length(T1), length(M0b)]);

y = (y2-y1)./y2 *100; % percent difference

figure; 
hold on;
for i = 1:length(M0b)
    plot(x,y(:,i),'Color', cm_data(end-length(M0b)+i-1,: ),'LineWidth',2)
end
xlabel( 'Input T_{1,obs} (ms)', 'FontSize', FontSize, 'FontWeight', 'bold' )
ylabel( 'Sim. T_{1,obs} Diff (%)' , 'FontSize', FontSize, 'FontWeight', 'bold');
xlim([lLim, hLim]); ylim([lLim2, hLim2])
xticks(ticV);
title('IR - MP2RAGE', 'FontSize', FontSize, 'FontWeight', 'bold')
ax = gca;  ax.FontSize = FontSize; 
hold on
plot(x,refy,'--','Color',[0,0,0],'LineWidth',1.5)

% Add lines for typical GM and WM T1s
xline(850)
text(850 + 20, hLim2 - 5, strcat(num2str(850), " ms"),'FontSize', FontSize-2)

xline(1400)
text(1400 + 20, hLim2 - 5, strcat(num2str(point(1)), " ms"),'FontSize', FontSize-2)

hold off

saveas(gcf,fullfile(outputPathFig,[savePrefix,'_T1plot_diff_IR_vs_mp2rage.png']));














