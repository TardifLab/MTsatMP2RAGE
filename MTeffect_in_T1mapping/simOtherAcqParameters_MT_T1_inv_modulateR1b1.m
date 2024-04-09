%% The goal here is that we want to look at different protocols to compare
% them! 

% add necessary paths:
setupSimPaths;

clear all;
clc;

outputPath = 'E:\GitHub\MTsatMP2RAGE\MTeffect_in_T1mapping\savedOutputs\matrices';
outputPathFig = 'E:\GitHub\MTsatMP2RAGE\MTeffect_in_T1mapping\savedOutputs';
savePrefix ='SeqCompare'; 

% Notes on savePrefix
% ModulateTP -> not realistic for all values to be maintained, so we use
% linear models to vary them.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% VFA First
 
Params.B0 = 3;
Params.MTC = 0; % Magnetization Transfer Contrast
Params.TissueType = 'WM';
Params.numExcitation = 1;
Params.echoSpacing = 7.7/1000;

Params = DefaultCortexTissueParams(Params);
Params = CalcImagingParams(Params);
Params = CalcVariableImagingParams(Params);

Params.PerfectSpoiling = 1;
Params.R1b = 1; % The assumed value for Sled and Pike
Params.Ra = 1; % The assumed value for Sled and Pike
Params.R = []; % clear it to avoid errors

%% To test the MT effect, we will simulate for a range of bound pool parameters
% Compared to the main script, here we will vary the values based on Sled
% and Pike 2001. 


T1 = linspace(600, 2500, 50);
% Then vectorize for easy looping:
T1_m = T1(:)./1000;

% We actually input Ra_obs in, so invert:
Raobs_m = 1./T1_m;
[M0b, kf, Ra, T2a, T2b, kr] = qMT_LinearModelVals( Raobs_m);
simNum = length(Raobs_m);


%% Baudrexel et al 2018 -> VFA2
vfa_sig = zeros( simNum, 2);
% Image 1 - PDw
Params.TR = 16.4/1000;
Params.flipAngle = 4;

tic % About 18 seconds for 450 options
for i = 1 : simNum
    [vfa_sig(i,1), ~, ~] = BlochSimFlashSequence_v2( Params,...
        'M0b', M0b(i),'T2a', T2a(i), 'T2b', T2b(i),'kf', kf(i), 'kr', kr(i));
end

% Image 2 - T1w
Params.TR = 16.4/1000;
Params.flipAngle = 24;
for i = 1: simNum
        [vfa_sig(i,2), ~, ~] = BlochSimFlashSequence_v2( Params,...
           'M0b', M0b(i),'T2a', T2a(i), 'T2b', T2b(i),'kf', kf(i), 'kr', kr(i) );
end
toc

vfa_T1_2 = Helms_VFA_T1_2008( vfa_sig(:,1), vfa_sig(:,2), 4, 24, 16.4, 16.4, 1);

save( fullfile( outputPath, [savePrefix,'_VFA_sim_Baudrexel.mat']), "vfa_sig");
save( fullfile( outputPath, [savePrefix,'_vfa_T1_Baudrexel.mat']), "vfa_T1_2");

%% Weiskopf et al 2013 - MPM -> VFA3
vfa_sig = zeros( simNum, 2);
% Image 1 - PDw
Params.TR = 23.7/1000;
Params.flipAngle = 6;
tic % About 18 seconds for 450 options
for i = 1 : simNum
    [vfa_sig(i,1), ~, ~] = BlochSimFlashSequence_v2( Params,...
       'M0b', M0b(i),'T2a', T2a(i), 'T2b', T2b(i),'kf', kf(i), 'kr', kr(i));
end

% Image 2 - T1w
Params.TR = 18.7/1000;
Params.flipAngle = 20;
for i = 1: simNum
        [vfa_sig(i,2), ~, ~] = BlochSimFlashSequence_v2( Params,...
            'M0b', M0b(i),'T2a', T2a(i), 'T2b', T2b(i),'kf', kf(i), 'kr', kr(i) );
end
toc

vfa_T1_3 = Helms_VFA_T1_2008( vfa_sig(:,1), vfa_sig(:,2), 6, 20, 23.7, 18.7, 1);

save( fullfile( outputPath, [savePrefix,'_VFA_sim_Weiskopf.mat']), "vfa_sig");
save( fullfile( outputPath, [savePrefix,'_vfa_T1_Weiskopf.mat']), "vfa_T1_3");

%%  Helms et al 2008 - Quantitative FlashMRI -> VFA4
vfa_sig = zeros( simNum, 2);
% Image 1 - PDw
Params.TR = 30/1000;
Params.flipAngle = 7;
tic % About 18 seconds for 450 options
for i = 1 : simNum
    [vfa_sig(i,1), ~, ~] = BlochSimFlashSequence_v2( Params,...
        'M0b', M0b(i),'T2a', T2a(i), 'T2b', T2b(i),'kf', kf(i), 'kr', kr(i));
end

% Image 2 - T1w
Params.TR = 30/1000;
Params.flipAngle = 20;
for i = 1: simNum
        [vfa_sig(i,2), ~, ~] = BlochSimFlashSequence_v2( Params,...
            'M0b', M0b(i),'T2a', T2a(i), 'T2b', T2b(i),'kf', kf(i), 'kr', kr(i) );
end
toc

vfa_T1_4 = Helms_VFA_T1_2008( vfa_sig(:,1), vfa_sig(:,2), 7, 20, 30, 30, 1);

save( fullfile( outputPath, [savePrefix,'_VFA_sim_Helms.mat']), "vfa_sig");
save( fullfile( outputPath, [savePrefix,'_vfa_T1_Helms.mat']), "vfa_T1_4");



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%        
%% MP2RAGE

%% Marques et al 2010 -3T -> mp2rage_T1_2
clear Params

Params.B0 = 3;
Params.TissueType = 'WM';
Params = DefaultCortexTissueParams(Params);
Params.PerfectSpoiling = 1;
Params.WExcDur = 0.1/1000;
Params.PerfectSpoiling = 1;
Params.R1b = 1; % The assumed value for Sled and Pike
Params.Ra = 1; % The assumed value for Sled and Pike

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

Params.TR = 6.25;
Params.flipAngle = [7,5];
Params.numExcitation = 160;
Params.echoSpacing = 5.8/1000;
Params.Readout = 'linear';
Params.TI = [800, 2200]./1000;
Params.DummyEcho = 0;
Params.PerfectSpoiling = 1;
Params.RFspoiling = 0;

MP2_sig = zeros( simNum, 2);

tic % roughly 3 mins
for i = 1: simNum
    [MP2_sig(i,1),MP2_sig(i,2), ~, ~] = BlochSim_MP2RAGESequence_6vec( Params,...
       'M0b', M0b(i),'T2a', T2a(i), 'T2b', T2b(i),'kf', kf(i), 'kr', kr(i) );
end
toc

MP2RAGE.B0          = Params.B0;                  % In Tesla
MP2RAGE.TR          = Params.TR;                  % MP2RAGE TR in seconds
MP2RAGE.TRFLASH     = Params.echoSpacing;             % TR of the GRE readout
MP2RAGE.TIs         = Params.TI;   % Inversion times - time between middle of refocusing pulse and excitatoin of the k-space center encoding
MP2RAGE.NZslices    = Params.numExcitation;            % Slices Per Slab * [PartialFourierInSlice-0.5  0.5]
MP2RAGE.FlipDegrees = Params.flipAngle;              % Flip angle of the two readouts in degrees

[mp2rage_T1_2, ~, ~] = MP2RAGE_dictionaryMatching(MP2RAGE, MP2_sig(:,1), MP2_sig(:,2),...
    ones(size(MP2_sig(:,2))), [0.0005, 0.005], 0);
mp2rage_T1_2 = mp2rage_T1_2*1000;

save( fullfile(outputPath,[savePrefix,'_MP2_sim_Marques_Optimal.mat']),"MP2_sig");
save( fullfile(outputPath,[savePrefix,'_mp2rage_T1_Marques_Optimal.mat']),"mp2rage_T1_2");




%% Marques et al 2010 -3T short TR -> mp2rage_T1_3
    Params.TR = 5;
    Params.flipAngle = [7,5];
    Params.numExcitation = 160;
    Params.echoSpacing = 5.8/1000;
    Params.Readout = 'linear';
    Params.TI = [700, 2500]./1000;
    Params.DummyEcho = 0;
    Params.PerfectSpoiling = 1;

    MP2_sig = zeros( simNum, 2);

tic % roughly 3 mins
for i = 1: simNum
    [MP2_sig(i,1),MP2_sig(i,2), ~, ~] = BlochSim_MP2RAGESequence_6vec( Params,...
        'M0b', M0b(i),'T2a', T2a(i), 'T2b', T2b(i),'kf', kf(i), 'kr', kr(i) );
end
toc

MP2RAGE.B0          = Params.B0;                  % In Tesla
MP2RAGE.TR          = Params.TR;                  % MP2RAGE TR in seconds
MP2RAGE.TRFLASH     = Params.echoSpacing;             % TR of the GRE readout
MP2RAGE.TIs         = Params.TI;   % Inversion times - time between middle of refocusing pulse and excitatoin of the k-space center encoding
MP2RAGE.NZslices    = Params.numExcitation;            % Slices Per Slab * [PartialFourierInSlice-0.5  0.5]
MP2RAGE.FlipDegrees = Params.flipAngle;              % Flip angle of the two readouts in degrees

[mp2rage_T1_3, ~, ~] = MP2RAGE_dictionaryMatching(MP2RAGE, MP2_sig(:,1), MP2_sig(:,2),...
    ones(size(MP2_sig(:,2))), [0.0005, 0.005], 0);
mp2rage_T1_3 = mp2rage_T1_3*1000;

save( fullfile(outputPath,[savePrefix,'_MP2_sim_Marques_shortTR.mat']),"MP2_sig");
save( fullfile(outputPath,[savePrefix,'_mp2rage_T1_Marques_shortTR.mat']),"mp2rage_T1_3");

%% Mussard et al 2020 -fully sampled data -> mp2rage_T1_4
    Params.TR = 5;
    Params.flipAngle = [4,5];
    Params.numExcitation = 195;
    Params.echoSpacing = 5.8/1000; % not provided. Using from current study
    Params.Readout = 'linear';
    Params.TI = [700, 2500]./1000;
    Params.DummyEcho = 0;
    Params.PerfectSpoiling = 1;

MP2_sig = zeros( simNum, 2);

tic % roughly 3 mins
for i = 1: simNum
    [MP2_sig(i,1),MP2_sig(i,2), ~, ~] = BlochSim_MP2RAGESequence_6vec( Params,...
        'M0b', M0b(i),'T2a', T2a(i), 'T2b', T2b(i),'kf', kf(i), 'kr', kr(i) );
end
toc

MP2RAGE.B0          = Params.B0;                  % In Tesla
MP2RAGE.TR          = Params.TR;                  % MP2RAGE TR in seconds
MP2RAGE.TRFLASH     = Params.echoSpacing;             % TR of the GRE readout
MP2RAGE.TIs         = Params.TI;   % Inversion times - time between middle of refocusing pulse and excitatoin of the k-space center encoding
MP2RAGE.NZslices    = Params.numExcitation;            % Slices Per Slab * [PartialFourierInSlice-0.5  0.5]
MP2RAGE.FlipDegrees = Params.flipAngle;              % Flip angle of the two readouts in degrees

[mp2rage_T1_4, ~, ~] = MP2RAGE_dictionaryMatching(MP2RAGE, MP2_sig(:,1), MP2_sig(:,2),...
    ones(size(MP2_sig(:,2))), [0.0005, 0.005], 0);
mp2rage_T1_4 = mp2rage_T1_4*1000;

save( fullfile(outputPath,[savePrefix,'_MP2_sim_Mussard.mat']),"MP2_sig");
save( fullfile(outputPath,[savePrefix,'_mp2rage_T1_Mussard.mat']),"mp2rage_T1_4");





%% Load the values from the rest of the manuscript for reference:
savePrefixOld ='R1B_1s'; % 'Sled2001', 'LorentzT270',

load( fullfile( outputPath, [savePrefixOld,'_vfa_T1_rev2.mat']));
load( fullfile(outputPath,[savePrefixOld,'_mp2rage_T1_rev2.mat']));
load( fullfile(outputPath,[savePrefixOld,'_IR_T1_biLong_rev2.mat']));



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Make Plots!
% We are interested how the measured T1 varies as a function of the bound
% pool fraction. So plot x-axis as true T1, y-axis = measured T1, and then
% do separate lines for each bound pool fraction.


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% plots using diff against T1 long


%% Set up values:
FontSize = 20;
lLim = 600;
hLim = 1800;
ticV = lLim:400:hLim;

ylLim = -40;
yhLim = 40;
yticV = ylLim:20:yhLim;

% x data = T1:
x = T1(:);
refy = x;
cm_data = 	[0, 0.4470, 0.9410; ...
    0.4660, 0.7740, 0.1880;...
    0.9290, 0.6940, 0.1250;...
    0.9350, 0.0780, 0.1840];

%% VFA methods

y = [(vfa_T1-IR_T1_biL)./IR_T1_biL,...
    (vfa_T1_2-IR_T1_biL)./IR_T1_biL,...
    (vfa_T1_3-IR_T1_biL)./IR_T1_biL,...
    (vfa_T1_4-IR_T1_biL)./IR_T1_biL]*100;


figure; 
hold on;
for i = 1:4
    plot(x,y(:,i),'Color', cm_data(i,: ),'LineWidth',2)
end
xlabel( 'Input T_{1,obs} (ms)', 'FontSize', FontSize, 'FontWeight', 'bold' )
ylabel( '% Diff. to T_{1,Long}' , 'FontSize', FontSize, 'FontWeight', 'bold');
xlim([lLim, hLim]); ylim([ylLim, yhLim])
xticks(ticV); yticks(yticV)
%title('VFA', 'FontSize', FontSize, 'FontWeight', 'bold')
ax = gca;  ax.FontSize = FontSize; 
hold on

legend('Current','Baudrexel et al','Weiskopf et al','Helms et al', 'FontSize',FontSize-5, 'location', 'best')
hold off

saveas(gcf,fullfile(outputPathFig,[savePrefix,'_VFA_T1diff_comparison.png']));



%% MP2RAGE methods

y = [(mp2rage_T1-IR_T1_biL)./IR_T1_biL,...
    (mp2rage_T1_2-IR_T1_biL)./IR_T1_biL,...
    (mp2rage_T1_3-IR_T1_biL)./IR_T1_biL,...
    (mp2rage_T1_4-IR_T1_biL)./IR_T1_biL]*100;

figure; 
hold on;
for i = 1:4
    plot(x,y(:,i),'Color', cm_data(i,: ),'LineWidth',2)
end
xlabel( 'Input T_{1,obs} (ms)', 'FontSize', FontSize, 'FontWeight', 'bold' )
ylabel( '% Diff. to T_{1,Long}' , 'FontSize', FontSize, 'FontWeight', 'bold');
xlim([lLim, hLim]); ylim([ylLim, yhLim])
xticks(ticV); yticks(yticV)
%title('VFA', 'FontSize', FontSize, 'FontWeight', 'bold')
ax = gca;  ax.FontSize = FontSize; 
hold on
plot(x,refy,':','Color',[0,0,0],'LineWidth',1.5)

legend('Current','Marques et al Optimal','Marques et al Short TR','Mussard et al', 'FontSize', FontSize-5, 'location', 'best')
hold off

saveas(gcf,fullfile(outputPathFig,[savePrefix,'_MP2RAGE_T1diff_comparison.png']));







%% These don't make sense without varying Ra more.  
% %% Set up values:
% FontSize = 20;
% x = T1(:);
% 
% lLim = 600;
% hLim = 1800;
% ticV = lLim:400:hLim;
% refy = x;
% cm_data = 	[0, 0.4470, 0.9410; ...
%     0.4660, 0.7740, 0.1880;...
%     0.9290, 0.6940, 0.1250;...
%     0.9350, 0.0780, 0.1840];

% %% VFA methods
% 
% y = [vfa_T1, vfa_T1_2, vfa_T1_3, vfa_T1_4];
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
% legend('Current','Baudrexel et al','Weiskopf et al','Helms et al', 'FontSize',FontSize-5, 'location', 'best')
% hold off
% 
% saveas(gcf,fullfile(outputPathFig,[savePrefix,'_VFA_T1plot_comparison.png']));
% 
% 
% 
% %% MP2RAGE methods
% 
% y = [mp2rage_T1, mp2rage_T1_2, mp2rage_T1_3, mp2rage_T1_4];
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
% legend('Current','Marques et al Optimal','Marques et al Short TR','Mussard et al', 'FontSize', FontSize-5, 'location', 'best')
% hold off
% 
% saveas(gcf,fullfile(outputPathFig,[savePrefix,'_MP2RAGE_T1plot_comparison.png']));















