%% Params to implement for IR:
Params.B0 = 3;
Params.TissueType = 'WM';
Params = DefaultCortexTissueParams(Params);
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
Params.Inv.PulseOpt.nSamples = 512;
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

%%
% Temp notes on sim with adiabatic:


alpha = 90; % degrees
Trf = 3.072/1000;
PulseOpt.TBW = 4;  
nSamples = 512;

excPulse = GetPulse(alpha, 0, Trf, 'sinc', PulseOpt);
%ViewPulse(excPulse, 'omega');


t = 0:Trf/(nSamples-1):Trf;
rf =excPulse.b1(t);
ftP = fftshift(fft(rf)); ftP = ftP/(max(abs(ftP)));
fshift = (-nSamples/2:nSamples/2-1)/Trf; % center our frequency around 0

figure; tiledlayout(1,2); nexttile;
    plot(t*1000, abs(rf), 'LineWidth', 3); 
    xlabel('Time(ms)'); ylabel('B_1 (μT)');
    title('B_1');ax = gca; ax.FontSize = 20;
    
    nexttile;
    plot( fshift,abs(ftP), 'LineWidth', 3);
    xlabel('Freq (Hz)'); ylabel('Energy'); xlim([-3000, 3000])
    title('Frequency');ax = gca; ax.FontSize = 20;
    set(gcf,'Position',[100 100 1200 600])

% Need the B1rms:
dt = Trf/(nSamples-1);
RF = (alpha*pi/180)* b1/sum(b1) /(2*pi*gamma*dt);
figure; plot(t*1000, (RF).^2, 'LineWidth', 3); 

B1rms90 = 1/Trf*trapz(t,RF.^2);


B1_time = excPulse.omega(t)/(2*pi*gam);
P3 = (1/Trf) * trapz(t,B1_time.^2);
B1rms = sqrt(P3);

% plot(t*1000, abs(b1), 'LineWidth', 3); 



%% BLoch Sim
Params.Ra = 1;
Params.R2a = 1000/80;
Params.PulseOpt.nSamples = 512;
B = [0 0 1]'; %, Params.R1b*Params.M0b, 0

t = 0:Trf/(nSamples-1):Trf;
rf_pulse =1000*excPulse.amp*excPulse.b1(t);
plot(t*1000, abs(rf_pulse), 'LineWidth', 3); 

M_return = blochSimAdiabaticPulse_1pool( rf_pulse, Trf, 0,...
                                            Params, B, B);


















































%% Define the adiabatic inversion pulse, non-selectice:
Params.B0 = 3;
Params.TissueType = 'WM';
Params = DefaultCortexTissueParams(Params);


Params.M0a = 1;
Params.M0b = 0.1;
Params.PulseOpt.beta = 672;  
Params.PulseOpt.n = 1; 
Params.PulseOpt.mu = 5;
Params.PulseOpt.A0 = 17.5;
Params.PulseOpt.nSamples = 512;
dispFig = 1;
Trf = 10.24/1000;
shape = 'hsn';

[rf_pulse, Params] = GetAdiabaticPulse( Trf, shape, dispFig, Params);

%% IR Define the adiabatic inversion pulse, slice-selectice:
Params.B0 = 3;
Params.TissueType = 'WM';
Params = DefaultCortexTissueParams(Params);


Params.M0a = 1;
Params.M0b = 0.1;
Params.PulseOpt.beta = 672;  
Params.PulseOpt.n = 1; 
Params.PulseOpt.mu = 5;
Params.PulseOpt.A0 = 17.5;
Params.PulseOpt.nSamples = 512;
dispFig = 1;
Trf = 10.24/1000;
shape = 'hsn';

[rf_pulse, Params] = GetAdiabaticPulse( Trf, shape, dispFig, Params);

% Need the B1rms:



%% IR Define the 90 degree sinc pulse, slice-selectice:

alpha = 90; % degrees
Trf = 3.072/1000;
PulseOpt.TBW = 4;  
nSamples = 512;

excPulse = GetPulse(alpha, 0, Trf, 'sinc', PulseOpt);
%ViewPulse(excPulse, 'omega');


t = 0:Trf/(nSamples-1):Trf;
rf =excPulse.b1(t);
ftP = fftshift(fft(rf)); ftP = ftP/(max(abs(ftP)));
fshift = (-nSamples/2:nSamples/2-1)/Trf; % center our frequency around 0

figure; tiledlayout(1,2); nexttile;
plot(t*1000, abs(rf), 'LineWidth', 3); 
xlabel('Time(ms)'); ylabel('B_1 (μT)');
title('B_1');ax = gca; ax.FontSize = 20;

nexttile;
plot( fshift,abs(ftP), 'LineWidth', 3);
xlabel('Freq (Hz)'); ylabel('Energy'); xlim([-3000, 3000])
title('Frequency');ax = gca; ax.FontSize = 20;
set(gcf,'Position',[100 100 1200 600])

% Need the B1rms:
gamma = 2*pi*42.576;
b1 = excPulse.omega(t) /gamma; % in microTesla

% plot(t*1000, abs(b1), 'LineWidth', 3); 


p2 = trapz(t,rf.^2)/ Params.pulseDur; % See Soustelle et al 2022 for def.
B1peak = sqrt( Psat/p2);

if B1peak > B1peak_limit % conform to hardware constraint
    B1peak = B1peak_limit;
end

Params.satFlipAngle = trapz( t, rf*B1peak)*gam * 360;
Params.satB1peak = B1peak*1e6; % convert to microTesla




%% IR Define the adiabatic refocusing pulse, slice-selectice:
% Looks like HSn, but amplitude only?



Params.PulseOpt.beta = 672;  
Params.PulseOpt.n = 1; 
Params.PulseOpt.mu = 5;
Params.PulseOpt.A0 = 17.5;
Params.PulseOpt.nSamples = 512;
dispFig = 1;
Trf = 10.24/1000;
shape = 'hsn';

[rf_pulse, Params] = GetAdiabaticPulse( Trf, shape, dispFig, Params);



    case 'hard';      pulse_fcn = @hard_pulse;  
    case 'sinc';      pulse_fcn = @sinc_pulse;        
    case 'sinchann';  pulse_fcn = @sinchann_pulse;        
    case 'sincgauss'; pulse_fcn = @sincgauss_pulse;        
    case 'gaussian';  pulse_fcn = @gaussian_pulse;        
    case 'gausshann'; pulse_fcn = @gausshann_pulse;    




alpha = 180; % degrees
Trf = 3.000/1000;
PulseOpt.TBW = 4;  
nSamples = 512;

excPulse = GetPulse(alpha, 0, Trf, 'gaussian');
figure; ViewPulse(excPulse, 'b1');


t = 0:Trf/(nSamples-1):Trf;
rf =excPulse.b1(t);
ftP = fftshift(fft(rf)); ftP = ftP/(max(abs(ftP)));
fshift = (-nSamples/2:nSamples/2-1)/Trf; % center our frequency around 0

figure; tiledlayout(1,2); nexttile;
plot(t*1000, abs(rf), 'LineWidth', 3); 
xlabel('Time(ms)'); ylabel('B_1 (μT)');
title('B_1');ax = gca; ax.FontSize = 20;

nexttile;
plot( fshift,abs(ftP), 'LineWidth', 3);
xlabel('Freq (Hz)'); ylabel('Energy'); xlim([-3000, 3000])
title('Frequency');ax = gca; ax.FontSize = 20;
set(gcf,'Position',[100 100 1200 600])


% Need the B1rms:
























% You can Bloch sim using:
b1Rel = 1;
freqOff = 0;
M_return = blochSimAdiabaticPulse( b1Rel* rf_pulse,...
                Trf, PulseOpt, freqOff);

% If you have (Ma x,y,z + Mbz) then you need to convert both to be 3dim
% Assume Bpool XY are 0. If we have dipolar order, we will assume 0.
Mstart = [1,2,3,4,5];
Msim = zeros(6,1);
Msim([1,3]) = Mstart(1:2);
Msim(5:6) = Mstart(3:4);

% Convert them back
Mstart(1:2) = Msim([1,3]);
Mstart(3:4) = Msim(5:6);

%% Trial run!

Params.B0 = 3;
Params.MTC = 0; % Magnetization Transfer Contrast
Params.TissueType = 'WM';
Params.numExcitation = 1;
Params.echoSpacing = 7.7/1000;

Params = DefaultCortexTissueParams(Params);
Params = CalcImagingParams(Params);
Params = CalcVariableImagingParams(Params);

Params.PerfectSpoiling = 1;

Params.InvPulseDur = 10.24/1000;
Params.invShape = 'hsn';
Params.PulseOpt = [];
Params.simAdiabatic = 1;
Params.CalcVector = 1;

% Old
Params.TR = 5000/1000;
Params.flipAngle = [4,5];
Params.numExcitation = 175;
Params.echoSpacing = 7.7/1000;
Params.Readout = 'linear';
Params.TI = [940, 2830]./1000;
Params.DummyEcho = 0;
Params.PerfectSpoiling = 1;


[outSig1, outSig2, M, time_vect] = BlochSim_MP2RAGESequence(Params);

figure; plot(time_vect,M(3,:), 'LineWidth',3);
hold on
plot(time_vect,M(4,:), 'LineWidth',3);
xlim([0 3*Params.TR])
% figure; plot(time_vect,M(5,:), 'LineWidth',3);














