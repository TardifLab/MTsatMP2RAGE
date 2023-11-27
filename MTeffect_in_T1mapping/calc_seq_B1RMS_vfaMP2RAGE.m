%% Calculate B1rms of sequences
% Using definition from Soustelle et al 2022 (ihMT)

% Start with the VFA
wExcDur = 0.1/1000; % ms
% PDw
TR = 27/1000;
flipAngle = 6;

excB1rms = getExcPulseB1rms(flipAngle, wExcDur);
    % B1rms over tr = pulse * sqrt(number pulses*pulseWidth/TR)
B1rmsTR = excB1rms * sqrt( 1 *wExcDur/ TR);  % -> 0.24 uT

% T1w
TR = 15/1000;
flipAngle = 20;

excB1rms = getExcPulseB1rms(flipAngle, wExcDur);
B1rmsTR = excB1rms * sqrt( 1 *wExcDur/ TR);  % -> 1.07 uT

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% MT-weighted image
TR = 27/1000;
flipAngle = 6;
wExcDur = 0.1/1000; % ms

excB1rms = getExcPulseB1rms(flipAngle, wExcDur);

% Sat pulse
B1rmsSat = 3.26;
tSat = 12/1000;
B1rmsTR = sqrt( B1rmsSat.^2 .*tSat./TR + excB1rms.^2 .*wExcDur./TR);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% MP2RAGE
Params.B0 = 3;
Params.TissueType = 'WM';
Params = DefaultCortexTissueParams(Params);
Params.Ra = 1;

Params.M0a = 1;
Params.M0b = 0.1;
Params.PulseOpt.beta = 672;  
Params.PulseOpt.n = 1; 
Params.PulseOpt.mu = 5;
Params.PulseOpt.A0 = 17.5;
Params.PulseOpt.nSamples = 512;

wExcDur = 0.1/1000; % ms
invDur = 10.24/1000; 

TR = 5000/1000;
flipAngle = [4,5];
numExcitation = 208;
echoSpacing = 7.7/1000;
TI = [940, 2830]./1000;

excB1rms1 = getExcPulseB1rms(flipAngle(1), wExcDur);
excB1rms2 = getExcPulseB1rms(flipAngle(2), wExcDur);

ExcitationBlock = numExcitation*echoSpacing;
TA = TI(1) - ExcitationBlock/2 - invDur;
TB = (TI(2) - ExcitationBlock/2) - (TI(1) + ExcitationBlock/2);
TC = TR - (TI(2) + ExcitationBlock/2);

if TR ~= (TA+TB+TC+invDur + 2*ExcitationBlock); error('times dont match'); end


% Block 1 -post inverstion (TA), RAGE1, time between (TB)
t1= TA + ExcitationBlock + TB;
B1rmsBlock1 = excB1rms1 * sqrt( numExcitation *wExcDur/ t1);  

% Block 2 - RAGE2, time delay (TC)
t2= ExcitationBlock + TC;
B1rmsBlock2 = excB1rms2 * sqrt( numExcitation *wExcDur/ t2); 

% Block 3 - inversion pulse:
invPulse = hyperbolicSecant_pulse( invDur, Params, dispFigure);

% rms
InvRMS = sqrt( (1/invDur) * trapz( linspace(0, invDur, length(invPulse)), abs(invPulse).^2) ); % 9.42 uT

% Combine!
B1rmsTR = sqrt( InvRMS.^2 .*invDur./TR + ...
                B1rmsBlock1.^2 .*t1./TR + ...
                B1rmsBlock2.^2 .*t2./TR);

% Apparently 0.31 uT?

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% cs - MP2RAGE
Params.B0 = 3;
Params.TissueType = 'WM';
Params = DefaultCortexTissueParams(Params);
Params.Ra = 1;

Params.M0a = 1;
Params.M0b = 0.1;
Params.PulseOpt.beta = 672;  
Params.PulseOpt.n = 1; 
Params.PulseOpt.mu = 5;
Params.PulseOpt.A0 = 17.5;
Params.PulseOpt.nSamples = 512;

wExcDur = 0.1/1000; % ms
invDur = 10.24/1000; 

TR = 5000/1000;
flipAngle = [4,5];
numExcitation = 175;
echoSpacing = 7.7/1000;
TI = [940, 2830]./1000;

excB1rms1 = getExcPulseB1rms(flipAngle(1), wExcDur);
excB1rms2 = getExcPulseB1rms(flipAngle(2), wExcDur);

ExcitationBlock = numExcitation*echoSpacing;
TA = TI(1) - ExcitationBlock/2 - invDur;
TB = (TI(2) - ExcitationBlock/2) - (TI(1) + ExcitationBlock/2);
TC = TR - (TI(2) + ExcitationBlock/2);

if TR ~= (TA+TB+TC+invDur + 2*ExcitationBlock); error('times dont match'); end


% Block 1 -post inverstion (TA), RAGE1, time between (TB)
t1= TA + ExcitationBlock + TB;
B1rmsBlock1 = excB1rms1 * sqrt( numExcitation *wExcDur/ t1);  % -> 1.07 uT

% Block 2 - RAGE2, time delay (TC)
t2= ExcitationBlock + TC;
B1rmsBlock2 = excB1rms2 * sqrt( numExcitation *wExcDur/ t2);  % -> 1.07 uT

% Block 3 - inversion pulse:
invPulse = hyperbolicSecant_pulse( invDur, Params, 1);

% rms
InvRMS = sqrt( (1/invDur) * trapz( linspace(0, invDur, length(invPulse)), abs(invPulse).^2) ); % 9.42 uT

% Combine!
B1rmsTR = sqrt( InvRMS.^2 .*invDur./TR + ...
                B1rmsBlock1.^2 .*t1./TR + ...
                B1rmsBlock2.^2 .*t2./TR);




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% IR

%% IR Define the adiabatic inversion pulse, slice-selectice:
Params.B0 = 3;
Params.TissueType = 'WM';
Params = DefaultCortexTissueParams(Params);
Params.Ra = 1;

Params.M0a = 1;
Params.M0b = 0.1;
Params.PulseOpt.beta = 672;  
Params.PulseOpt.n = 1; 
Params.PulseOpt.mu = 5;
Params.PulseOpt.A0 = 13.726;
Params.PulseOpt.nSamples = 512;
dispFig = 1;
Trf = 10.24/1000;
shape = 'hsn';

[rf_pulse, Params] = GetAdiabaticPulse( Trf, shape, dispFig, Params);

% Need the B1rms:
t = linspace(0, Trf, length(rf_pulse));
InvRMS = sqrt( (1/Trf) * trapz( t, abs(rf_pulse.^2) ) ); % 7.39 uT


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


% Bloch Sim to make sure I have 90 deg pulse:
Params.Ra = 1;
Params.R2a = 1000/80;
Params.PulseOpt.nSamples = 512;
B = [0 0 1]'; %, Params.R1b*Params.M0b, 0

t = 0:Trf/(nSamples-1):Trf;
rf_pulse =1000*excPulse.amp*excPulse.b1(t);
plot(t*1000, abs(rf_pulse), 'LineWidth', 3); 

M_return = blochSimAdiabaticPulse_1pool( rf_pulse, Trf, 0,...
                                            Params, B, B);
% Need the B1rms:
B1rms90 = sqrt( 1/Trf*trapz(t,rf_pulse.^2) ); % 4.13 uT



%% IR Define the adiabatic refocusing pulse, slice-selectice:
% Looks like Gaussian

alpha = 180; % degrees
Trf = 3.000/1000;
nSamples = 512;

excPulse = GetPulse(alpha, 0, Trf, 'gaussian');
% figure; ViewPulse(excPulse, 'b1');

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


% Bloch Sim to make sure I have 180 deg pulse:
Params.Ra = 1;
Params.R2a = 1000/80;
Params.PulseOpt.nSamples = 512;
B = [0 1 0]'; %, Params.R1b*Params.M0b, 0

t = 0:Trf/(nSamples-1):Trf;
rf_pulse =1000*excPulse.amp*excPulse.b1(t);
% figure; plot(t*1000, abs(rf_pulse), 'LineWidth', 3); 

M_return = blochSimAdiabaticPulse_1pool( rf_pulse, Trf, 0,...
                                            Params, B, B);
% Need the B1rms:
B1rms180 = sqrt( 1/Trf*trapz(t,rf_pulse.^2) ); % B1rms = 5.67;


%% Combine!
TR = 2;
B1rmsTR = sqrt( InvRMS.^2 .*invDur./TR + ...
                B1rms90.^2 .*(3.072/1000)./TR + ...
                B1rms180.^2 .*(3.000/1000)./TR); % 0.5951 uT


