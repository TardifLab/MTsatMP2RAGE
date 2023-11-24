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

%% MT-weighted image
TR = 27/1000;
flipAngle = 6;
wExcDur = 0.1/1000; % ms

excB1rms = getExcPulseB1rms(flipAngle, wExcDur);

% Sat pulse
B1rmsSat = 3.26;
tSat = 12/1000;
B1rmsTR = sqrt( B1rmsSat.^2 .*tSat./TR + excB1rms.^2 .*wExcDur./TR);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% MP2RAGE
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
invPulse = hyperbolicSecant_pulse( invDur);

% rms
dt = invDur./length(invPulse);
InvRMS = sqrt( (1/invDur) * trapz( linspace(0, invDur, length(invPulse)), invPulse.^2) );

% Combine!
B1rmsTR = sqrt( InvRMS.^2 .*invDur./TR + ...
                B1rmsBlock1.^2 .*t1./TR + ...
                B1rmsBlock2.^2 .*t2./TR);

% Apparently 0.31 uT?

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% cs - MP2RAGE
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
invPulse = hyperbolicSecant_pulse( invDur, [], 1);

% rms
dt = invDur./length(invPulse);
InvRMS = sqrt( (1/invDur) * trapz( linspace(0, invDur, length(invPulse)), invPulse.^2) );

% Combine!
B1rmsTR = sqrt( InvRMS.^2 .*invDur./TR + ...
                B1rmsBlock1.^2 .*t1./TR + ...
                B1rmsBlock2.^2 .*t2./TR);










% %% Spin Echo 180 then 90...
% 
% % Block 3 - inversion pulse:
% invPulse = hyperbolicSecant_pulse( invDur);
% 
% % rms
% dt = invDur./length(invPulse);
% InvRMS = sqrt( (1/invDur) * trapz( linspace(0, invDur, length(invPulse)), invPulse.^2) );
% 
% % Combine!
% B1rmsTR = sqrt( InvRMS.^2 .*invDur./TR + ...
%                 B1rmsBlock1.^2 .*t1./TR + ...
%                 B1rmsBlock2.^2 .*t2./TR);


