% Sim MTsat SNR
SNR = 70; 

% see the function in Batch_sim/functions/CR_add_ihMTsat_SNR.m
Params.B0 = 3;
Params.MTC = 1; % Magnetization Transfer Contrast
Params.TissueType = 'GM';
Params.numExcitation = 1;
Params.echoSpacing = 7.7/1000;

Params = DefaultCortexTissueParams(Params);
Params = CalcImagingParams(Params);
Params.PerfectSpoiling = 1;



Params.delta = 2000;
Params.flipAngle = 6; % excitation flip angle water.
Params.TR = 27/1000; % total repetition time = MT pulse train and readout.
Params.numSatPulse = 1;
Params.TurboFactor = 1;
Params.pulseDur = 4/1000; %duration of 1 MT pulse in seconds
Params.satFlipAngle = 220; % microTesla
Params.satTrainPerBoost = 1;
Params.TR_MT = 0;    
Params.freqPattern = 'single'; % options: 'single', 'dualAlternate', 'dualContinuous'
Params.pulseGapDur = 0.3/1000; %ms gap between MT pulses in train % C.R. new, shift from 1ms to 0.5
Params.DummyEcho = 0;
Params.numExcitation = Params.TurboFactor + Params.DummyEcho; % number of readout lines/TR
Params.boosted = 0;
Params.SatPulseShape = 'gaussian'; % options: 'hanning', 'gaussian', 'square' % CR adjustes
Params.PulseOpt.bw = 2./Params.pulseDur;

Params.TD_MT =  Params.TR_MT - Params.numSatPulse* (Params.pulseDur + Params.pulseGapDur) ;   
Params = CalcVariableImagingParams(Params);




simMT_signal = BlochSimFlashSequence_v2(Params);  


noiseLvl = simMT_signal/SNR;

T1 = 1./ Params.Raobs;


MP2RAGE.B0=3;           % in Tesla
MP2RAGE.TR=5;           % MP2RAGE TR in seconds 
MP2RAGE.TRFLASH=6.4e-3; % TR of the GRE readout
MP2RAGE.TIs=[940e-3 2830e-3];% inversion times - time between middle of refocusing pulse and excitatoin of the k-space center encoding
MP2RAGE.NZslices=176;% C.R. just need the turbofactor here.
MP2RAGE.FlipDegrees=[4 5];% Flip angle of the two readouts in degrees



% Generate noise vectors (30,000 x 4)
additive_noise = normrnd( 1,noiseLvl ,30000,5) - 1;
%std(additive_noise) % Below 30k, you don't get the same value

% can check visually the result:
%histogram(additive_noise(:,1))

% Use the T1 and M0 == 1 to get simulated signal values for MP2RAGE
MP2RAGE_Signal = 1*MPRAGEfunc(2, MP2RAGE.TR, MP2RAGE.TIs, MP2RAGE.NZslices, MP2RAGE.TRFLASH, MP2RAGE.FlipDegrees, 'normal',T1);

% Make 2 vectors, 1 for each inversion time, and add the random noise
inv1.img = additive_noise(:,1) + MP2RAGE_Signal(1);
inv2.img = additive_noise(:,2) + MP2RAGE_Signal(2);

% Generate Uni Image
Uni.img=real(inv1.img).*conj(inv2.img)./(abs(inv1.img).^2+abs(inv2.img).^2); % taken from line 68 of MP2RAGE_lookuptable.m
B1.img = ones(size(inv2.img));

% uncenter and multiply up
Uni.img = double(Uni.img+0.5)*4095;

% Calculate noisy vectors for M0 and T1
[ T1_n, ~, M0_n] = CR_T1B1correctpackageTFL_withM0( B1, Uni, inv2, MP2RAGE, [], 1);

T1_n = T1_n.img; % want T1 in ms for lookup table
M0_n = M0_n.img;

% View results if desired
% histogram(M0_n)
% histogram(T1_n)

%% For MTsat

noise1 = additive_noise(:,3);

% Generate noisy signals
Single_n = simMT_signal + noise1;

    
MTsat_sim_S = calcMTsatThruLookupTablewithDummyV3( Single_n,   [], T1_n, [], M0_n,...
    Params.echoSpacing * 1000, Params.numExcitation, Params.TR * 1000, ...
    Params.flipAngle, Params.DummyEcho);
    


MT_noise = std(MTsat_sim_S);
MTsat_SNR = mean(MTsat_sim_S)./ MT_noise;


T1_SNR = mean(T1_n)./ std(T1_n);
M0_SNR = mean(M0_n)./ std(M0_n);



