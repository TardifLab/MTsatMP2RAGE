function Params = CR_getSeqParams_MNIprot( seqStr)

%filled with values from simulation results:
Params.B0 = 3;
Params.MTC = 1; % Magnetization Transfer Contrast
Params.TissueType = 'GM';
Params = DefaultCortexTissueParams(Params);
Params = CalcImagingParams(Params);

Params.WExcDur = 0.1/1000; % duration of water pulse
Params.PerfectSpoiling = 1;

%% MT params
if  strcmp( seqStr, 'csMP2RAGE_MTsat') || strcmp( seqStr, 'MP2RAGE_MTsat') || ...
        strcmp( seqStr, 'VFA_MTsat')

    Params.delta = 2200; % Saturation pulse offset in Hz
    Params.flipAngle = 5; % excitation flip angle water.
    Params.TR = 25/1000; % total repetition time = MT pulse train and readout.
    Params.numSatPulse = 1;
    Params.TurboFactor = 1;
    Params.pulseDur = 12.8/1000; %duration of 1 MT pulse in seconds
    Params.satFlipAngle = 540; % degrees
    Params.satTrainPerBoost = 1;
    Params.TR_MT = 0; 
    Params.freqPattern = 'single'; % options: 'single', 'dualAlternate', 'dualContinuous'
    Params.pulseGapDur = 0.3/1000; %ms gap between MT pulses in train % C.R. new, shift from 1ms to 0.5
    Params.DummyEcho = 0;
    Params.numExcitation = 1; % number of readout lines/TR
    Params.boosted = 0;
    Params.satTrainPerBoost = 1; 
    Params.TR_MT = 0; 
    Params.echospacing = 7.66/1000;
    
elseif strcmp( seqStr, 'meGRE_MTsat')
    Params.delta = 2000; % Saturation pulse offset in Hz
    Params.flipAngle = 6; % excitation flip angle water.
    Params.TR = 23/1000; % total repetition time = MT pulse train and readout.
    Params.numSatPulse = 1;
    Params.TurboFactor = 1;
    Params.pulseDur = 4.096/1000; %duration of 1 MT pulse in seconds
    Params.satFlipAngle = 220; % degrees
    Params.satTrainPerBoost = 1;
    Params.TR_MT = 0; 
    Params.freqPattern = 'single'; % options: 'single', 'dualAlternate', 'dualContinuous'
    Params.pulseGapDur = 0.3/1000; %ms gap between MT pulses in train % C.R. new, shift from 1ms to 0.5
    Params.DummyEcho = 0;
    Params.numExcitation = 1; % number of readout lines/TR
    Params.boosted = 0;
    Params.satTrainPerBoost = 1; 
    Params.TR_MT = 0; 
    Params.echospacing = 7.66/1000;

else
    error('Input Sequence String does not match the provided folder names. No Match Found')
end 

% MT parameters that will be consistent:
Params.SatPulseShape = 'gaussian';
Params.TD_MT =  0 ;   
Params = CalcVariableImagingParams(Params);



































