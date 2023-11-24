% Temp notes on sim with adiabatic:

% Definte the adiabatic pulse:

Params.PulseOpt.beta = 672;  
Params.PulseOpt.n = 1; 
Params.PulseOpt.mu = 5;
Params.PulseOpt.A0 = 17.5;
Params.PulseOpt.nSamples = 512;
dispFig = 0;
Trf = 10.24/1000;
shape = 'hsn';


[rf_pulse, PulseOpt] = GetAdiabaticPulse( Trf, shape, dispFig, Params.PulseOpt);


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














