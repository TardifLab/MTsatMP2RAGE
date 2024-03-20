function [outSig1, M, time_vect] = BlochSim_SpinEcho_IR_Sequence(Params, varargin)

%% Overview of how the code works with applicable function calls:

% Sequence timing:
% TI1 = time between 180 inversion pulse, and middle of first excitation train(TF/2)
% ET1 = evolution time, between 180 inversion and first excitation
% EBT = excitation block timing (Turbofactor * echospacing) 
% TD = time delay after last excitation pulse (and echospacing) until next
%   inversion
% TR = ET + EBT + TD

% We will assume that we will do full Bloch simulation of all pulses, so
% that the time steps will have the pulse durations subtracted

%% Need to define defaults for:
% Params.Inv.Trf
% Params.InversionEfficiency
% Params.flipAngle, Params.flipAngle

Params.CalcVector = 1;

%% Use name-value pairs to override other variables set. Great for parfor loops!
for i = 1:2:length(varargin)
    if ischar(varargin{i})
        Params.(varargin{i}) = varargin{i+1};
    end
end

if ~isfield(Params,'kf')
    Params.kf = (Params.R*Params.M0b); 
end
if ~isfield(Params,'kr')
    Params.kr = (Params.R*Params.M0a);
end

if isempty(Params.Ra) % allow you to specify either Ra or Raobs
    Params.Ra = Params.Raobs - ((Params.R * Params.M0b * (Params.R1b - Params.Raobs)) / (Params.R1b - Params.Raobs + Params.R));
    if isnan(Params.Ra)
        Params.Ra = 1;
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Build sequence, then convert to loop structure.
% play sequence for 5 seconds, then fill sampling table.

num2avgOver = 10; % you can get some variation in signal, so keep the last few and average

% Equivalent of 5 seconds of imaging to steady state, then record data. 
loops = ceil(6/Params.TR) + num2avgOver;

%% Standard Stuff
M0 = [0 0 0 0 Params.M0a, Params.M0b]'; 

%% Timing Variables
% INV ... Time1...EXC...Time2...Refocus...Time3...Acquire...TD... REPEAT

% Inversion time (TI) is the time between the inversion time and the
% excitation pulse
% TE/2 is the time between the excitation pulse and the refocus pulse
% TE/2 is the time between the refocus pulse and acquisition.

% Inversion time to calculate relaxation. Center on the pulses
Time1 = Params.TI - Params.Inv.Trf/2 - Params.Exc.Trf/2 ; 

% Evolution between Excitation and refocusing, and refocus and acquisiton
Time2 = Params.TE/2 - Params.Exc.Trf/2 - Params.refocus.Trf/2;

% Evolution between refocus and acquisiton
Time3 = Params.TE/2 - Params.refocus.Trf/2;
 
% Relaxation after acquisition: (need to account for half the inversion
% pulse
TD = Params.TR - Params.TI - Params.TE-Params.Inv.Trf/2 ;

if Params.TR ~= (Params.Inv.Trf + Time1 + Params.ExcPulseDur + Time2+...
        Params.RefocusPulseDur + Time3+ TD)
    error('sequence timing does not add up to TR value')
end

%% Calculate the B1 waveforms of the RF Pulses:

% Inversion
[inv_pulse, ~] = GetAdiabaticPulse( Params.Inv.Trf, Params.Inv.shape, 0, Params.Inv);

% 90
Pulse = GetPulse(Params.Exc.alpha, 0, Params.Exc.Trf, Params.Exc.shape, Params.Exc.PulseOpt);
t = 0:Params.Exc.Trf/(Params.Exc.nSamples-1):Params.Exc.Trf;
exc_pulse = 1000*Pulse.amp*Pulse.b1(t);

% Refocus
Pulse = GetPulse(Params.refocus.alpha, 0, Params.refocus.Trf, Params.refocus.shape, Params.refocus.PulseOpt);
t = 0:Params.refocus.Trf/(Params.refocus.nSamples-1):Params.refocus.Trf;
refocus_pulse = 1000*Pulse.amp*Pulse.b1(t);


%% Need to determine a sufficient number of isochromats. 
% can use this function to use more if needed.
% Params.N_spin = DetermineNumberIsoChromat(Params, TD)

% if Params.PerfectSpoiling % number of spins wont matter in this case
%     Params.N_spin = 1;
% else
%     Params.N_spin = 201;
% end
if ~Params.PerfectSpoiling
    error('Isochromate simulation currently not implemented here')
end
Params.N_spin = 1;


%% Setup Matrices
if Params.CalcVector == 1
    M = zeros(6,loops*20); % time-dependent magnetization storage for plotting
    M(:,1) = M0;
    time_vect = zeros( loops*20,1);
    idx = 2;
end

M_t = repmat(M0,1, Params.N_spin); % vector for a single time point

Sig_vec1 = zeros(num2avgOver, 1 );
rep = 1; % to count over the number to average over

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Start of sequence loop
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for i = 1:loops
    
    %% Start by applying an inversion pulse

    % bloch sim and return magnetization
    M_t = blochSimAdiabaticPulse( inv_pulse, Params.Inv ,...
                                  0, Params, M_t, []);

    M_t(1:4,:) = 0; % Typically gradients used to crush this

    if Params.CalcVector == 1
        M(:,idx) = mean(M_t,2); 
        time_vect(idx) = time_vect(idx-1)+ Params.Inv.Trf;
        idx = idx+1;
    end % For viewing; 

    %% Spin evolution over time time 1
    M_t = XYmag_Spoil_6vec(Params, M_t, Time1, 0, 0);

    if Params.CalcVector == 1
        M(:,idx) = mean(M_t,2); 
        time_vect(idx) = time_vect(idx-1) + Time1;
        idx = idx+1;
    end % For viewing; 


    %% Excitation Pulse
    % bloch sim and return magnetization
    M_t = blochSimAdiabaticPulse( exc_pulse, Params.Exc ,...
                                0, Params, M_t, []);

    if Params.CalcVector == 1
        M(:,idx) = mean(M_t,2); 
        time_vect(idx) = time_vect(idx-1) + Params.Exc.Trf;
        idx = idx+1;
    end % For viewing;    
    
    %% Spin evolution over time time 2
    M_t = XYmag_Spoil_6vec(Params, M_t, Time2, 0, 0);

    if Params.CalcVector == 1
        M(:,idx) = mean(M_t,2); 
        time_vect(idx) = time_vect(idx-1) + Time2;
        idx = idx+1;
    end % For viewing; 
   
    
    %% Refocus Pulse
    % bloch sim and return magnetization
    M_t = blochSimAdiabaticPulse( refocus_pulse, Params.refocus ,...
                                0, Params, M_t, []);

    if Params.CalcVector == 1
        M(:,idx) = mean(M_t,2); 
        time_vect(idx) = time_vect(idx-1) + Params.refocus.Trf ;
        idx = idx+1;
    end % For viewing;    
    
    %% Spin evolution over time time3
    M_t = XYmag_Spoil_6vec(Params, M_t, Time3, 0, 0);

    if Params.CalcVector == 1
        M(:,idx) = mean(M_t,2); 
        time_vect(idx) = time_vect(idx-1) + Time3;
        idx = idx+1;
    end % For viewing; 
   
  
    %% Store the magnetization of each excitation pulse after 5 seconds prep
    if i > loops-num2avgOver

        Sig_vec1(rep) = TransverseMagnetizationMagnitude_6vec(M_t);

       if (i == loops) % if simulation is done...             
           outSig1 = mean(Sig_vec1); 
       end
        rep = rep+1;
    end % End 'if SS_reached'   

    %% Spin evolution over time TD
    M_t = XYmag_Spoil_6vec(Params, M_t, TD, 0, 0);

    if Params.CalcVector == 1
        M(:,idx) = mean(M_t,2); 
        time_vect(idx) = time_vect(idx-1) + TD;
        idx = idx+1;
    end % For viewing;    
    
end

M(:,idx:end) = [];
time_vect(idx:end) = [];

%% Debug and view
% figure;
% plot(time_vect, sqrt(sum(M([1,3],:).^2)))
% % 
% figure;
% plot(time_vect, M(5,:),'LineWidth',3); xlim([0 6])
% hold on; plot(time_vect, sqrt(sum(M([1,3],:).^2)),'LineWidth',3)


% warning('') % Clear last warning message
% [warnMsg, warnId] = lastwarn;
% if ~isempty(warnMsg)
% return
% end


