function [outSig1, outSig2, M, time_vect] = BlochSim_MP2RAGESequence_6vec(Params, varargin)

% This version removes dipolar, and add XY for bound pool 
%% Overview of how the code works with applicable function calls:

% Sequence timing:
% TI1 = time between 180 inversion pulse, and middle of first excitation train(TF/2)
% ET1 = evolution time, between 180 inversion and first excitation
% EBT = excitation block timing (Turbofactor * echospacing) 
% TI2 = time between 180 inversion pulse, and middle of second excitation train(TF/2)
% ET2 = evolution time, between first and second excitation blocks
% TD = time delay after last excitation pulse (and echospacing) until next
%   inversion
% TR = ET + EBT + TD

% if centric, dummy echoes == 2, else ==0. 

%% Need to define defaults for:
% Params.InvPulseDur
% Params.InversionEfficiency
% Params.flipAngle, Params.flipAngle

Params.CalcVector = 1;

%% Use name-value pairs to override other variables set. Great for parfor loops!
for i = 1:2:length(varargin)
    if ischar(varargin{i})
        Params.(varargin{i}) = varargin{i+1};
    end
end

if length(Params.flipAngle) < 2
    Params.flipAngle = [Params.flipAngle, Params.flipAngle];
    disp('Only one flip angle entered, assuming it is used for both readouts');
end

if ~isfield(Params,'kf')
    Params.kf = (Params.R*Params.M0b); 
end
if ~isfield(Params,'kr')
    Params.kr = (Params.R*Params.M0a);
end

if ~isfield(Params,'Ra') || isempty(Params.Ra) % allow you to specify either Ra or Raobs
    Params.Ra = Params.Raobs - ((Params.R * Params.M0b * (Params.R1b - Params.Raobs)) / (Params.R1b - Params.Raobs + Params.R));
    if isnan(Params.Ra)
        Params.Ra = 1;
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Build sequence, then convert to loop structure.
% play sequence for 5 seconds, then fill sampling table.

num2avgOver = 3; % you get some variation in signal, so keep the last few and average

% Equivalent of 5 seconds of imaging to steady state, then record data. 
loops = ceil(6/Params.TR) + num2avgOver;

%% Standard Stuff
M0 = [0 0 0 0 Params.M0a, Params.M0b]'; 

if Params.echoSpacing == 0 && Params.numExcitation > 1
    error( 'Please define Params.echoSpacing');
end   


% if centric, dummy echoes == 2, else ==0. 
if strcmp(Params.Readout, 'centric')
    error('Currently only supporting linear readouts.')
end

% With the above, we are assuming the signal value is in the middle of the
% readout:
readNum = ceil(Params.numExcitation/2);

%% Timing Variables
% INV ... ET1... EBT... ET2... EBT...TD... REPEAT

% excitation block timing (Turbofactor * echospacing) 
EBT = (Params.numExcitation+Params.DummyEcho)*( Params.echoSpacing);

% Inversion time needs to be specified by user
TI1 = Params.TI(1); 
TI2 = Params.TI(2); 

% evolution time, between 180 inversion and first excitation
% Incorporate the length of the pulses.  
ET1 = TI1 - EBT/2  - Params.Inv.Trf/2; % take off half the inversion pulse
ET2 = (TI2 - EBT/2) - (TI1 + EBT/2); 

if ET2 < 0 || ET1 < 0
    error('not enough time to permit acquisition blocks.')
end

% time delay after last excitation pulse (and echospacing) until next inversion
TD = Params.TR - (TI2 + EBT/2) - Params.Inv.Trf/2; % take off half the inversion pulse

if Params.TR ~= (ET1 + EBT + ET2 + EBT+ TD + Params.Inv.Trf)
    error('sequence timing does not add up to TR value')
end

%% If simulating adiabatic pulse, get the pulse:
[inv_pulse, ~] = GetAdiabaticPulse( Params.Inv.Trf,  Params.Inv.shape, 0, Params.Inv);

% Exc1
Pulse = GetPulse(Params.flipAngle(1), 0, Params.Exc.Trf, Params.Exc.shape);
t = 0:Params.Exc.Trf/(Params.Exc.nSamples-1):Params.Exc.Trf;
exc_pulse1 = -1*1000*Pulse.amp*Pulse.b1(t);
% The negative 1 factor is to apply in right direction based on how I wrote
% out other bloch equations

% Exc2
Pulse = GetPulse(Params.flipAngle(2), 0, Params.Exc.Trf, Params.Exc.shape);
t = 0:Params.Exc.Trf/(Params.Exc.nSamples-1):Params.Exc.Trf;
exc_pulse2 = -1*1000*Pulse.amp*Pulse.b1(t);

%% Need to determine a sufficient number of isochromats. 
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
Sig_vec2 = zeros(num2avgOver, 1 );
rep = 1; % to count over the number to average over

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Start of sequence loop
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for i = 1:loops    
    %% Start by applying an inversion pulse
    M_t = blochSimAdiabaticPulse( inv_pulse, Params.Inv ,...
                                    0, Params, M_t, []);

    M_t(1:4,:) = 0; % Typically gradients used to crush this

    if Params.CalcVector == 1
        M(:,idx) = mean(M_t,2); 
        time_vect(idx) = time_vect(idx-1)+ Params.Inv.Trf;
        idx = idx+1;
    end % For viewing; 

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% Spin evolution over time ET1
    M_t = XYmag_Spoil_6vec(Params, M_t, ET1, 0, 0);

    if Params.CalcVector == 1
        M(:,idx) = mean(M_t,2); 
        time_vect(idx) = time_vect(idx-1) + ET1;
        idx = idx+1;
    end % For viewing; 

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% Excitation Block  
    for j = 1: Params.numExcitation
        
        % Calculate rotation matrix for excitation-specific phase
        M_t = blochSimAdiabaticPulse( exc_pulse1, Params.Exc ,...
                                0, Params, M_t, []);
    
        if Params.CalcVector == 1
            M(:,idx) = mean(M_t,2); 
            time_vect(idx) = time_vect(idx-1)+ Params.Exc.Trf;
            idx = idx+1;
        end % For viewing;  

        %% Store the magnetization of each excitation pulse after 5 seconds prep
        if i > loops-num2avgOver && (j == readNum) 

           Sig_vec1(rep) = TransverseMagnetizationMagnitude_6vec(M_t);
    
           if (i == loops) && (j == readNum) % if simulation is done...             
               outSig1 = mean(Sig_vec1); 
           end
        end % End 'if SS_reached'   
        
        %% Evolution and Apply Spoiling
        M_t = XYmag_Spoil_6vec( Params, M_t, Params.echoSpacing- Params.Exc.Trf, 0, 0);
        M_t(1:4,:) = 0;

        if Params.CalcVector == 1
            M(:,idx) = mean(M_t,2); 
            time_vect(idx) = time_vect(idx-1)+ Params.echoSpacing- Params.Exc.Trf;
            idx = idx+1;
        end % For viewing;      
    end % End '1: Params.numExcitation' 

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% Spin evolution between two excitation trains, over time ET2
    M_t = XYmag_Spoil_6vec(Params, M_t, ET2, 0, 0);

    if Params.CalcVector == 1
        M(:,idx) = mean(M_t,2); 
        time_vect(idx) = time_vect(idx-1) + ET2;
        idx = idx+1;
    end % For viewing; 

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% Excitation Block - 2nd image
    for j = 1: Params.numExcitation        
        % Calculate rotation matrix for excitation-specific phase
        M_t = blochSimAdiabaticPulse( exc_pulse2, Params.Exc ,...
                                0, Params, M_t, []);
    
        if Params.CalcVector == 1
            M(:,idx) = mean(M_t,2); 
            time_vect(idx) = time_vect(idx-1) + Params.Exc.Trf;
            idx = idx+1;
        end % For viewing;  

        %% Store the magnetization of each excitation pulse after 5 seconds prep
        if i > loops-num2avgOver && (j == readNum) 

            Sig_vec2( rep ) = TransverseMagnetizationMagnitude_6vec(M_t);
    
           if (i == loops) && (j == readNum) % if simulation is done...             
               outSig2 = mean(Sig_vec2); % output 1xTurbofactor vector
               
               if Params.CalcVector == 1
                   M(:,idx:end) = [];
                   time_vect(idx:end) = [];
               end             
               return;
           end
            
           % increase repetition index
           if j == readNum
               rep = rep+1;
           end
        end % End 'if SS_reached'   

        %% Evolution and Apply Spoiling
        M_t = XYmag_Spoil_6vec( Params, M_t, Params.echoSpacing- Params.Exc.Trf, 0, 0);
        M_t(1:4,:) = 0;
            
        if Params.CalcVector == 1
            M(:,idx) = mean(M_t,2); 
            time_vect(idx) = time_vect(idx-1)+ Params.echoSpacing- Params.Exc.Trf;
            idx = idx+1;
        end % For viewing;      
        
    end % End '1: Params.numExcitation' 

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% Spin evolution over time TD
    M_t = XYmag_Spoil_6vec(Params, M_t, TD, 0, 0);

    if Params.CalcVector == 1
        M(:,idx) = mean(M_t,2); 
        time_vect(idx) = time_vect(idx-1)+TD;
        idx = idx+1;
    end % For viewing; 

end


% MP2RAGE.B0          = Params.B0;                  % In Tesla
% MP2RAGE.TR          = Params.TR;                  % MP2RAGE TR in seconds
% MP2RAGE.TRFLASH     = Params.echoSpacing;             % TR of the GRE readout
% MP2RAGE.TIs         = Params.TI;   % Inversion times - time between middle of refocusing pulse and excitatoin of the k-space center encoding
% MP2RAGE.NZslices    = Params.numExcitation;            % Slices Per Slab * [PartialFourierInSlice-0.5  0.5]
% MP2RAGE.FlipDegrees = Params.flipAngle;  
% B1.img = 1;
% MP2RAGEimg.img = calculate_UNI_from_sims( outSig1, outSig2);
% MP2RAGEINV2img.img = outSig2;
% [ T1corr, ~, ~] = CR_T1B1correctpackageTFL_withM0( B1, MP2RAGEimg, MP2RAGEINV2img, MP2RAGE, [], 0.96);
% mp2rage_T1 = T1corr.img;
% disp(['Input T1: ',num2str(1/Params.Raobs*1000),'; Output T1: ',num2str(mp2rage_T1)]);

%% Debug and view
% figure;
% plot(time_vect, sqrt(sum(M([1,3],:).^2))); xlim([0 15])
% % 
% figure;
% plot(time_vect, M(5,:)); xlim([0 15])

% Sig= CR_FLASH_solver(5, Params.echoSpacing, 1, 1, 1./Params.Raobs)


% figure;
% plot(time_vect, M(4,:))

% warning('') % Clear last warning message
% [warnMsg, warnId] = lastwarn;
% if ~isempty(warnMsg)
% return
% end

% 
% figure;
% plot(M_t(1,:),'-r','LineWidth',1)
% hold on
% plot(M_t(2,:),'--b')

