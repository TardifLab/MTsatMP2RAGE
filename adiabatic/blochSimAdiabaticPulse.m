function M_return = blochSimAdiabaticPulse( rf_pulse, PulseParams, delta,...
                                            Params, M_start, B)
%% Sim Adiabatic Pulse

% 'rf_pulse' is a 1xnSamples vector that stores the B1 in microtesla over
%            time
% PulseParams structure containing:
%   Trf -> is the pulse duration in seconds
%   nSamples -> number of samples in the pulse. Typically 512 (or multiple
%               of 256)
% 'Params' stores the tissue parameters for the simluation
% 'delta' is used to see the offset frequency to do a frequency sweep 
%  'M_start' is the magnetization vector just before the start of the pulse
%            (6x1)
%  'B' is the thermal equilibrium magnetization vector (6x1)
%
% Written by Christopher Rowley 2023
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    if ~exist('B','var') || isempty(B)
        B = [0, 0, 0, 0, Params.Ra*Params.M0a, Params.R1b*Params.M0b]';    
    end
    
    if ~isfield(Params,'kf')
        Params.kf = (Params.R*Params.M0b); 
    end
    if ~isfield(Params,'kr')
        Params.kr = (Params.R*Params.M0a);
    end

    nSamples = PulseParams.nSamples;
    
    % Follow equations in Murase 2011 - Pull values from Params struct
    R2a = 1/Params.T2a; % 80 ms
    R2b = 1/Params.T2b; % 12 us
    R1a = Params.Ra; % 1000 ms
    R1b = Params.R1b; % 1000 ms
    kr = Params.kr;
    kf = Params.kf;
    
    Mt = zeros(6, nSamples+1);
    I = eye(6); % identity matrix      
    
    Mt(:,1) = M_start; % start mag
    
    dt = PulseParams.Trf/nSamples;
    
    for t = 1:nSamples
    
        w1 = 2*pi *42.577478518 * rf_pulse(t);  % assume B1 is in microTesla, and drop the 10^6 from gamma. w1 in rad/s
    
        A_rf =[ -(R2a+kf),          kr, -2*pi*delta,          0, -imag(w1),        0; ...       % Water X
                       kf,   -(R2b+kr),          0, -2*pi*delta,       0,  -imag(w1); ...       % Bound X
              2*pi*delta,           0,  -(R2a+kf),         kr, -real(w1),        0;...        % Water Y
                        0, 2*pi*delta,         kf,  -(R2b+kr),       0,  -real(w1);...        % Bound Y
                imag(w1),           0,  real(w1),          0, -(R1a+kf),    kr;...        % Water Z
                        0,   imag(w1),          0,  real(w1),      kf,  -(R1b+kr)];        % Bound Z 
    
    
        % Apply
        Mt(:,t+1) = expm(A_rf*dt) * Mt(:,t) + (expm(A_rf*dt) - I)* (A_rf\B);

    end
    
    M_return = Mt(:,end);
return;

% figure; tiledlayout(1,2); nexttile;
% plot(0:dt:Trf, Mt(5,:), 'LineWidth',3); 
% xlabel('Time(s)'); ylabel('M_z Water'); ax = gca; ax.FontSize = 20;
% nexttile;
% plot(0:dt:Trf, Mt(6,:), 'LineWidth',3); 
% xlabel('Time(s)'); ylabel('M_z Bound'); ax = gca; ax.FontSize = 20;
% set(gcf,'Position',[100 100 800 500])
% 
% % Proportion saturation
% disp(['Relative Sat on MTpool from Inversion Pulse:', num2str(Mt(6,end)/Mt(6,1))]);



