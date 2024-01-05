function M_return = blochSimAdiabaticPulse_1pool( rf_pulse, Trf, delta,...
                                            Params, M_start, B)
%% Sim Adiabatic Pulse

% 'rf_pulse' is a 1xnSamples vector that stores the B1 in microtesla over
%            time
% 'Trf' is the pulse duration in seconds
% 'PulseOpt' contains pulse details. Only PulseOpt.nSamples is needed here
% 'delta' is used to see the offset frequency to do a frequency sweep 
%  'M_start' is the magnetization vector just before the start of the pulse
%            (3x1)
%  'B' is the thermal equilibrium magnetization vector (3x1)
%
% Written by Christopher Rowley 2023
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    if ~exist('B','var') || isempty(B)
        B = [0, 0, 0, 0, Params.M0a, Params.M0b]';    
    end
    
    nSamples = Params.PulseOpt.nSamples;
       
    %% If you wanted to do single pool:
    Mt = zeros(3, nSamples+1);
    Mt(:,1) = M_start; % start mag = [0 0 1];
    I = eye(3); % identity matrix      
    
    % We will numerically evaluate this over time 
    dt = Trf/nSamples;
    R2a = Params.R2a; %1000/80; % 80 ms
    R1a = Params.Ra; % 1; % 1000 ms
    
    for t = 1:nSamples
    
        w1 = 2*pi *42.577478518 * rf_pulse(t);  % assume B1 is in microTesla, and drop the 10^6 from gamma. w1 in rad/s
        % Generate RF matrix
        A_rf =[ -R2a,   -2*pi*delta, -imag(w1); ...       % Water X
                2*pi*delta,    -R2a, -real(w1);...        % Water Y
                imag(w1),  real(w1), -R1a];   %  Water Z
        % Apply
        AExp = expm(A_rf*dt);
        AEnd = (AExp - I)* (A_rf\B);
        Mt(:,t+1) = pagemtimes(AExp, Mt(:,t)) + AEnd;
    
    end
    M_return = Mt(:,end);
    figure; plot(0:dt:Trf, Mt(3,:), 'LineWidth',3);
    hold on
    plot(0:dt:Trf, Mt(1,:), 'LineWidth',3);
    plot(0:dt:Trf, Mt(2,:), 'LineWidth',3);
    legend('M_z', 'M_x', 'M_y');

return;
