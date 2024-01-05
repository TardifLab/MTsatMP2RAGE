function [rf_pulse, Params] = GetAdiabaticPulse( Trf, shape, dispFig, Params)

% This need to take in parameter related to the adiabatic pulse, and return
% the B1. 
% Returns the complex RF pulse, and the parameters.
% dispFig is a flag to return the BlochSimulation results. 

%GetAdiabaticPulse Generate an RF pulse structure.
%   Pulse = GetAdiabaticPulse(alpha, delta, Trf, shape, PulseOpt)
%
%   A note on adiabatic pulses:
%   B1(t) = A(t) * exp( -1i *integral(omega1(t')) dt' )
%   where A(t) is the envelope, omega1 is the frequency sweep
%
%   --args--
%   alpha: Flip angle (in degrees).
%   delta: Off-resonance frequency (in Hz);
%   Trf: RF pulse duration (in seconds)
%   shape: String. Represents the shape of the RF envelope, and each have
%          their own unique functions handles associated to them.
%
%           -values-
%           'hard'
%           'sinc'
%           'sinchann'
%           'gaussian'
%           'gausshann'
%           'sincgauss'
%
%   --optional args--
%   PulseOpt: Struct. Contains optional parameters for pulse shapes. See
%             pulse shape objective function files for more information.
%
%   See also VIEWPULSE.
%


if (nargin < 4)
    Params.PulseOpt = struct;
    Params.M0a = 1;
    Params.M0b = 0.1;
    Params.Ra = 1;
end

switch shape
    % case 'hard';      pulse_fcn = @hard_pulse;  
    % case 'sinc';      pulse_fcn = @sinc_pulse;        
    % case 'sinchann';  pulse_fcn = @sinchann_pulse;        
    % case 'sincgauss'; pulse_fcn = @sincgauss_pulse;        
    % case 'gaussian';  pulse_fcn = @gaussian_pulse;        
    % case 'gausshann'; pulse_fcn = @gausshann_pulse;    
    % case 'fermi';     pulse_fcn = @fermi_pulse;
    case 'hsn'       
        [rf_pulse, Params] = hyperbolicSecant_pulse( Trf, Params, dispFig);
end







% b1     =  @(t) pulse_fcn(t,Trf,PulseOpt);
% if moxunit_util_platform_is_octave
%     amp    =  2*pi*alpha / ( 360 * gamma * quad(@(t) (b1(t)), 0, Trf) );
% else
%     amp    =  2*pi*alpha / ( 360 * gamma * integral(@(t) (b1(t)), 0, Trf) );
% end
% % amp    =  2*pi*alpha / ( 360 * gamma * integral(@(t) abs(b1(t)), 0, Trf,'ArrayValued',true) );
% omega  =  @(t) (gamma*amp*pulse_fcn(t,Trf,PulseOpt));
% omega2 =  @(t) (gamma*amp*pulse_fcn(t,Trf,PulseOpt)).^2;
% 
% Pulse.pulse_fcn = pulse_fcn;  % Fcn handle to pulse shape function
% Pulse.b1     =   b1;          % Fcn handle to pulse envelope amplitude
% Pulse.amp    =   amp;         % Pulse max amplitude
% Pulse.omega  =   omega;       % Fcn handle to pulse omega1
% Pulse.omega2 =   omega2;      % Fcn handle to pulse omega1^2 (power)
% Pulse.alpha  =   alpha;       % Flip angle
% Pulse.delta  =   delta;       % Pulse offset
% Pulse.Trf    =   Trf;         % Pulse duration
% Pulse.shape  =   shape;       % Pulse shape string
% Pulse.opt    =   PulseOpt;    % Additional options (e.g. TBW for sinc time-bandwidth window, bw for gaussian bandwidth)


end
