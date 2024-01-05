function [rf_pulse, Params] = hyperbolicSecant_pulse( Trf, Params, dispFigure)

%   hyperbolicSecant_pulse Adiabatic hyperbolic secant RF pulse function.
%   pulse = hyperbolicSecant_pulse(t, Trf, PulseOpt)
%
%   B1(t) = A(t) * exp( -1i *integral(omega1(t')) dt' )
%   where A(t) is the envelope, omega1 is the frequency sweep
%
%   For the case of a hyperbolic secant pulse:
%   A(t) = A0 * sech(Beta*t)
%   omega1(t) = -mu*Beta*tanh(Beta*t)
%   A0 is the peak amplitude in microTesla
%   Beta is a frequency modulation parameter in rad/s
%   mu is a phase modulation parameter (dimensionless)
%
%   The pulse is defined to be 0 outside the pulse window (before 
%   t = 0 or after t=Trf). (HSn, n = 1-8+) 
%
%   --args--
%   t: Function handle variable, represents the time.
%   Trf: Duration of the RF pulse in seconds.
%
%   --optional args--
%   PulseOpt: Struct. Contains optional parameters for pulse shapes.
%   PulseOpt.Beta: frequency modulation parameter
%   PulseOpt.n: time modulation - Typical 4 for non-selective, 1 for slab
% 
%   Reference: Matt A. Bernstein, Kevin F. Kink and Xiaohong Joe Zhou.
%              Handbook of MRI Pulse Sequences, pp. 110, Eq. 4.10, (2004)
%
%              Tannús, A. and M. Garwood (1997). "Adiabatic pulses." 
%              NMR in Biomedicine 10(8): 423-434.
%
%
% To be used with qMRlab
% Written by Christopher Rowley 2023


if ~exist('dispFigure','var') || isempty(dispFigure) || ~isfinite(dispFigure)
    dispFigure = 0;      
end


% Function to fill default values;
Params.PulseOpt = defaultHyperbolicSecParams(Params.PulseOpt);

nSamples = Params.PulseOpt.nSamples;  
t = linspace(0, Trf, nSamples);

% Amplitude
A_t =  Params.PulseOpt.A0* sech(Params.PulseOpt.beta* ( (t - Trf/2)).^Params.PulseOpt.n);
A_t((t < 0 | t>Trf)) = 0;
% disp( ['Average B1 of the pulse is:', num2str(mean(A_t))]) 


% Frequency modulation function 
% Carrier frequency modulation function w(t):
omega1 = -Params.PulseOpt.mu.*Params.PulseOpt.beta .* ...
            tanh(Params.PulseOpt.beta .* (t - Trf/2))./(2*pi); % 2pi to convert from rad/s to Hz

% Phase modulation function phi(t):
phi = Params.PulseOpt.mu .* log(sech(Params.PulseOpt.beta .* (t - Trf/2)) );

% Put together complex RF pulse waveform:
rf_pulse = A_t .* exp(1i .* phi);

%% Can do Bloch Sim to get inversion profile and display figure if interested:

if dispFigure
    M_start = [0, 0, 0, 0, Params.M0a, Params.M0b]';
    b1Rel = 0.5:0.1:1.5;
    freqOff = -2000:200:2000;
    [b1m, freqm] = ndgrid(b1Rel, freqOff);
    
    Mza = zeros(size(b1m));
    Mzb = zeros(size(b1m));
    
    for i = 1:length(b1Rel)
        for j = 1:length(freqOff)
        
            M_return = blochSimAdiabaticPulse( b1Rel(i)*rf_pulse, Params.Inv,  ...
                            freqOff(j), Params, M_start, []);

            Mza(i,j) = M_return(5);
            Mzb(i,j) = M_return(6);
        end
    end

    figure; tiledlayout(2,2)
    nexttile; plot(t*1000, A_t, 'LineWidth', 3); 
    xlabel('Time(ms)'); ylabel('B_1 (μT)')
    title('Amplitude Function');ax = gca; ax.FontSize = 20;
    
    nexttile; plot(t*1000, omega1, 'LineWidth', 3);
    xlabel('Time(ms)'); ylabel('Frequency (Hz)');
    title('Frequency Modulation function');ax = gca; ax.FontSize = 20;
    
    nexttile; surf(b1m, freqm, Mza);
    xlabel('Rel. B1'); ylabel('Freq (Hz)'); zlabel('M_{za}');ax = gca; ax.FontSize = 20;
    
    nexttile; surf(b1m, freqm, Mzb);
    xlabel('Rel. B1'); ylabel('Freq (Hz)'); zlabel('M_{zb}');ax = gca; ax.FontSize = 20;
    
    set(gcf,'Position',[100 100 1200 1000])
end

return; 








































