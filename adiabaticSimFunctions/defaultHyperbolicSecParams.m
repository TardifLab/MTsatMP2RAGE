function PulseOpt = defaultHyperbolicSecParams(PulseOpt)

% Function designed to be used with the adiabatic pulse code
% Fills in default values if they are not user-defined

if(~isfield(PulseOpt,'beta') || isempty(PulseOpt.beta) || ~isfinite(PulseOpt.beta))
    % Default beta value in rad/s
    PulseOpt.beta = 672;       
end

if(~isfield(PulseOpt,'n') || isempty(PulseOpt.n) || ~isfinite(PulseOpt.n))
    % Default sech exponent
    PulseOpt.n = 1;       
end

if(~isfield(PulseOpt,'mu') || isempty(PulseOpt.mu) || ~isfinite(PulseOpt.mu))
    % Default phase modulation parameter
    PulseOpt.mu = 5;       
end

if(~isfield(PulseOpt,'A0') || isempty(PulseOpt.A0) || ~isfinite(PulseOpt.A0))
    % Peak B1 of the pulse in microTesla
    PulseOpt.A0 = 17.5;       
end

if(~isfield(PulseOpt,'nSamples') || isempty(PulseOpt.nSamples) || ~isfinite(PulseOpt.nSamples))
    % Peak B1 of the pulse in microTesla
    PulseOpt.nSamples = 512;       
end

return;