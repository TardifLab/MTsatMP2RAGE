function M_out = SpinEvolution_Relaxation_6vec( Params, M_in, t)

% This version has 6 vectors. XYZ for free and bound. In order of XX,YY,ZZ
% Intended to spoil the magnetization at the end of a magnetizaton
% saturation pulse train.
% Calculate spin evolution in the absence of RF and gradients

% Input is:
% Params structure that stores a bunch of variables
% M_in is a 6x Number_spin vector
% t is the time over which this step occurs

B = [0 0 0 0 Params.Ra*Params.M0a Params.R1b*Params.M0b]';
I = eye(6);

% Tissue Parameters:
R2a = 1/Params.T2a; % 80 ms
R2b = 1/Params.T2b; % 12 us
R1a = Params.Ra; % 1000 ms
R1b = Params.R1b; % 1000 ms
kr = Params.kr;
kf = Params.kf;

% Evolution Matrix
w1 = 0;  % no external RF field
delta = 0;
E =[ -(R2a+kf),          kr, -2*pi*delta,          0, -imag(w1),        0; ...       % Water X
               kf,   -(R2b+kr),          0, -2*pi*delta,       0,  -imag(w1); ...       % Bound X
      2*pi*delta,           0,  -(R2a+kf),         kr, -real(w1),        0;...        % Water Y
                0, 2*pi*delta,         kf,  -(R2b+kr),       0,  -real(w1);...        % Bound Y
        imag(w1),           0,  real(w1),          0, -(R1a+kf),    kr;...        % Water Z
                0,   imag(w1),          0,  real(w1),      kf,  -(R1b+kr)];        % Bound Z 


% Apply
M_out = expm(E*t) * M_in + (expm(E*t) - I)* (E\B);



