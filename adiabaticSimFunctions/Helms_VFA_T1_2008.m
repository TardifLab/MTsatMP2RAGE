function T1 = Helms_VFA_T1_2008( PD_sig, T1w_sig, PD_fa, T1w_fa, PD_TR, T1w_TR, B1)

% Input:
%   PD and T1w are structs with the fields:
%     {PD,T1w}_sig   signal values
%     {PDw,T1w}_fa   (nominal flip angle, deg)
%     {PDw,T1w}_TR   (repetition time, s or ms)
%     B1   (array of B1 ratios: actual fa / nominal fa)

% Output will be in the units of the input TR.

% See original function by Helms group: 
% https://github.com/hMRI-group/hMRI-toolbox/blob/master/hmri_calc_R1.m


PD_fa = deg2rad(PD_fa) * B1;
T1w_fa = deg2rad(T1w_fa) * B1;

R1=0.5*( PD_sig.*PD_fa/PD_TR - T1w_sig.*T1w_fa/T1w_TR )...
    ./( T1w_sig./T1w_fa - PD_sig./PD_fa );

T1 = 1./R1;
T1(isnan(T1)) = 0; 





