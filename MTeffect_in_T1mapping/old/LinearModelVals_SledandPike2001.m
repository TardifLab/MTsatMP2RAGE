function [M0b, kf, Ra, T2a, T2b, kr] = LinearModelVals_SledandPike2001( R1obs)
% use Table 3 in Sled and Pike 2001 (Magnetic Resonance in Medicine) to get
% linear models to estimate the other tissue parameters from input R1obs.
% The following values were used to make the models:

% R1 = [ 0.99, 1.72, 0.95, 1.7];
% M0b = [0.056, 0.152, 0.072, 0.161];
% kf = [2.2, 4.6, 2.4, 4.3];
% Ra = [0.99, 1.8, 0.93, 1.8];
% T2a = [55, 31, 56, 37]*10^(-3);
% T2b = [9.7, 11.8, 11.1, 12.3]*10^(-6);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% M0b = 0.124 * R1obs - 0.05594;
M0b = 0.06732 * R1obs.^1.569;
kf = 2.898 * R1obs - 0.5081;
Ra = 1.135* R1obs - 0.1412;
T2a = -0.02915* R1obs + 0.08381;
T2b = 2.166e-6* R1obs + 8.323e-6;

kr = kf./M0b;


