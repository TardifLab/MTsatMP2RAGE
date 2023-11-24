function [M0b, kf, Ra, T2a, T2b, kr] = qMT_LinearModelVals( R1obs)

% Goal is to provide a rough estimate of how other parameters change with 
% R1obs, since there is a high spatial correlation between these metricsto get
% linear models to estimate the other tissue parameters from input R1obs.
% The following values were used to make the models:

% Goal is to extract a bunch of parameters and build some models

% % Parameters from Sled and Pike 2001:
% R1 = [ 0.99, 1.72, 0.95, 1.7];
% M0b = [0.056, 0.152, 0.072, 0.161];
% kf = [2.2, 4.6, 2.4, 4.3];
% Ra = [0.99, 1.8, 0.93, 1.8];
% T2a = [55, 31, 56, 37]*10^(-3);
% T2b = [9.7, 11.8, 11.1, 12.3]*10^(-6);

% % From Garcia et al 2015: - cancer study- ref ROI in non-affected WM
% % Table 1
% R1 = [1/0.806, 1/0.873, 1/0.839];
% kf = [3.3, 3.1, 3.1];
% M0b = [0.133, 0.122, 0.134];

% % Soustelle et al 2022 - qMT with dipolar contributions - average for WM
% % and DGM - last section
% R1 = [1/1.070, 1/1.458];
% M0b = [0.155, 0.094];

% % Jang et al 2023 - B1 corrected T1 and qMT
% % Taken from table 2 with first order correction
% R1 = [1/1.043, 1/1.064, 1/1.1121, 1/1.088, 1/1.522, 1/1.372, 1/1.407];
% M0b = [0.153, 0.149, 0.139, 0.1489, 0.0643, 0.0658, 0.088];
% kf = [2.758, 2.57, 2.196, 1.956, 1.322, 1.072, 1.358];

% % Mossahebi et al 2014 - Methods from fit to WM and GM
% R1 = [0.97, 0.71];
% M0b = [0.1536, 0.088];
% kf = [2.7, 1.57];
% T2b = [9.84, 10.21]*10^(-6);

% % Wood et al 2020 - qMT with SSFP ellipse
% % Values taken from end of discussion for bilateral WM ROI
% R1 = 1/0.794;
% M0b = 0.113;
% kf =5; % they mentioned a large number of fitting errors, maybe not include.
% T2a = 54*10^(-3);

% % Dortch et al 2011 qMT, extracted from Figure 4 (approximates
% R1= [1.04, 0.68];
% M0b = [0.114, 0.075];
% kf = [1.254, 1.125];

% % Bagnato et al 2019; SIR inMS
% % Taken from table 3, NWM
% R1 = 0.89;
% M0b = 0.125;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
M0b =  0.0681* R1obs+    0.0483;
kf = 3.38 * R1obs - 0.934;
Ra = R1obs;
T2a = -0.03* R1obs + 0.0863;
T2b = 2.12e-6* R1obs + 8.34e-6;

kr = kf./M0b;


