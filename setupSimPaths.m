function setupSimPaths

% Edit this function to contain the paths of your folders and run at the
% start of scripts. 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Add This Repository Code base directory:
% https://github.com/TardifLab/MTsatMP2RAGE
baseDir = 'E:\GitHub\MTsatMP2RAGE\';
 
%% Add in Simulation Code:
% https://github.com/TardifLab/OptimizeIHMTimaging
simCodeDir = 'E:\GitHub\OptimizeIHMTimaging\';

%% qMRlab code directory:
% http://qmrlab.org/
% https://github.com/qMRLab/qMRLab
qMRlabDirectory = 'E:\GitHub\qMRLab-master';

%% Misc Code for a variety of things:
% https://github.com/christopherrowley/NeuroImagingMatlab
% The MP2RAGE script are all included here as well!
miscDir =  'E:\GitHub\NeuroImagingMatlab';



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
addpath(genpath(baseDir))
addpath(genpath(simCodeDir))
addpath(genpath(qMRlabDirectory)) 
addpath(genpath(miscDir))

 



