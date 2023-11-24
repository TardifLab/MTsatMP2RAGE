function setupSimPaths_MTsatMP2RAGE

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
qMRlabDirectory = 'E:\GitHub\qMRLab';

%% MP2RAGE scripts:
% https://github.com/JosePMarques/MP2RAGE-related-scripts
MP2RAGEDirectory = 'E:\GitHub\MP2RAGE-related-scripts';

%% Misc Code for a variety of things:
% https://github.com/christopherrowley/NeuroImagingMatlab
miscDir =  'E:\GitHub\NeuroImagingMatlab';



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
addpath(genpath(baseDir))
addpath(genpath(simCodeDir))
addpath(genpath(qMRlabDirectory)) 
addpath(genpath(MP2RAGEDirectory))
addpath(genpath(miscDir))

 
cd(baseDir); % makes saving things easier!


