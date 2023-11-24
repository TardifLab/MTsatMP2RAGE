%% Test extra parameter options for T1 mapping
% Here we have included 3 other options for VFA and MP2RAGE

%% VFA -
% Baudrexel et al 2018
    % Image 1 - PDw
        Params.TR = 16.4/1000;
        Params.flipAngle = 4;
    % Image 2 - T1w
        Params.TR = 16.4/1000;
        Params.flipAngle = 24;

% Weiskopf et al 2013 - MPM
    % Image 1 - PDw
        Params.TR = 23.7/1000;
        Params.flipAngle = 6;
    % Image 2 - T1w
        Params.TR = 18.7/1000;
        Params.flipAngle = 20;

% Helms et al 2008 - Quantitative FLASH MRI
    % Image 1 - PDw
        Params.TR = 30/1000;
        Params.flipAngle = 7;
    % Image 2 - T1w
        Params.TR = 30/1000;
        Params.flipAngle = 20;

%% MP2RAGE

% Marques et al 2010 -3T
    Params.TR = 6.25;
    Params.flipAngle = [7,5];
    Params.numExcitation = 160;
    Params.echoSpacing = 5.8/1000;
    Params.Readout = 'linear';
    Params.TI = [800, 2200]./1000;
    Params.DummyEcho = 0;
    Params.PerfectSpoiling = 1;


% Marques et al 2010 -3T short TR
    Params.TR = 5;
    Params.flipAngle = [7,5];
    Params.numExcitation = 160;
    Params.echoSpacing = 5.8/1000;
    Params.Readout = 'linear';
    Params.TI = [700, 2500]./1000;
    Params.DummyEcho = 0;
    Params.PerfectSpoiling = 1;


% Mussard et al 2020 -fully sampled data
    Params.TR = 5;
    Params.flipAngle = [4,5];
    Params.numExcitation = 195;
    Params.echoSpacing = 7.7/1000; % not provided. Using from current study
    Params.Readout = 'linear';
    Params.TI = [700, 2500]./1000;
    Params.DummyEcho = 0;
    Params.PerfectSpoiling = 1;








