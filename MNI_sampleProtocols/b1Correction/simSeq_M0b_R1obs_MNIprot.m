%% Simulate sequence and generate fitequation to cover the spectrum of MTsat
% results for varying B1 field, R1obs and M0b. 
% Please consult README document first to be sure you have downloaded all
% necessary packages.

% Written by Christopher Rowley 2023
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Make sure the necessary paths are set and change to the correct working
% directory
setupSimPaths_MTsatMP2RAGE;

OutputDir =  'MNI_sampleProtocols/simOutputs';

seqStr = {'meGRE_MTsat', 'csMP2RAGE_MTsat', 'MP2RAGE_MTsat','VFA_MTsat'};

b1 = linspace(0, 1.5, 30); % relative B1 field
M0b = 0:0.03:0.18; 
T1obs = horzcat(0.6:0.06:2,2.1:0.4:4); %600ms to 4500ms to cover WM to CSF. 
Raobs = 1./T1obs;

for i = 3:length(seqStr)

    Params = CR_getSeqParams_MNIprot( seqStr{i} );

    % Loop variables:
    Params.M0b =  []; % going to loop over this
    Params.Raobs = [];
    Params.Ra = [];
    GRE_sigs = zeros(size(b1,2),size(M0b,2),size(Raobs,2));

    tic
    for z = 1:size(b1,2) % took nearly 5 hours for matrix 25x41x33.
            
        for j = 1:size(M0b,2)
            Params.M0b = M0b(j);
            
            for k = 1:size(Raobs,2)
                Params.Raobs = Raobs(k);
                GRE_sigs(z,j,k) = BlochSimFlashSequence_v2(Params, 'satFlipAngle' , Params.satFlipAngle*b1(z));       
            end
        end

        disp(z/size(b1,2) *100)  % print percent done...
        toc
    end


    %% MTsat calculation
    %reformat Aapp and R1app matrices for 3D calculation
    Aapp = ones(size(GRE_sigs));
    T1app = repmat(T1obs,[7,1,size(b1,2)]);
    T1app = permute(T1app,[3,1,2]);
    
    MTsat_sim_S = calcMTsatThruLookupTablewithDummyV3( GRE_sigs,   [], T1app* 1000, [], Aapp,...
        Params.echoSpacing * 1000, Params.numExcitation, Params.TR * 1000, Params.flipAngle, Params.DummyEcho);

    MTsatValue_fn = fullfile(OutputDir, strcat('MTsat_sim_',seqStr{i},'.mat')); 
    save(MTsatValue_fn,'MTsat_sim_S')

    %% Clean up then fit:
    MTsat_sim_S(MTsat_sim_S < 0) = NaN;
 
    % Single
    [fit_SS_eqn, fit_SS_eqn_sprintf, fit_SSsat, numTerms] = CR_generateFitequationsV2(M0b(1:6), b1, Raobs, MTsat_sim_S(:,1:6,:));

    % put into one variable for export
    fitValues.fitvals_coeff = fit_SSsat.Coefficients;
    fitValues.fit_SS_eqn = fit_SS_eqn;
    fitValues.fit_SS_eqn_sprintf = fit_SS_eqn_sprintf;
    fitValues.Params = Params; % export params to reference later if desired
    fitValues.numTerms = numTerms; % for fitting later...

    fitValue_fn = fullfile(OutputDir, strcat('fitValues_',seqStr{i},'.mat')); 
    save(fitValue_fn,'fitValues')

    img_fn = fullfile(OutputDir, strcat('simFig_',seqStr{i},'.png')); 
    CR_generateFitSimFiguresRelB1(M0b(1:6), b1, Raobs, MTsat_sim_S(:,1:6,:), fit_SS_eqn, img_fn)

    str = ['Done simulations for: ', seqStr{i}];
    disp(str)
    toc
end







