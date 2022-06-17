%% MPT model of study/test MST with old-similar-new responses

clear; close all;
preLoad = true;
printFigures = true;
dataDir = 'data/';
figureDir = 'figures/';
storageDir = '../storage/';

%% Data
dataName = 'msttDataWithin';
load([dataDir dataName], 'd');

%% Constants
load pantoneColors pantone;

%% MCMC sampling from graphical model

% which engine to use
engine = 'jags';

% graphical model script
modelName = 'mptStudyTestOSN_3';

% parameters to monitor
params = {...
    'muRho', 'muPsi', 'muDelta', 'muGammaO', 'muGammaS', ...
    'rho', 'psi', 'delta', 'gammaO', 'gammaS', 'gammaN'...
    'rhoRep', 'psiRep', 'deltaRep', 'gammaORep', 'gammaSRep', 'gammaNRep', ...
    'z', 'phi', ...
    'ppO', 'ppS', ...
    };

% MCMC properties
nChains    = 8;     % number of MCMC chains
nBurnin    = 5e3;   % number of discarded burn-in samples
nSamples   = 1e3;   % number of collected samples
nThin      = 10;     % number of samples between those collected
doParallel = 1;     % whether MATLAB parallel toolbox parallizes chains

% assign MATLAB variables to the observed nodes
data = struct(...
    'y'              , d.decisionOSN                  , ...
    'truth'          , d.truthOSN                          , ...
    'lureBin'        , d.lureBinOSN                       , ...
    'nParticipants'  , d.nParticipants                 , ...
    'p'              , d.participantOSN                   , ...
    'nLures'         , d.nLures             , ...
    'nTotalTrials'   , length(d.trialOSN)                 );

% generator for initialization
generator = @()struct('muRho', rand);

%% Sample using Trinity
fileName = sprintf('%s_%s_%s.mat', modelName, dataName, engine);

if preLoad && isfile(sprintf('%s/%s', storageDir, fileName))
    fprintf('Loading pre-stored samples for model %s on data %s\n', modelName, dataName);
    load(sprintf('%s/%s', storageDir, fileName), 'chains', 'stats', 'diagnostics', 'info');
else
   tic; % start clock
    [stats, chains, diagnostics, info] = callbayes(engine, ...
        'model'           , sprintf('%s_%s.txt', modelName, engine)   , ...
        'data'            , data                                      , ...
        'outputname'      , 'samples'                                 , ...
        'init'            , generator                                 , ...
        'datafilename'    , modelName                                 , ...
        'initfilename'    , modelName                                 , ...
        'scriptfilename'  , modelName                                 , ...
        'logfilename'     , sprintf('tmp/%s', modelName)              , ...
        'nchains'         , nChains                                   , ...
        'nburnin'         , nBurnin                                   , ...
        'nsamples'        , nSamples                                  , ...
        'monitorparams'   , params                                    , ...
        'thin'            , nThin                                     , ...
        'workingdir'      , sprintf('tmp/%s', modelName)              , ...
        'verbosity'       , 0                                         , ...
        'saveoutput'      , true                                      , ...
        'parallel'        , doParallel                                );
    fprintf('%s took %f seconds!\n', upper(engine), toc); % show timing

    fprintf('Saving samples for model %s on data %s\n', modelName, dataName);
    save(sprintf('%s/%s', storageDir, fileName), 'chains', 'stats', 'diagnostics', 'info', 'yPred');

    % convergence of each parameter
    disp('Convergence statistics:')
    grtable(chains, 1.05)

    % basic descriptive statistics
    disp('Descriptive statistics for all chains:')
    codatable(chains);

end

%% Analysis
rho = get_matrix_from_coda(chains, 'rho');
psi = get_matrix_from_coda(chains, 'psi');
delta = get_matrix_from_coda(chains, 'delta');

zMean = codatable(chains, 'z', @mean);
zMode = codatable(chains, 'z', @mode);
ppO = get_matrix_from_coda(chains, 'ppO');
ppS = get_matrix_from_coda(chains, 'ppS');

% user constants
markerSize = 8;
fontSize = 14;
CIbounds = [2.5 97.5];
whiskerWidth = 0.25;
modelColor = pantone.ClassicBlue;
dataColor = pantone.Tangerine;
dataSimColor = pantone.Greenery;
dataOldColor = pantone.DuskBlue;

binWidth = 0.01; violinScale = 5;
violinColor = pantone.ClassicBlue;

% derived
binsC = binWidth/2:binWidth:1-binWidth/2;
binsE = 0:binWidth:1;

% which participants in which order
keep = find(zMode == 1);
[~, sortIdx] = sort(d.participantCorrectOSN(keep), 'descend');
sortIdx = keep(sortIdx);
[nRows, nCols] = subplotArrange(length(sortIdx));

%% fit
F = figure; clf; hold on;
setFigure(F, [0.2 0.2 0.6 0.6], '');

for plotIdx = 1:length(sortIdx)

    part = sortIdx(plotIdx);

    subplot(nRows, nCols, plotIdx); hold on;
    set(gca, ...
        'xlim'       , [0.5 d.nLures+0.5]         , ...
        'xtick'      , 1:d.nLures  , ...
        'xticklabel'  , [], ...
        'XTickLabelRotation', 0, ...
        'ytick'      , 0:0.2:1 , ...
        'ylim'       , [0 1]         , ...
        'box'        , 'off'                 , ...
        'tickdir'    , 'out'                 , ...
        'layer'      , 'top'                 , ...
        'ticklength' , [0.025 0]              , ...
        'layer'      , 'top'                 , ...
        'fontsize'   , fontSize              );
    if plotIdx == 16
        set(gca, 'xticklabel', {'\delta_1', '\delta_2', '\delta_3', '\delta_4', '\delta_5'});
    else
        set(gca, 'yticklabel', []);
    end
    set(gca, 'pos', get(gca, 'pos') + [0 0.05 0 0])
    Raxes(gca, 0.01, 0.01);

    text(1, 1, sprintf('P%d', part), ...
        'fontsize', fontSize, ...
        'hor', 'cen');

    % accuracy for each stimulus type
    match = find(d.participantON == part);

    for idx = 1:d.nLures
        samples = chains.(sprintf('delta_%d_%d', part, idx))(:);
        count = histcounts(samples, binsE, 'normalization', 'probability');
        for i = 1:length(count)
            if count(i) > 0
                plot(idx + violinScale*count(i).*[-1 1], [binsC(i) binsC(i)], '-', ...
                    'color', violinColor, ...
                    'linewidth', 1);
            end
        end
    end
end

[~, AH(1)] = suplabel('Lure Discriminability');
[~, AH(2)] = suplabel('Posterior Probability', 'y');
set(AH, 'fontsize', fontSize + 4);

% print
if printFigures
    warning off;
    print(sprintf('%s/deltas_%s_%s.png', figureDir, dataName, modelName), '-dpng');
    print(sprintf('%s/deltas_%s_%s.eps', figureDir, dataName, modelName), '-depsc');
    warning on;
end



%% fit
 F = figure; clf; hold on;
setFigure(F, [0.2 0.2 0.6 0.6], '');

for plotIdx = 1:length(sortIdx)

    part = sortIdx(plotIdx);

    subplot(nRows, nCols, plotIdx); hold on;
    set(gca, ...
        'xlim'       , [0.5 d.nLures+2+0.5]         , ...
        'xtick'      , 1:d.nLures+2  , ...
        'xticklabel'  , [], ...
        'XTickLabelRotation', 0, ...
        'ytick'      , 0:0.2:1 , ...
        'ylim'       , [0 1]         , ...
        'box'        , 'off'                 , ...
        'tickdir'    , 'out'                 , ...
        'layer'      , 'top'                 , ...
        'ticklength' , [0.025 0]              , ...
        'layer'      , 'top'                 , ...
        'fontsize'   , fontSize              );
    if plotIdx == 16
        set(gca, 'xticklabel', {'O', 'L_1', 'L_2', 'L_3', 'L_4', 'L_5', 'N'});
    else
        set(gca, 'yticklabel', []);
    end
    set(gca, 'pos', get(gca, 'pos') + [0 0.05 0 0])
    Raxes(gca, 0.01, 0.01);

    text(1, 1, sprintf('P%d', part), ...
        'fontsize', fontSize, ...
        'hor', 'cen');

    % old responses for each stimulus type
    match = find(d.participantOSN == part);
    yO(1) = nansum(d.decisionOSN(match) == 1 & d.truthOSN(match) == 1)/nansum(d.truthOSN(match) == 1);
    for idx = 1:d.nLures
        yO(idx+1) = nansum(d.decisionOSN(match) == 1 & d.truthOSN(match) == 3 & d.lureBinOSN(match) == idx)/nansum(d.truthOSN(match) == 3 & d.lureBinOSN(match) == idx);
    end
    yO(d.nLures+2) = nansum(d.decisionOSN(match) == 1 & d.truthOSN(match) == 2)/nansum(d.truthOSN(match) == 2);

    H(1) = plot(1:d.nLures+2, yO, 'o-',...
        'markersize', 8, ...
        'color', dataOldColor, ...
        'markeredgecolor', 'w', ...
        'markerfacecolor', dataOldColor, ...
        'linewidth', 2);

    for idx = 1:d.nLures+2

        H(2) = plot(1:d.nLures+2, ppO(part,:), '+',...
            'markersize', 10, ...
            'color', modelColor, ...
            'linewidth', 2);
% 
%         bounds = prctile(chains.(sprintf('ppO_%d_%d', part, idx))(:), CIbounds);
%         mn = ppO(part, idx);
%         plot(idx+[-whiskerWidth +whiskerWidth], [bounds(1) bounds(1)], '-', ...
%             'color', modelColor, ...
%             'linewidth', 0.5);
%         plot(idx+[-whiskerWidth +whiskerWidth], [bounds(2) bounds(2)], '-', ...
%             'color', modelColor, ...
%             'linewidth', 0.5);
     end

    % similar responses for each stimulus type
    match = find(d.participantOSN == part);
    yS(1) = nansum(d.decisionOSN(match) == 3 & d.truthOSN(match) == 1)/nansum(d.truthOSN(match) == 1);
    for idx = 1:d.nLures
        yS(idx+1) = nansum(d.decisionOSN(match) == 3 & d.truthOSN(match) == 3 & d.lureBinOSN(match) == idx)/nansum(d.truthOSN(match) == 3 & d.lureBinOSN(match) == idx);
    end
    yS(d.nLures+2) = nansum(d.decisionOSN(match) == 3 & d.truthOSN(match) == 2)/nansum(d.truthOSN(match) == 2);

    H(3) = plot(1:d.nLures+2, yS, 'o:',...
        'markersize', 8, ...
        'color', dataSimColor, ...
        'markeredgecolor', 'w', ...
        'markerfacecolor', dataSimColor, ...
        'linewidth', 2);

    for idx = 1:d.nLures+2

        H(4) = plot(1:d.nLures+2, ppS(part,:), 'x',...
            'markersize', 10, ...
            'color', modelColor, ...
            'linewidth', 2.5);

%         bounds = prctile(chains.(sprintf('ppS_%d_%d', part, idx))(:), CIbounds);
%         mn = ppS(part, idx);
%         plot(idx+[-whiskerWidth +whiskerWidth], [bounds(1) bounds(1)], ':', ...
%             'color', modelColor, ...
%             'linewidth', 0.5);
%         plot(idx+[-whiskerWidth +whiskerWidth], [bounds(2) bounds(2)], ':', ...
%             'color', modelColor, ...
%             'linewidth', 0.5);
     end

 %   return

end

L =  legend(H, 'old data', 'old model', 'sim data', 'sim model', ...
    'box', 'off', ...
    'location', 'northeast', ...
    'fontsize', fontSize);
set(L, 'pos', get(L, 'pos') + [0.175 0 0 0]);

[~, AH(1)] = suplabel('Stimulus Type');
[~, AH(2)] = suplabel('Probability Old/Similar', 'y');
set(AH, 'fontsize', fontSize + 4);

% print
if printFigures
    warning off;
    print(sprintf('%s/fit_%s_%s.png', figureDir, dataName, modelName), '-dpng');
    print(sprintf('%s/fit_%s_%s.eps', figureDir, dataName, modelName), '-depsc');
    warning on;
end
