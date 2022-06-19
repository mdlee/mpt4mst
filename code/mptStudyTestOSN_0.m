%% MPT model of study/test with OSN
% basic model

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
modelName = 'mptStudyTestOSN_0';

% parameters to monitor
params = {...
    'rho', 'psi', 'delta', 'gammaO', 'gammaS', ...
    'ppO', 'ppS'};

% MCMC properties
nChains    = 8;     % number of MCMC chains
nBurnin    = 1e3;   % number of discarded burn-in samples
nSamples   = 1e3;   % number of collected samples
nThin      = 1;     % number of samples between those collected
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
    save(sprintf('%s/%s', storageDir, fileName), 'chains', 'stats', 'diagnostics', 'info');

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
gammaO = get_matrix_from_coda(chains, 'gammaO');
gammaS = get_matrix_from_coda(chains, 'gammaS');
delta = get_matrix_from_coda(chains, 'delta');
ppO = get_matrix_from_coda(chains, 'ppO');
ppS = get_matrix_from_coda(chains, 'ppS');

% user constants
markerSize = 8;
fontSize = 14;
CIbounds = [25 75];
whiskerWidth = 0.25;
modelColor = pantone.ClassicBlue;
dataColor = pantone.Tangerine;
dataSimColor = pantone.DuskBlue;
dataOldColor = pantone.Tangerine;

binWidth = 0.02; violinScale = 5;
violinColor = pantone.ClassicBlue;

% derived
binsC = binWidth/2:binWidth:1-binWidth/2;
binsE = 0:binWidth:1;

% which participants in which order
keep = 1:d.nParticipants;
[~, sortIdx] = sort(d.participantCorrectON(keep), 'descend');
sortIdx = keep(sortIdx);
[nRows, nCols] = subplotArrange(length(sortIdx));

%% parameters

% joint posterior figure
F = figure; clf; hold on;
setFigure(F, [0.2 0.2 0.5*1.1 0.6*1.1], '');
set(gcf, 'renderer', 'opengl');

subplot(221); hold on;
set(gca, ...
    'xlim'       , [0 1]               , ...
    'xtick'      , 0:0.2:1          , ...
    'ylim'       , [0 1]         , ...
    'ytick'      , 0:0.2:1 , ...
    'box'        , 'off'                     , ...
    'tickdir'    , 'out'                     , ...
    'xticklabelrotation', 0, ...
    'layer'      , 'top'                     , ...
    'ticklength' , [0.02 0]                  , ...
    'layer'      , 'top'                     , ...
    'fontsize'   , fontSize                  );
axis square;
xlabel('Memory (\rho)', 'fontsize', fontSize+8);
H = ylabel('Memory Absent (\psi)', 'fontsize', fontSize+8);
set(H, 'pos', get(H, 'pos') + [-0.06 0 0]);
set(gca, 'pos', get(gca, 'pos') + [0 0.025 0 0]);
Raxes(gca, 0.015, 0.01);

% draw posterior densities
for idx = 1:d.nParticipants
    rhoBounds = prctile(chains.(sprintf('rho_%d', idx))(:), CIbounds);
    tauBounds = prctile(chains.(sprintf('psi_%d', idx))(:), CIbounds);
    
    plot([rho(idx) rho(idx)], tauBounds, '-', ...
        'color', pantone.GlacierGray);
    plot(rhoBounds, [psi(idx) psi(idx)], '-', ...
        'color', pantone.GlacierGray);
end

for idx = 1:d.nParticipants
    H = plot(rho(idx), psi(idx), 'o', ...
        'markerfacecolor', pantone.ClassicBlue, ...
        'markeredgecolor', 'w');
    text(rho(idx), psi(idx), sprintf(' %s', char(64+find(sortIdx == idx))), 'fontsize', fontSize);
end

subplot(222); hold on;
set(gca, ...
    'xlim'       , [0 1]               , ...
    'xtick'      , 0:0.2:1          , ...
    'ylim'       , [0 1]         , ...
    'ytick'      , 0:0.2:1 , ...
    'box'        , 'off'                     , ...
    'tickdir'    , 'out'                     , ...
    'xticklabelrotation', 0, ...
    'layer'      , 'top'                     , ...
    'ticklength' , [0.02 0]                  , ...
    'layer'      , 'top'                     , ...
    'fontsize'   , fontSize                  );
axis square;
xlabel('Memory (\rho)', 'fontsize', fontSize+8);
ylabel('Guess Old (\gamma_o)', 'fontsize', fontSize+8);
set(gca, 'pos', get(gca, 'pos') + [0 0.025 0 0]);
Raxes(gca, 0.015, 0.01);

% draw posterior densities
for idx = 1:d.nParticipants
    rhoBounds = prctile(chains.(sprintf('rho_%d', idx))(:), CIbounds);
    tauBounds = prctile(chains.(sprintf('gammaO_%d', idx))(:), CIbounds);
    
    plot([rho(idx) rho(idx)], tauBounds, '-', ...
        'color', pantone.GlacierGray);
    plot(rhoBounds, [gammaO(idx) gammaO(idx)], '-', ...
        'color', pantone.GlacierGray);
end

for idx = 1:d.nParticipants
    H = plot(rho(idx), gammaO(idx), 'o', ...
        'markerfacecolor', pantone.ClassicBlue, ...
        'markeredgecolor', 'w');
    text(rho(idx), gammaO(idx), sprintf(' %s', char(64+find(sortIdx == idx))), 'fontsize', fontSize);
    
end

subplot(223); hold on;
set(gca, ...
    'xlim'       , [0 1]               , ...
    'xtick'      , 0:0.2:1          , ...
    'ylim'       , [0 1]         , ...
    'ytick'      , 0:0.2:1 , ...
    'box'        , 'off'                     , ...
    'tickdir'    , 'out'                     , ...
     'xticklabelrotation', 0, ...
   'layer'      , 'top'                     , ...
    'ticklength' , [0.02 0]                  , ...
    'layer'      , 'top'                     , ...
    'fontsize'   , fontSize                  );
axis square;
xlabel('Memory (\rho)', 'fontsize', fontSize+8);
ylabel('Guess Similar (\gamma_s)', 'fontsize', fontSize+8);
set(gca, 'pos', get(gca, 'pos') + [0 0.025 0 0]);
Raxes(gca, 0.015, 0.01);

% draw posterior densities
for idx = 1:d.nParticipants
    rhoBounds = prctile(chains.(sprintf('rho_%d', idx))(:), CIbounds);
    tauBounds = prctile(chains.(sprintf('gammaS_%d', idx))(:), CIbounds);
    
    plot([rho(idx) rho(idx)], tauBounds, '-', ...
        'color', pantone.GlacierGray);
    plot(rhoBounds, [gammaS(idx) gammaS(idx)], '-', ...
        'color', pantone.GlacierGray);
end

for idx = 1:d.nParticipants
    H = plot(rho(idx), gammaS(idx), 'o', ...
        'markerfacecolor', pantone.ClassicBlue, ...
        'markeredgecolor', 'w');
    text(rho(idx), gammaS(idx), sprintf(' %s', char(64+find(sortIdx == idx))), 'fontsize', fontSize);
    
end

subplot(224); hold on;
set(gca, ...
    'xlim'       , [0 1]               , ...
    'xtick'      , 0:0.2:1          , ...
    'ylim'       , [0 1]         , ...
    'ytick'      , 0:0.2:1 , ...
    'box'        , 'off'                     , ...
    'tickdir'    , 'out'                     , ...
     'xticklabelrotation', 0, ...
   'layer'      , 'top'                     , ...
    'ticklength' , [0.02 0]                  , ...
    'layer'      , 'top'                     , ...
    'fontsize'   , fontSize                  );
axis square;
xlabel('Guess Old (\gamma_o)', 'fontsize', fontSize+8);
ylabel('Guess Similar (\gamma_s)', 'fontsize', fontSize+8);
set(gca, 'pos', get(gca, 'pos') + [0 0.025 0 0]);
Raxes(gca, 0.015, 0.01);

% draw posterior densities
for idx = 1:d.nParticipants
    rhoBounds = prctile(chains.(sprintf('gammaO_%d', idx))(:), CIbounds);
    tauBounds = prctile(chains.(sprintf('gammaS_%d', idx))(:), CIbounds);
    
    plot([gammaO(idx) gammaO(idx)], tauBounds, '-', ...
        'color', pantone.GlacierGray);
    plot(rhoBounds, [gammaS(idx) gammaS(idx)], '-', ...
        'color', pantone.GlacierGray);
end

for idx = 1:d.nParticipants
    H = plot(gammaO(idx), gammaS(idx), 'o', ...
        'markerfacecolor', pantone.ClassicBlue, ...
        'markeredgecolor', 'w');
    text(gammaO(idx), gammaS(idx), sprintf(' %s', char(64+find(sortIdx == idx))), 'fontsize', fontSize);
    
end

% print
if printFigures
    warning off;
    print(sprintf('%s/parameters_%s_%s.png', figureDir, dataName, modelName), '-dpng');
    print(sprintf('%s/parameters_%s_%s.eps', figureDir, dataName, modelName), '-depsc');
    warning on;
end

%% delta
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
    if plotIdx == (nRows*nCols-nCols+1)
        set(gca, 'xticklabel', {'\delta_1', '\delta_2', '\delta_3', '\delta_4', '\delta_5'});
    else
        set(gca, 'yticklabel', []);
    end
    set(gca, 'pos', get(gca, 'pos') + [0 0.025 0 0])
    Raxes(gca, 0.01, 0.01);

 text(1, 1.1, sprintf('%s', char(64+plotIdx)), ...
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
                    'linewidth', 2);
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
setFigure(F, [0.1 0.2 0.8 0.6], '');

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
    if plotIdx == (nRows*nCols - nCols+1)
        set(gca, 'xticklabel', {'O', 'L_1', 'L_2', 'L_3', 'L_4', 'L_5', 'N'});
    else
        set(gca, 'yticklabel', []);
    end
    set(gca, 'pos', get(gca, 'pos') + [-0.025 0.025 0 0])
    Raxes(gca, 0.01, 0.01);

 text(1, 1.1, sprintf('%s', char(64+plotIdx)), ...
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
    end

end

L =  legend(H, 'old data', 'old model', 'sim data', 'sim model', ...
    'box', 'off', ...
    'location', 'northeast', ...
    'fontsize', fontSize);
set(L, 'pos', get(L, 'pos') + [0.125 0 0 0]);

[~, AH(1)] = suplabel('Stimulus Type');
set(AH(1), 'pos', get(AH(1), 'pos') + [0 0.025 0]);
[~, AH(2)] = suplabel('Probability Old/Similar', 'y');
set(AH(2), 'pos', get(AH(2), 'pos') + [0.025 0 0]);
set(AH, 'fontsize', fontSize + 12);

% print
if printFigures
    warning off;
    print(sprintf('%s/fit_%s_%s.png', figureDir, dataName, modelName), '-dpng');
    print(sprintf('%s/fit_%s_%s.eps', figureDir, dataName, modelName), '-depsc');
    warning on;
end
