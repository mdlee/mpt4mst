%% MPT model of study/test MST with old-new responses
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
modelName = 'mptStudyTestON_0';

% parameters to monitor
params = {...
    'rho', 'psi', 'delta', 'gamma', ...
    'pp'};

% MCMC properties
nChains    = 8;     % number of MCMC chains
nBurnin    = 1e3;   % number of discarded burn-in samples
nSamples   = 1e3;   % number of collected samples
nThin      = 1;     % number of samples between those collected
doParallel = 1;     % whether MATLAB parallel toolbox parallizes chains

% assign MATLAB variables to the observed nodes
truth = d.truthON;
truth(find(d.lureON == 1)) = 3;
data = struct(...
    'y'              , 2 - d.decisionON                  , ...
    'truth'          , truth                           , ...
    'lureBin'        , d.lureBinON                       , ...
    'nParticipants'  , d.nParticipants                 , ...
    'p'              , d.participantON                   , ...
    'nLures'         , d.nLures             , ...
    'nTotalTrials'   , length(d.trialON)                 );

% generator for initialization
generator = @()struct('muRho', rand);

%% Sample using Trinity
fileName = sprintf('%s_%s_%s.mat', modelName, dataName, engine);

if preLoad && isfile(sprintf('%s/%s', storageDir, fileName))
    fprintf('Loading pre-stored samples for model %s on data %s\n', modelName, dataName);
    load(sprintf('%s/%s', storageDir, fileName), 'chains', 'stats', 'diagnostics', 'info', 'yPred');
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
gamma = get_matrix_from_coda(chains, 'gamma');
delta = get_matrix_from_coda(chains, 'delta');
pp = get_matrix_from_coda(chains, 'pp');

% user constants
markerSize = 8;
fontSize = 14;
CIbounds = [25 75];
whiskerWidth = 0.25;
modelColor = pantone.ClassicBlue;
dataColor = pantone.Tangerine;

binWidth = 0.01; violinScale = 5;
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
setFigure(F, [0.1 0.2 0.5 0.5], '');
set(gcf, 'renderer', 'opengl');

AX1 = subplot(121); hold on;
set(gca, ...
    'xlim'       , [0 1]               , ...
    'xtick'      , 0:0.2:1          , ...
    'ylim'       , [0 1]         , ...
    'ytick'      , 0:0.2:1 , ...
    'box'        , 'off'                     , ...
    'tickdir'    , 'out'                     , ...
    'layer'      , 'top'                     , ...
    'ticklength' , [0.02 0]                  , ...
    'layer'      , 'top'                     , ...
    'fontsize'   , fontSize                  );
axis square;
xlabel('Memory (\rho)', 'fontsize', fontSize+8);
ylabel('Memory Absent (\psi)', 'fontsize', fontSize+8);
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

AX2 = subplot(122); hold on;
set(gca, ...
    'xlim'       , [0 1]               , ...
    'xtick'      , 0:0.2:1          , ...
    'ylim'       , [0 1]         , ...
    'ytick'      , 0:0.2:1 , ...
    'box'        , 'off'                     , ...
    'tickdir'    , 'out'                     , ...
    'layer'      , 'top'                     , ...
    'ticklength' , [0.02 0]                  , ...
    'layer'      , 'top'                     , ...
    'fontsize'   , fontSize                  );
axis square;
xlabel('Memory (\rho)', 'fontsize', fontSize+8);
ylabel('Guessing Old (\gamma)', 'fontsize', fontSize+8);
set(gca, 'pos', get(gca, 'pos') + [0.025 0.025 0 0]);
Raxes(gca, 0.015, 0.01);

% draw posterior densities
for idx = 1:d.nParticipants
    rhoBounds = prctile(chains.(sprintf('rho_%d', idx))(:), CIbounds);
    tauBounds = prctile(chains.(sprintf('gamma_%d', idx))(:), CIbounds);
    
    plot([rho(idx) rho(idx)], tauBounds, '-', ...
        'color', pantone.GlacierGray);
    plot(rhoBounds, [gamma(idx) gamma(idx)], '-', ...
        'color', pantone.GlacierGray);
end

for idx = 1:d.nParticipants
    H = plot(rho(idx), gamma(idx), 'o', ...
        'markerfacecolor', pantone.ClassicBlue, ...
        'markeredgecolor', 'w');
    text(rho(idx), gamma(idx), sprintf(' %s', char(64+find(sortIdx ==idx))), 'fontsize', fontSize);  
end

% print
if printFigures
    warning off;
    print(sprintf('%s/parameters_%s_%s.png', figureDir, dataName, modelName), '-dpng');
    print(sprintf('%s/parameters_%s_%s.eps', figureDir, dataName, modelName), '-depsc');
    warning on;
end

% annotate
axes(AX1);
rectangle('position', [0.025 0.025 0.15 0.15], ...
    'facecolor', [pantone.Greenery 0.5], ...
    'linewidth', 15, ...
    'edgecolor', [pantone.Custard 0.5]);

axes(AX2);
rectangle('position', [0.25 0.62 0.55 0.25], ...
    'facecolor', [pantone.ClassicBlue 0.5], ...
    'edgecolor', 'none');
rectangle('position', [0.45 0.15 0.6 0.45], ...
    'facecolor', [pantone.DuskBlue 0.5], ...
    'edgecolor', 'none');
rectangle('position', [0 0.2 0.2 0.3], ...
    'facecolor', [pantone.Greenery 0.5], ...
    'edgecolor', 'none')
rectangle('position', [0 0.8 0.2 0.3], ...
    'facecolor', [pantone.Custard 0.5], ...
    'edgecolor', 'none')

% print
if printFigures
    warning off;
    print(sprintf('%s/parametersAnnotate_%s_%s.png', figureDir, dataName, modelName), '-dpng');
    print(sprintf('%s/parametersAnnotate_%s_%s.eps', figureDir, dataName, modelName), '-depsc');
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
                    'linewidth', 1);
            end
        end
    end
end

[~, AH(1)] = suplabel('Lure Discriminability');
set(AH(1), 'pos', get(AH(1), 'pos') + [0 0.025 0]);
[~, AH(2)] = suplabel('Posterior Probability', 'y');
set(AH(2), 'pos', get(AH(2), 'pos') + [0.025 0 0]);

set(AH, 'fontsize', fontSize + 8);

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
        'fontsize'   , fontSize+2              );
    if plotIdx == (nRows*nCols-nCols+1)
        set(gca, 'xticklabel', {'O', 'L_1', 'L_2', 'L_3', 'L_4', 'L_5', 'N'});
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
    y(1) = nanmean(d.decisionON(match(find(d.truthON(match) == 1))));
    for idx = 1:d.nLures
        y(idx+1) = nanmean(d.decisionON(match(find(d.lureON(match) == 1 & d.lureBinON(match) == idx))));
    end
    y(d.nLures+2) = nanmean(d.decisionON(match(find(d.truthON(match) == 2 & d.lureON(match) == 0))));

    H(1) = plot(1:d.nLures+2, 2-y, 'o-',...
        'markersize', 8, ...
        'color', dataColor, ...
        'markeredgecolor', 'w', ...
        'markerfacecolor', dataColor, ...
        'linewidth', 2);

    for idx = 1:d.nLures+2
        bounds = prctile(chains.(sprintf('pp_%d_%d', part, idx))(:), CIbounds);
        mn = pp(part, idx);

        H(2) = plot(1:d.nLures+2, pp(part,:), '+',...
            'markersize', 10, ...
            'color', modelColor, ...
            'linewidth', 2);

        plot(idx+[-whiskerWidth +whiskerWidth], [bounds(1) bounds(1)], '-', ...
            'color', modelColor, ...
            'linewidth', 0.5);
        plot(idx+[-whiskerWidth +whiskerWidth], [bounds(2) bounds(2)], '-', ...
            'color', modelColor, ...
            'linewidth', 0.5);

    end
end

L =  legend(H, 'data', 'model', ...
    'box', 'off', ...
    'location', 'east', ...
    'fontsize', fontSize);
set(L, 'pos', get(L, 'pos') + [0.01 0 0 0]);

[~, AH(1)] = suplabel('Item Type');
set(AH(1), 'pos', get(AH(1), 'pos') + [0 0.025 0]);
[~, AH(2)] = suplabel('Probability Old', 'y');
set(AH(2), 'pos', get(AH(2), 'pos') + [0.025 0 0]);
set(AH, 'fontsize', fontSize + 12);

% print
if printFigures
    warning off;
    print(sprintf('%s/fit_%s_%s.png', figureDir, dataName, modelName), '-dpng');
    print(sprintf('%s/fit_%s_%s.eps', figureDir, dataName, modelName), '-depsc');
    warning on;
end
