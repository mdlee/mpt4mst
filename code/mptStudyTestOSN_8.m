%% MPT model of study/test with OSN
% v8 add sigmoid delta and assume contaminants removed

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
modelName = 'mptStudyTestOSN_8';

% parameters to monitor
params = {...
    'muRho', 'muPsi', 'muDelta', 'muGammaO', 'muGammaS', ...
    'rho', 'psi', 'delta', 'gammaO', 'gammaS', 'beta', 'tau', ...
    'ppO', 'ppS', ...
    };

% MCMC properties
nChains    = 8;     % number of MCMC chains
nBurnin    = 1e3;   % number of discarded burn-in samples
nSamples   = 1e3;   % number of collected samples
nThin      = 1;     % number of samples between those collected
doParallel = 1;     % whether MATLAB parallel toolbox parallizes chains

% remove contaminants
fileName = sprintf('%s_%s_%s.mat', 'mptStudyTestOSN_7', dataName, engine);
load(sprintf('%s/%s', storageDir, fileName), 'chains', 'stats', 'diagnostics', 'info');
zMode = codatable(chains, 'z', @mode);
keep = find(ismember(zMode, [1 2]));

% assign MATLAB variables to the observed nodes
data = struct(...
    'y'              , d.decisionOSN(ismember(d.participantOSN, keep))                  , ...
    'truth'          , d.truthOSN(ismember(d.participantOSN, keep))                          , ...
    'lureBin'        , d.lureBinOSN(ismember(d.participantOSN, keep))                       , ...
    'nParticipants'  , d.nParticipants                , ...
    'p'              , d.participantOSN(ismember(d.participantOSN, keep))                   , ...
    'nLures'         , d.nLures             , ...
    'nTotalTrials'   , length(d.trialOSN(ismember(d.participantOSN, keep)))                 );


% generator for initialization
generator = @()struct('muPsi', rand);

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

    yPred = get_matrix_from_coda(chains, 'yPred');

    prefix = 'yPred';
    fn = fieldnames(chains);
    tf = strncmp(fn, prefix, length(prefix));
    chains = rmfield(chains, fn(tf));

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
ppO = get_matrix_from_coda(chains, 'ppO', @mean);
ppS = get_matrix_from_coda(chains, 'ppS', @mean);
beta = codatable(chains, 'beta', @mean);
tau = codatable(chains, 'tau', @mean);
betaMedian = codatable(chains, 'beta', @median);
tauMedian = codatable(chains, 'tau', @median);
rho = get_matrix_from_coda(chains, 'rho');

[~, sortIdxFull] = sort(d.participantCorrectON, 'descend');

% Compare OSN inferences with ON
modelName = 'mptStudyTestON_8';
fileName = sprintf('%s_%s_%s.mat', modelName, dataName, engine);
load(sprintf('%s/%s', storageDir, fileName), 'chains', 'stats', 'diagnostics', 'info');

betaON = codatable(chains, 'beta', @mean);
tauON = codatable(chains, 'tau', @mean);
rhoONtmp = get_matrix_from_coda(chains, 'rho', @mean);
z = get_matrix_from_coda(chains, 'z', @mode);
for idx = 1:length(z)
    rhoON(idx) = rhoONtmp(find(sortIdxFull == idx), z(find(sortIdxFull == idx)));
end

modelName = 'mptStudyTestOSN_8';
fileName = sprintf('%s_%s_%s.mat', modelName, dataName, engine);
load(sprintf('%s/%s', storageDir, fileName), 'chains', 'stats', 'diagnostics', 'info');

tmp = corrcoef(rho(keep), rhoON(keep));
fprintf('The correlation of rho between ON and OSN tasks is %1.2f\n', tmp(1, 2));
tmp = corrcoef(tau(keep), tauON(keep));
fprintf('The correlation of tau between ON and OSN tasks is %1.2f\n', tmp(1, 2));

% Compare \rho with REC and \tau with LDI
tmp = corrcoef(rho(keep), d.REC(keep));
fprintf('The correlation of rho with REC is %1.2f\n', tmp(1, 2));
tmp = corrcoef(tau(keep), d.LDI(keep));
fprintf('The correlation of tau with LDI is %1.2f\n', tmp(1, 2));

return

%% Figures

% user constants
markerSize = 8;
fontSize = 14;
modelColor = pantone.ClassicBlue;
dataColor = pantone.Tangerine;
dataSimColor = pantone.DuskBlue;
dataOldColor = pantone.Tangerine;
violinColor = pantone.DuskBlue; imposeColor = pantone.ClassicBlue;
binWidth = 0.02; violinScale = 5;
barWidth = [1 0.6 0.3 0.3 0.3];
CIbounds = [25 75];
whiskerWidth = 0.25;
nShowSamples = 50;
xx = 0:0.01:1;

% derived
binsC = binWidth/2:binWidth:1-binWidth/2;
binsE = 0:binWidth:1;

% which participants in which order
[~, sortIdxFull] = sort(d.participantCorrectON, 'descend');
sortIdx = intersect(sortIdxFull, keep, 'stable');
[nRows, nCols] = subplotArrange(length(sortIdx));
%nRows = 4; nCols = 5;

%% rho tau breakdown

select = [4 3 2 21];
selectSubplot = [3 4 7 8];

% joint posterior figure
F = figure; clf; hold on;
setFigure(F, [0.1 0.2 0.6 0.5], '');
set(gcf, 'renderer', 'opengl');

AX1 = subplot(2, 4, [1 2 5 6]); hold on;
set(gca, ...
    'xlim'       , [0 1]               , ...
    'xtick'      , 0:0.2:1          , ...
    'ylim'       , [0 1]         , ...
    'ytick'      , 0:0.2:1 , ...
    'box'        , 'off'                     , ...
    'tickdir'    , 'out'                     , ...
    'layer'      , 'top'                     , ...
    'ticklength' , [0.01 0]                  , ...
    'layer'      , 'top'                     , ...
    'fontsize'   , fontSize                  );
axis square;
xlabel('Memory (\rho)', 'fontsize', fontSize+8);
ylabel('Threshold (\tau)', 'fontsize', fontSize+8);
set(gca, 'pos', get(gca, 'pos') + [0 0.025 0 0]);
Raxes(gca, 0.015, 0.01);

% draw posterior densities
for idx = 1:length(sortIdx)
    part = sortIdx(idx);
    rhoBounds = prctile(chains.(sprintf('rho_%d', part))(:), CIbounds);
    tauBounds = prctile(chains.(sprintf('tau_%d', part))(:), CIbounds);

    plot([rho(part) rho(part)], tauBounds, '-', ...
        'color', pantone.GlacierGray);
    plot(rhoBounds, [tau(part) tau(part)], '-', ...
        'color', pantone.GlacierGray);
end

for idx = 1:length(sortIdx)
    part = sortIdx(idx);
    H = plot(rho(part), tau(part), 'o', ...
        'markeredgecolor', pantone.ClassicBlue, ...
        'markerfacecolor', 'w');
    if ismember(part, select)
        set(H, 'markeredgecolor', 'w', ...
            'markerfacecolor', pantone.ClassicBlue, ...
            'markersize', 8);
        text(rho(part), tau(part), sprintf(' %s',  char(64+find(sortIdxFull == part))), 'fontsize', fontSize);
    end
end

for selectIdx = 1:length(select)
    part = select(selectIdx);
    subplot(2, 4, selectSubplot(selectIdx)); hold on;
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
        'ticklength' , [0.03 0]              , ...
        'layer'      , 'top'                 , ...
        'fontsize'   , fontSize+2              );
    if selectIdx == 3
        set(gca, 'xticklabel', {'O', 'L_1', 'L_2', 'L_3', 'L_4', 'L_5', 'N'});
    else
        set(gca, 'yticklabel', []);
    end
    set(gca, 'pos', get(gca, 'pos').*[1 1 1 0.875] + [0.025 0.04 0 0])
    Raxes(gca, 0.01, 0.01);

    text(1, 1.1, sprintf('%s', char(64+find(sortIdxFull == part))), ...
        'fontsize', fontSize, ...
        'hor', 'cen');

    % old responses for each stimulus type
    match = find(d.participantOSN == part);
    yO(1) = nansum(d.decisionOSN(match) == 1 & d.truthOSN(match) == 1)/nansum(d.truthOSN(match) == 1);
    for idx = 1:d.nLures
        yO(idx+1) = nansum(d.decisionOSN(match) == 1 & d.truthOSN(match) == 3 & d.lureBinOSN(match) == idx)/nansum(d.truthOSN(match) == 3 & d.lureBinOSN(match) == idx);
    end
    yO(d.nLures+2) = nansum(d.decisionOSN(match) == 1 & d.truthOSN(match) == 2)/nansum(d.truthOSN(match) == 2);

    H(1) = plot(1:d.nLures+2, yO, 's-',...
        'markersize', 10, ...
        'color', dataOldColor, ...
        'markeredgecolor', 'w', ...
        'markerfacecolor', dataOldColor, ...
        'linewidth', 2);

    % similar responses for each stimulus type
    match = find(d.participantOSN == part);
    yS(1) = nansum(d.decisionOSN(match) == 3 & d.truthOSN(match) == 1)/nansum(d.truthOSN(match) == 1);
    for idx = 1:d.nLures
        yS(idx+1) = nansum(d.decisionOSN(match) == 3 & d.truthOSN(match) == 3 & d.lureBinOSN(match) == idx)/nansum(d.truthOSN(match) == 3 & d.lureBinOSN(match) == idx);
    end
    yS(d.nLures+2) = nansum(d.decisionOSN(match) == 3 & d.truthOSN(match) == 2)/nansum(d.truthOSN(match) == 2);

    H(3) = plot(1:d.nLures+2, yS, 'o--',...
        'markersize', 6, ...
        'color', dataSimColor, ...
        'markerfacecolor', 'w', ...
        'markeredgecolor', dataSimColor, ...
        'linewidth', 2);

end



if printFigures
    warning off;
    print(sprintf('%s/rhoTau_%s_%s.png', figureDir, dataName, modelName), '-dpng');
    print(sprintf('%s/rhoTau_%s_%s.eps', figureDir, dataName, modelName), '-depsc');
    warning on;
end

%% fitOSN
F = figure; clf; hold on;
setFigure(F, [0.1 0.2 0.7 0.7], '');

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

    text(1, 1.1, sprintf('%s', char(64+find(sortIdxFull == part))), ...
        'fontsize', fontSize, ...
        'hor', 'cen');

    % old responses for each stimulus type
    match = find(d.participantOSN == part);
    yO(1) = nansum(d.decisionOSN(match) == 1 & d.truthOSN(match) == 1)/nansum(d.truthOSN(match) == 1);
    for idx = 1:d.nLures
        yO(idx+1) = nansum(d.decisionOSN(match) == 1 & d.truthOSN(match) == 3 & d.lureBinOSN(match) == idx)/nansum(d.truthOSN(match) == 3 & d.lureBinOSN(match) == idx);
    end
    yO(d.nLures+2) = nansum(d.decisionOSN(match) == 1 & d.truthOSN(match) == 2)/nansum(d.truthOSN(match) == 2);

    H(1) = plot(1:d.nLures+2, yO, 's-',...
        'markersize', 10, ...
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

    H(3) = plot(1:d.nLures+2, yS, 'o--',...
        'markersize', 6, ...
        'color', dataSimColor, ...
        'markerfacecolor', 'w', ...
        'markeredgecolor', dataSimColor, ...
        'linewidth', 2);

    for idx = 1:d.nLures+2
        H(4) = plot(1:d.nLures+2, ppS(part,:), 'x',...
            'markersize', 10, ...
            'color', modelColor, ...
            'linewidth', 2);
    end

end

L =  legend(H, 'old data', 'old model', 'sim data', 'sim model', ...
    'box', 'off', ...
    'location', 'northeast', ...
    'fontsize', fontSize);
set(L, 'pos', get(L, 'pos') + [0.125 0 0 0]);

[~, AH(1)] = suplabel('Item Type');
[~, AH(2)] = suplabel('Probability Old/Similar', 'y');
set(AH(1), 'pos', get(AH(1), 'pos')  + [0 0.02 0]);
set(AH(2), 'pos', get(AH(2), 'pos')  + [0.02 0 0]);
set(AH, 'fontsize', fontSize + 12);

% print
if printFigures
    warning off;
    print(sprintf('%s/fitOSN_%s_%s.png', figureDir, dataName, modelName), '-dpng');
    print(sprintf('%s/fitOSN_%s_%s.eps', figureDir, dataName, modelName), '-depsc');
    warning on;
end
