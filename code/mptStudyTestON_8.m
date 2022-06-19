%% MPT model of study/test MST with old-new responses
% hierarchical latent mixture model
% logistic model of discriminability

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
modelName = 'mptStudyTestON_8';

% parameters to monitor
params = {...
    'muRho', 'muPsi', 'muDelta', 'muGamma', ......
    'rho', 'psi', 'delta', 'gamma', 'beta', 'tau', ...
    'z', 'phi', ...
    'pp', ...
    };

% MCMC properties
nChains    = 8;     % number of MCMC chains
nBurnin    = 5e3;   % number of discarded burn-in samples
nSamples   = 5e3;   % number of collected samples
nThin      = 5;     % number of samples between those collected
doParallel = 1;     % whether MATLAB parallel toolbox parallizes chains

% recode
truth = d.truthON;
truth(find(d.lureON == 1)) = 3;

% remove contaminants
fileName = sprintf('%s_%s_%s.mat', 'mptStudyTestON_7', dataName, engine);
load(sprintf('%s/%s', storageDir, fileName), 'chains', 'stats', 'diagnostics', 'info');
zMode = codatable(chains, 'z', @mode);
keep = find(ismember(zMode, [1 2]));

% assign MATLAB variables to the observed nodes
data = struct(...
    'y'              , 2 - d.decisionON(ismember(d.participantON, keep))                  , ...
    'truth'          , truth(ismember(d.participantON, keep))                          , ...
    'lureBin'        , d.lureBinON(ismember(d.participantON, keep))                       , ...
    'nParticipants'  , d.nParticipants                , ...
    'p'              , d.participantON(ismember(d.participantON, keep))                   , ...
    'nLures'         , d.nLures             , ...
    'nTotalTrials'   , length(d.trialON(ismember(d.participantON, keep)))                 );

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

zMean = codatable(chains, 'z', @mean);
zMode = codatable(chains, 'z', @mode);
pp = get_matrix_from_coda(chains, 'pp', @mean);
beta = codatable(chains, 'beta', @mean);
tau = codatable(chains, 'tau', @mean);
betaMedian = codatable(chains, 'beta', @median);
tauMedian = codatable(chains, 'tau', @median);
gamma = get_matrix_from_coda(chains, 'gamma', @mean);
rhoTmp = get_matrix_from_coda(chains, 'rho', @mean);
rho = nan(d.nParticipants, 1);
for idx = 1:d.nParticipants
    if ismember(zMode(idx), [1 2])
        rho(idx) = rhoTmp(idx, zMode(idx));
    end
end

nComponents = 5;

% user constants
markerSize = 8;
fontSize = 14;
modelColor = pantone.ClassicBlue;
dataColor = pantone.Tangerine;
violinColor = pantone.DuskBlue; imposeColor = pantone.ClassicBlue;
modelNames = {'low \rho', 'high \rho', 'guess', 'old', 'new'};
modelColors = {pantone.ClassicBlue; pantone.DuskBlue; pantone.Greenery; pantone.Custard; pantone.Tangerine};
binWidth = 0.02; violinScale = 5;
barWidth = [1 0.6 0.3 0.3 0.3];
barOffset = [0 0 0.005 0 -0.005];
CIbounds = [25 75];
whiskerWidth = 0.25;
nShowSamples = 50;
xx = 0:0.01:1;

rng(1);

% derived
binsC = binWidth/2:binWidth:1-binWidth/2;
binsE = 0:binWidth:1;

% which participants in which order
[~, sortIdx] = sort(d.participantCorrectON, 'descend');
[nRows, nCols] = subplotArrange(length(keep));

%% Analysis with empirical measures

% Compare \rho with dpTF and dpTL
tmp = corrcoef(rho(keep), d.dpTF(keep));
fprintf('The correlation of rho with dprime for old vs new is %1.2f\n', tmp(1, 2));
tmp = corrcoef(rho(keep), d.dpTL(keep));
fprintf('The correlation of rho with dprime for old vs lure is %1.2f\n', tmp(1, 2));

%% rho tau breakdown

select = [11 12 5 2 ];
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
    rhoBounds = prctile(chains.(sprintf('rho_%d_%d', part, zMode(part)))(:), CIbounds);
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
        % end
        text(rho(part), tau(part), sprintf(' %s', char(64+find(sortIdx == part))), 'fontsize', fontSize);
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

    text(1, 1.1, sprintf('%s', char(64+find(sortIdx == part))), ...
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
    end

end

% print
if printFigures
    warning off;
    print(sprintf('%s/rhoTau_%s_%s.png', figureDir, dataName, modelName), '-dpng');
    print(sprintf('%s/rhoTau_%s_%s.eps', figureDir, dataName, modelName), '-depsc');
    warning on;
end

%% delta

F = figure; clf; hold on;
setFigure(F, [0.2 0.2 0.6 0.6], '');

for plotIdx = 1:length(keep)

    part = sortIdx(keep(plotIdx));

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

    text(1, 1.1, sprintf('%s', char(64+find(sortIdx(keep) == part))), ...
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

    yy = 1./(1+exp(-betaMedian(part)*(xx - tauMedian(part))));
    plot(xx*d.nLures, yy, '-', ...
        'color', imposeColor, ...
        'linewidth', 4);

end

[~, AH(1)] = suplabel('Lure Level');
set(AH(1), 'pos', get(AH(1), 'pos') + [0 0.025 0]);
[~, AH(2)] = suplabel('Discriminability', 'y');
set(AH(2), 'pos', get(AH(2), 'pos') + [0.025 0 0]);

set(AH, 'fontsize', fontSize + 12);

% print
if printFigures
    warning off;
    print(sprintf('%s/deltas_%s_%s.png', figureDir, dataName, modelName), '-dpng');
    print(sprintf('%s/deltas_%s_%s.eps', figureDir, dataName, modelName), '-depsc');
    warning on;
end


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
    rhoBounds = prctile(chains.(sprintf('rho_%d_%d', part, zMode(part)))(:), CIbounds);
    tauBounds = prctile(chains.(sprintf('tau_%d', part))(:), CIbounds);

    plot([rho(part) rho(part)], tauBounds, '-', ...
        'color', pantone.GlacierGray);
    plot(rhoBounds, [tau(part) tau(part)], '-', ...
        'color', pantone.GlacierGray);
end

for idx = 1:length(sortIdx)
    part = sortIdx(idx);
    H = plot(rho(part), tau(part), 'o', ...
        'markerfacecolor', pantone.ClassicBlue, ...
        'markeredgecolor', 'w');
    text(rho(part), tau(part), sprintf(' %s', char(64+part)), 'fontsize', fontSize);
end

AX2 = subplot(122); hold on;
set(gca, ...
    'xlim'       , [0 20]               , ...
    'xtick'      , 0:5:20          , ...
    'ylim'       , [0 1]         , ...
    'ytick'      , 0:0.2:1 , ...
    'box'        , 'off'                     , ...
    'tickdir'    , 'out'                     , ...
    'layer'      , 'top'                     , ...
    'ticklength' , [0.01 0]                  , ...
    'layer'      , 'top'                     , ...
    'fontsize'   , fontSize                  );
axis square;
xlabel('Slope (\beta)', 'fontsize', fontSize+8);
ylabel('Threshold (\tau)', 'fontsize', fontSize+8);
set(gca, 'pos', get(gca, 'pos') + [0.025 0.025 0 0]);
Raxes(gca, 0.015, 0.01);

% draw posterior densities
for idx = 1:length(sortIdx)
    part = sortIdx(idx);
    rhoBounds = prctile(chains.(sprintf('beta_%d', part))(:), CIbounds);
    tauBounds = prctile(chains.(sprintf('tau_%d', part))(:), CIbounds);

    plot([beta(part) beta(part)], tauBounds, '-', ...
        'color', pantone.GlacierGray);
    plot(rhoBounds, [tau(part) tau(part)], '-', ...
        'color', pantone.GlacierGray);
end

for idx = 1:length(sortIdx)
    part = sortIdx(idx);
    H = plot(beta(part), tau(part), 'o', ...
        'markerfacecolor', pantone.ClassicBlue, ...
        'markeredgecolor', 'w');
    text(beta(part), tau(part), sprintf(' %s', char(64+find(sortIdx(keep) == part))), 'fontsize', fontSize);
end

% print
if printFigures
    warning off;
    print(sprintf('%s/parameters_%s_%s.png', figureDir, dataName, modelName), '-dpng');
    print(sprintf('%s/parameters_%s_%s.eps', figureDir, dataName, modelName), '-depsc');
    warning on;
end


%% fit
CIbounds = [2.5 97.5];
whiskerWidth = 0.25;
clear H;

F = figure; clf; hold on;
setFigure(F, [0.1 0.2 0.7 0.6], '');

for plotIdx = 1:length(keep)

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

    text(1, 1.1, sprintf('%s', char(64+find(sortIdx == part))), ...
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
set(L, 'pos', get(L, 'pos') + [0.08 0 0 0]);

[~, AH(1)] = suplabel('Stimulus Type');
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

