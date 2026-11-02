%% PICON-V pure-steam downward condensation validation of the mixture solver
%
% Compares the mixture solver condensation HTC models (CONDHTM) against the
% pure-steam tests of the PICON-V facility (Ahn et al., PICONV_data_Ahn.xlsx):
% 40 mm vertical tube, 3 m cooled length, downward steam flow, 10 measurement
% sections.
%
% Each run is simulated with the measured inner wall temperatures as TWALLIN
% boundary condition (:func:`piconv_mixture_case`), so the predicted local
% wall heat flux directly tests the condensation heat transfer coefficient.
%
% Selected runs: all pure-steam tests (ST_mode 0) excluding the run groups
% flagged in the workbook (42-112: coolant leakage and pressure offset;
% 161-196: suspected steam leakage). Sections with complete condensation
% are flagged invalid in the workbook (air is vented into the outlet region
% and forms a stagnant air pillow, which the pure-steam model does not
% represent). Runs with complete condensation are kept for their valid
% upstream sections; statistics are also given for runs with all sections
% valid, and the total heat is only compared for these runs.
% Missing case results are computed on the fly; they can be precomputed in
% parallel MATLAB sessions with piconv_mixture_case(runId).
%
% Only measurement sections flagged valid in the workbook are compared.

% Set simulateMissing = false before running to only post-process existing results.

scriptDir = fileparts(mfilename('fullpath'));
outputDir = fullfile(scriptDir,'outputs','validation_piconv_mixture');
if ~isfolder(outputDir), mkdir(outputDir); end

runIds = [1 113 119 124 129 5 6 7 8 9 134 140 145 13 150 15 16 17 18 19 ...
    197 203 23 208 213 218 224 229 234 235 236];
models = ["NUSSELT","SHAH","KIMMUDAWAR"];
colors = [0.00 0.45 0.74; 0.85 0.33 0.10; 0.47 0.67 0.19];

%% Simulations

if ~exist('simulateMissing','var'), simulateMissing = true; end

results = {};
for i = 1:numel(runIds)
    resFile = fullfile(outputDir,sprintf('PICONV_%03d.mat',runIds(i)));
    if ~isfile(resFile)
        if ~simulateMissing
            fprintf('Skipping PICON-V run %d (no result file)\n',runIds(i));
            continue
        end
        fprintf('Simulating PICON-V run %d ...\n',runIds(i));
        piconv_mixture_case(runIds(i),'outputDir',outputDir,'CONDHTM',models);
    end
    S = load(resFile);
    results{end+1} = S.result; %#ok<SAGROW>
end

%% Local comparison

DSEC = 0.3;                                                                % [m] Measurement section length
rows = struct('run',{},'P',{},'G',{},'allValid',{},'node',{},'L',{},'xexp',{},'HFexp',{},'HFpred',{},'xpred',{},'ok',{});
summary = table();

for i = 1:numel(results)
    r = results{i}; d = r.data; z = r.z;
    HFpred = nan(numel(d.L),numel(models));
    xpred  = nan(numel(d.L),numel(models));
    Qpred  = nan(1,numel(models));
    for m = 1:numel(models)
        s = r.models.(models(m));
        if isempty(s.q), continue; end                                     % Non-converged runs are kept and flagged by their state
        for k = 1:numel(d.L)
            inSec = z > d.L(k)-DSEC/2 & z <= d.L(k)+DSEC/2;
            HFpred(k,m) = -mean(s.q(inSec));                               % [W/m^2] Section-averaged removed heat flux
            xpred(k,m)  = interp1(z,s.xeq,d.L(k));                         % [-] Equilibrium quality at section center
        end
        Qpred(m) = sum(HFpred(:,m))*pi*d.D*DSEC;                           % [W]
    end
    for k = find(d.valid).'
        rows(end+1) = struct('run',r.runId,'P',d.P,'G',d.G,'allValid',all(d.valid),'node',d.node(k),'L',d.L(k), ...
            'xexp',d.x(k),'HFexp',d.HF(k),'HFpred',HFpred(k,:),'xpred',xpred(k,:), ...
            'ok',arrayfun(@(m) r.models.(m).ok, models)); %#ok<SAGROW>
    end
    states = arrayfun(@(m) r.models.(m).state, models);
    summary = [summary; table(r.runId, d.P/1E5, d.G, d.xin, d.xout, all(d.valid), d.Qw/1E3, ...
        Qpred(1)/1E3, Qpred(2)/1E3, Qpred(3)/1E3, states, ...
        'VariableNames',{'run','P_bar','G','x_in','x_out','all_valid','Qexp_kW', ...
        'Q_NUSSELT_kW','Q_SHAH_kW','Q_KIMMUDAWAR_kW','states'})]; %#ok<AGROW>
end

HFexp  = [rows.HFexp].';
HFpred = vertcat(rows.HFpred);
HFpred(~vertcat(rows.ok)) = NaN;                                           % Statistics on converged runs only
relErr = HFpred./HFexp - 1;
allValid = [rows.allValid].';
QErr   = summary{:,{'Q_NUSSELT_kW','Q_SHAH_kW','Q_KIMMUDAWAR_kW'}}./summary.Qexp_kW - 1;
QErr(~summary.all_valid,:) = NaN;                                          % Total heat only for runs with all sections valid
QErr(summary.states ~= "INITIALSTEPCONVERGED") = NaN;
converged = sum(summary.states == "INITIALSTEPCONVERGED").';

stats = table(models.', ...
    converged, ...
    sum(~isnan(relErr)).', ...
    mean(relErr,'omitnan').'*100, ...
    sqrt(mean(relErr.^2,'omitnan')).'*100, ...
    (sum(abs(relErr) <= 0.3)./sum(~isnan(relErr))).'*100, ...
    mean(relErr(allValid,:),'omitnan').'*100, ...
    sqrt(mean(relErr(allValid,:).^2,'omitnan')).'*100, ...
    sum(~isnan(QErr)).', ...
    mean(QErr,'omitnan').'*100, ...
    sqrt(mean(QErr.^2,'omitnan')).'*100, ...
    'VariableNames',{'Model','ConvergedRuns','Points','MeanErr_pct','RMSErr_pct','Within30_pct', ...
    'MeanErr_allValid_pct','RMSErr_allValid_pct','QRuns','MeanErrQ_pct','RMSErrQ_pct'});

disp(summary)
disp(stats)
writetable(summary,fullfile(outputDir,'summary.csv'));
writetable(stats,fullfile(outputDir,'statistics.csv'));

%% Figures

% Parity plot of local wall heat flux
fig = figure('Visible','off','Position',[100 100 1350 450]);
tl = tiledlayout(fig,1,numel(models),'TileSpacing','compact');
lim = [0 1.1*max(HFexp)]/1E3;
for m = 1:numel(models)
    ax = nexttile(tl); hold(ax,'on'); box(ax,'on'); axis(ax,'square');
    scatter(ax, HFexp/1E3, HFpred(:,m)/1E3, 18, [rows.G].', 'filled');
    plot(ax, lim, lim, 'k-'); plot(ax, lim, 1.3*lim, 'k--'); plot(ax, lim, 0.7*lim, 'k--');
    xlim(ax,lim); ylim(ax,lim);
    xlabel(ax,'Measured q_{wall} [kW/m^2]'); ylabel(ax,'Predicted q_{wall} [kW/m^2]');
    title(ax, sprintf('%s: mean %+.0f%%, RMS %.0f%%', models(m), stats.MeanErr_pct(m), stats.RMSErr_pct(m)));
end
cb = colorbar(ax); cb.Label.String = 'G [kg/m^2/s]';
exportgraphics(fig, fullfile(outputDir,'piconv_parity.png'), 'Resolution', 150);
close(fig);

% Axial profiles of wall heat flux and quality for a representative subset
showRuns = [119 9 150 17 19 229 15 197 218 234];
fig = figure('Visible','off','Position',[50 50 1600 1100]);
tl = tiledlayout(fig,4,5,'TileSpacing','compact');
iShow = find(ismember(cellfun(@(r) r.runId, results), showRuns));
for j = 1:numel(iShow)
    r = results{iShow(j)}; d = r.data;
    ax = nexttile(tl,j + 5*floor((j-1)/5)); hold(ax,'on'); box(ax,'on');
    for m = 1:numel(models)
        s = r.models.(models(m));
        if ~isempty(s.q), plot(ax, r.z, -s.q/1E3, 'Color', colors(m,:), 'LineWidth', 1.2, 'DisplayName', models(m)); end
    end
    plot(ax, d.L(d.valid), d.HF(d.valid)/1E3, 'ko', 'MarkerFaceColor', 'k', 'DisplayName', 'PICON-V');
    plot(ax, d.L(~d.valid), d.HF(~d.valid)/1E3, 'ko', 'DisplayName', 'PICON-V (excluded)');
    title(ax, sprintf('Run %d: %.1f bar, G = %.0f', r.runId, d.P/1E5, d.G));
    ylabel(ax,'q_{wall} [kW/m^2]'); xlim(ax,[0 3]); ylim(ax,[0 1.3*max(d.HF)/1E3]);
    if j == 1, legend(ax,'Location','best','FontSize',7); end

    ax = nexttile(tl,j + 5*ceil(j/5)); hold(ax,'on'); box(ax,'on');
    for m = 1:numel(models)
        s = r.models.(models(m));
        if ~isempty(s.x), plot(ax, r.z, s.xeq, 'Color', colors(m,:), 'LineWidth', 1.2); end
    end
    plot(ax, d.L, d.x, 'ko', 'MarkerFaceColor', 'k');
    ylabel(ax,'x [-]'); xlabel(ax,'z [m]'); xlim(ax,[0 3]);
end
exportgraphics(fig, fullfile(outputDir,'piconv_profiles.png'), 'Resolution', 120);
close(fig);

fprintf('PICON-V mixture validation results saved to %s\n', outputDir);
