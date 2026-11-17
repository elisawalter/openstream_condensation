%% Stand-alone assessment of condensation HTC correlations against PICON-V
%
% Evaluates condensation heat transfer correlations at the measured local
% conditions (pressure, mass flux, local quality) of the valid PICON-V
% pure-steam measurement sections and compares them with the measured local
% HTC. This isolates the closure from the solver: with TWALLIN boundary
% conditions the predicted local heat flux error equals the HTC error at the
% predicted local quality.
%
% Correlations:
%
% - NUSSELT, SHAH, SHAH2009, KIMMUDAWAR, NUSSELTSHAH — the CONDHTM models,
%   evaluated with :meth:`Solvers.Mixture.Mixture.CONDHTC` (same code as the solver)
% - DOBSONCHATO   — Dobson & Chato (1998) annular (reference only)
% - MAX_NU_KM     — max(NUSSELT, KIMMUDAWAR) blending (reference only)
%
% Statistics are given for all valid sections and for runs with all
% sections flagged valid (no complete condensation).

scriptDir = fileparts(mfilename('fullpath'));
addpath(fileparts(scriptDir));
dataFile  = fullfile(scriptDir,'PICONV_data_Ahn.xlsx');
outputDir = fullfile(scriptDir,'outputs','validation_piconv_mixture');
if ~isfolder(outputDir), mkdir(outputDir); end

% Pure-steam runs excluding the workbook-flagged groups (42-112, 161-196)
S = readcell(dataFile,'Sheet','PICONV_summary (118)');
S = S(4:end,:);
isPst = cellfun(@(v) isnumeric(v) && isequal(v,0), S(:,3));
S = S(isPst,:);
runs = cell2mat(S(:,2));
keep = ~(runs >= 42 & runs <= 112) & ~(runs >= 161 & runs <= 196);
S = S(keep,:);

L = readcell(dataFile,'Sheet','Local data(1180)');
L = L(4:end,:);

cp = CoolPropWrapper.CoolPropWrapper('WATER');
g  = 9.81;                                                                 % [m/s^2]

pts = table();
for i = 1:size(S,1)
    runId = S{i,2}; D = S{i,6}; P = S{i,10}*1E3; G = S{i,11};
    rows = cellfun(@(v) isequal(v,runId), L(:,2));
    loc  = cell2mat(L(rows,[5 7 8 13 15 16 17 20]));
    allValid = all(loc(:,3) == 1);
    loc  = loc(loc(:,3) == 1 & loc(:,4) > 0 & loc(:,4) < 1,:);
    n = size(loc,1);
    pts = [pts; table(repmat(runId,n,1), repmat(P,n,1), repmat(G,n,1), repmat(D,n,1), ...
        repmat(allValid,n,1), loc(:,1), loc(:,2), loc(:,4), loc(:,6), loc(:,7), ...
        'VariableNames',{'run','P','G','D','allValid','node','L','x','HF','HTCexp'})]; %#ok<AGROW>
end

% Saturation properties
np = height(pts);
[rl,rv,mul,muv,kl,prl,sig,pr] = deal(zeros(np,1));
Pc = cp.CoolProp.p_critical;                                              % [Pa] Critical pressure
for k = 1:np
    rl(k)  = cp.density('P',pts.P(k),'Q',0);
    rv(k)  = cp.density('P',pts.P(k),'Q',1);
    mul(k) = cp.viscosity('P',pts.P(k),'Q',0);
    muv(k) = cp.viscosity('P',pts.P(k),'Q',1);
    kl(k)  = cp.conductivity('P',pts.P(k),'Q',0);
    prl(k) = cp.prandtl('P',pts.P(k),'Q',0);
    sig(k) = cp.surfaceTension('P',pts.P(k),'Q',0);
    pr(k)  = pts.P(k)/Pc;
end

x = pts.x; G = pts.G; D = pts.D;
Ref = G.*(1-x).*D./mul;                                                    % [-] Film Reynolds number, 4*Gamma/mu
Xtt = ((1-x)./x).^0.9.*(rv./rl).^0.5.*(mul./muv).^0.1;                     % [-] Martinelli parameter
Jg  = x.*G./sqrt(g.*D.*rv.*(rl-rv));                                       % [-] Dimensionless vapor velocity
props = struct('RHOL',rl,'RHOV',rv,'MUL',mul,'MUV',muv,'KL',kl,'PRL',prl,'SIGMA',sig,'PR',pr);

H = struct();

% CONDHTM models, same code as the solver
for m = ["NUSSELT","SHAH","SHAH2009","KIMMUDAWAR","NUSSELTSHAH"]
    H.(m) = Solvers.Mixture.Mixture.CONDHTC(m, x, Ref, G, D, g, props);
end

% Dobson & Chato (1998) annular, reference only
H.DOBSONCHATO = 0.023.*Ref.^0.8.*prl.^0.4.*(1 + 2.22./Xtt.^0.89).*kl./D;

% Nusselt and Kim-Mudawar blending, reference only
H.MAX_NU_KM = max(H.NUSSELT, H.KIMMUDAWAR);

%% Statistics

names = string(fieldnames(H));
ratio = zeros(np,numel(names));
for m = 1:numel(names)
    ratio(:,m) = H.(names(m))./pts.HTCexp;
    pts.(names(m)) = H.(names(m));
end
pts.Ref = Ref; pts.Jg = Jg; pts.Xtt = Xtt;
relErr = ratio - 1;
twoPh = pts.allValid;
stats = table(names, ...
    mean(relErr).'*100, sqrt(mean(relErr.^2)).'*100, mean(abs(relErr) <= 0.3).'*100, ...
    mean(relErr(twoPh,:)).'*100, sqrt(mean(relErr(twoPh,:).^2)).'*100, ...
    'VariableNames',{'Model','MeanErr_pct','RMSErr_pct','Within30_pct','MeanErr_allValid_pct','RMSErr_allValid_pct'});
fprintf('%d valid two-phase points from %d runs\n', np, numel(unique(pts.run)));
disp(stats)
writetable(stats,fullfile(outputDir,'correlation_statistics.csv'));
writetable(pts,fullfile(outputDir,'correlation_points.csv'));

%% Figures

fig = figure('Visible','off','Position',[50 50 1600 800]);
tl = tiledlayout(fig,2,ceil(numel(names)/2),'TileSpacing','compact');
lim = [0 1.1*max([pts.HTCexp; ratio(:).*repmat(pts.HTCexp,numel(names),1)])]/1E3;
lim(2) = min(lim(2), 3*max(pts.HTCexp)/1E3);
for m = 1:numel(names)
    ax = nexttile(tl); hold(ax,'on'); box(ax,'on'); axis(ax,'square');
    scatter(ax, pts.HTCexp/1E3, H.(names(m))/1E3, 14, pts.G, 'filled');
    plot(ax, lim, lim, 'k-'); plot(ax, lim, 1.3*lim, 'k--'); plot(ax, lim, 0.7*lim, 'k--');
    xlim(ax,lim); ylim(ax,lim);
    xlabel(ax,'Measured h [kW/m^2/K]'); ylabel(ax,'Predicted h [kW/m^2/K]');
    title(ax, sprintf('%s: %+.0f%% / RMS %.0f%%', strrep(names(m),'_','\_'), stats.MeanErr_pct(m), stats.RMSErr_pct(m)));
end
cb = colorbar(ax); cb.Label.String = 'G [kg/m^2/s]';
exportgraphics(fig, fullfile(outputDir,'correlation_parity.png'), 'Resolution', 130);
close(fig);

% Error trends for the Nusselt film model
fig = figure('Visible','off','Position',[50 50 1400 400]);
tl = tiledlayout(fig,1,3,'TileSpacing','compact');
im = find(names == "NUSSELT");
trend = {pts.x,'x [-]'; pts.Jg,'J_g^* [-]'; pts.Ref,'Re_f [-]'};
for k = 1:3
    ax = nexttile(tl); hold(ax,'on'); box(ax,'on');
    scatter(ax, trend{k,1}, ratio(:,im), 14, pts.P/1E5, 'filled');
    yline(ax,1,'k-'); yline(ax,[0.7 1.3],'k--');
    xlabel(ax,trend{k,2}); ylabel(ax,'h_{Nusselt}/h_{exp} [-]');
end
cb = colorbar(ax); cb.Label.String = 'P [bar]';
exportgraphics(fig, fullfile(outputDir,'nusselt_error_trends.png'), 'Resolution', 130);
close(fig);
