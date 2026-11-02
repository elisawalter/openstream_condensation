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
% - NUSSELT_CODE  — CONDHTM = NUSSELT as implemented (Kutateladze/Labuntsov falling film)
% - SHAH_CODE     — CONDHTM = SHAH as implemented (Dittus-Boelter x (1+1.8/Xtt^0.8))
% - KM_CODE       — CONDHTM = KIMMUDAWAR as implemented (0.048 Re_f^0.69 Pr^0.34)
% - SHAH1979      — Shah (1979)
% - SHAH2009V     — Shah (2009), vertical tube regimes I/II
% - KM2013        — Kim & Mudawar (2013) annular, with phi_g/X_tt
% - DOBSONCHATO   — Dobson & Chato (1998) annular
% - MAX_NU_SHAH   — max(NUSSELT_CODE, SHAH1979), RELAP5-type blending

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
    runId = S{i,2}; D = S{i,6}; P = S{i,10}*1E3; G = S{i,11}; xout = S{i,14};
    rows = cellfun(@(v) isequal(v,runId), L(:,2));
    loc  = cell2mat(L(rows,[5 7 8 13 15 16 17 20]));
    loc  = loc(loc(:,3) == 1 & loc(:,4) > 0 & loc(:,4) < 1,:);
    n = size(loc,1);
    pts = [pts; table(repmat(runId,n,1), repmat(P,n,1), repmat(G,n,1), repmat(D,n,1), ...
        repmat(xout,n,1), loc(:,1), loc(:,2), loc(:,4), loc(:,6), loc(:,7), ...
        'VariableNames',{'run','P','G','D','xout','node','L','x','HF','HTCexp'})]; %#ok<AGROW>
end

% Saturation properties
np = height(pts);
[rl,rv,mul,muv,kl,prl,sig,pr] = deal(zeros(np,1));
Pc = 22.064E6;                                                             % [Pa] Critical pressure
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
Rev = G.*x.*D./muv;                                                        % [-] Vapor Reynolds number
Relo = G.*D./mul;                                                          % [-] Liquid-only Reynolds number
Xtt = ((1-x)./x).^0.9.*(rv./rl).^0.5.*(mul./muv).^0.1;                     % [-] Martinelli parameter
nul = mul./rl;
lNu = (nul.^2/g).^(1/3);                                                   % [m] Nusselt length scale

H = struct();

% NUSSELT as implemented
hstar = zeros(np,1);
iLam = Ref < 30; iWav = Ref >= 30 & Ref < 1600; iTur = Ref >= 1600;
hstar(iLam) = 1.47.*Ref(iLam).^(-1/3);
hstar(iWav) = Ref(iWav)./(1.08.*Ref(iWav).^1.22 - 5.2);
hstar(iTur) = Ref(iTur)./(8750 + 58.*prl(iTur).^(-0.5).*(Ref(iTur).^0.75 - 253));
H.NUSSELT_CODE = hstar.*kl./lNu;

% SHAH as implemented
hDB = 0.023.*Relo.^0.8.*prl.^0.4.*kl./D;
H.SHAH_CODE = hDB.*(1 + 1.8./Xtt.^0.8);

% KIMMUDAWAR as implemented
H.KM_CODE = 0.048.*max(1,Ref).^0.69.*prl.^0.34.*kl./D;

% Shah (1979)
H.SHAH1979 = hDB.*((1-x).^0.8 + 3.8.*x.^0.76.*(1-x).^0.04./pr.^0.38);

% Shah (2009), vertical tubes
nS  = 0.0058 + 0.557.*pr;
hI  = hDB.*(mul./(14.*muv)).^nS.*((1-x).^0.8 + 3.8.*x.^0.76.*(1-x).^0.04./pr.^0.38);
hNu = 1.32.*Ref.^(-1/3).*(rl.*(rl-rv).*g.*kl.^3./mul.^2).^(1/3);
Z   = (1./x-1).^0.8.*pr.^0.4;
Jg  = x.*G./sqrt(g.*D.*rv.*(rl-rv));
H.SHAH2009V = hI + hNu.*(Jg < 1./(2.4.*Z + 0.73));

% Kim & Mudawar (2013) annular, two-phase multiplier from Kim & Mudawar (2012)
Suvo = rv.*sig.*D./muv.^2;
Relo_ = Relo;
ff = 16./max(Ref,1);  ff(Ref >= 2000) = 0.079.*Ref(Ref >= 2000).^-0.25;
fg = 16./max(Rev,1);  fg(Rev >= 2000 & Rev < 2E4) = 0.079.*Rev(Rev >= 2000 & Rev < 2E4).^-0.25;
fg(Rev >= 2E4) = 0.046.*Rev(Rev >= 2E4).^-0.2;
X2 = (ff.*(1-x).^2./rl)./(fg.*x.^2./rv);
X  = sqrt(X2);
C  = zeros(np,1);
ltt = Ref >= 2000 & Rev >= 2000; lvt = Ref < 2000 & Rev >= 2000;
ltv = Ref >= 2000 & Rev < 2000;  lvv = Ref < 2000 & Rev < 2000;
C(ltt) = 0.39.*Relo_(ltt).^0.03.*Suvo(ltt).^0.10.*(rl(ltt)./rv(ltt)).^0.35;
C(lvt) = 8.7E-4.*Relo_(lvt).^0.17.*Suvo(lvt).^0.50.*(rl(lvt)./rv(lvt)).^0.14;
C(ltv) = 0.0015.*Relo_(ltv).^0.59.*Suvo(ltv).^0.19.*(rl(ltv)./rv(ltv)).^0.36;
C(lvv) = 3.5E-5.*Relo_(lvv).^0.44.*Suvo(lvv).^0.50.*(rl(lvv)./rv(lvv)).^0.48;
phig = sqrt(1 + C.*X + X2);
H.KM2013 = 0.048.*Ref.^0.69.*prl.^0.34.*phig./Xtt.*kl./D;

% Dobson & Chato (1998) annular
H.DOBSONCHATO = 0.023.*Ref.^0.8.*prl.^0.4.*(1 + 2.22./Xtt.^0.89).*kl./D;

% RELAP5-type blending
H.MAX_NU_SHAH = max(H.NUSSELT_CODE, H.SHAH1979);

%% Statistics

names = string(fieldnames(H));
ratio = zeros(np,numel(names));
for m = 1:numel(names)
    ratio(:,m) = H.(names(m))./pts.HTCexp;
    pts.(names(m)) = H.(names(m));
end
pts.Ref = Ref; pts.Jg = Jg; pts.Xtt = Xtt;
relErr = ratio - 1;
twoPh = pts.xout > 0;
stats = table(names, ...
    mean(relErr).'*100, sqrt(mean(relErr.^2)).'*100, mean(abs(relErr) <= 0.3).'*100, ...
    mean(relErr(twoPh,:)).'*100, sqrt(mean(relErr(twoPh,:).^2)).'*100, ...
    'VariableNames',{'Model','MeanErr_pct','RMSErr_pct','Within30_pct','MeanErr_xout_gt0_pct','RMSErr_xout_gt0_pct'});
fprintf('%d valid two-phase points from %d runs\n', np, numel(unique(pts.run)));
disp(stats)
writetable(stats,fullfile(outputDir,'correlation_statistics.csv'));
writetable(pts,fullfile(outputDir,'correlation_points.csv'));

%% Figures

fig = figure('Visible','off','Position',[50 50 1600 800]);
tl = tiledlayout(fig,2,4,'TileSpacing','compact');
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
im = find(names == "NUSSELT_CODE");
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
