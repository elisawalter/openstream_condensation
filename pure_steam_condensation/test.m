%% Downward pure-steam condensation integration test
%
% This script validates the downward-condensation paths recently added to
% the mixture solver:
% - Negative POWER/HFLUX handling
% - Condensation HTC models (NONE, NUSSELT, SHAH, KIMMUDAWAR)
% - Void model extension (ZIVI)
% - Bestion drift-velocity angle dependence in downflow

clearvars

scriptPath = mfilename('fullpath');
scriptDir = fileparts(scriptPath);
openstreamRoot = fileparts(scriptDir);
addpath(openstreamRoot);

import Inputs.*
import Solvers.*
import Solvers.Mixture.*

inputsDir = fullfile(openstreamRoot,'inputs');
outputDir = fullfile(scriptDir,'outputs','integration_downward_condensation');

if isfolder(outputDir)
    rmdir(outputDir,'s');
end
mkdir(outputDir);

tempModelFile = fullfile(outputDir,'models_downward_condensation.inp');
fclose(fopen(tempModelFile,'w'));

% Build model cases that exercise the newly implemented options.
modelCases = [ ...
    struct('id','PURECOND_NUSSELT_ZIVI'   ,'condhtm','NUSSELT'   ,'void','ZIVI'), ...
    struct('id','PURECOND_SHAH_ZIVI'      ,'condhtm','SHAH'      ,'void','ZIVI'), ...
    struct('id','PURECOND_KM_ZIVI'        ,'condhtm','KIMMUDAWAR','void','ZIVI'), ...
    struct('id','PURECOND_NONE_ZIVI'      ,'condhtm','NONE'      ,'void','ZIVI'), ...
    struct('id','PURECOND_NUSSELT_BESTION','condhtm','NUSSELT'   ,'void','BESTION') ...
    ];

for idx = 1:numel(modelCases)
    c = modelCases(idx);
    Model.writeInputFile( ...
        tempModelFile, c.id, ...
        'NNODES', 100, ...
        'FLUID', 'WATER', ...
        'VOID', c.void, ...
        'CONDHTM', c.condhtm, ...
        'TPFM', 'SLIP', ...
        'MOMENTLIQUID', 'FULL', ...
        'MOMENTGAS', 'FULL', ...
        'MOMENTFILM', 'FULL', ...
        'MOMENTDROP', 'FULL');
end

warningState = warning;
warning('off','all');
warningCleanup = onCleanup(@() warning(warningState));

results = struct('id',{},'state',{},'maxCondResidual',{},'hasTwoPhaseCondZone',{},'vfMean',{});
plotCases = struct('id',{},'mix',{},'fluid',{},'geom',{},'twall',{},'hwall',{});

for idx = 1:numel(modelCases)
    c = modelCases(idx);

    inputSet = InputSet( ...
        modelFilePath          = tempModelFile, ...
        modelID                = c.id, ...
        optionsFilePath        = fullfile(inputsDir,'options.inp'), ...
        optionsID              = 'DEFAULT', ...
        geometryFilePath       = fullfile(inputsDir,'geom.inp'), ...
        geometryID             = 'PURECOND_DOWN', ...
        bcFilePath             = fullfile(inputsDir,'purecond_downflow.inp'), ...
        sessionParentDir       = outputDir, ...
        overwriteSessionFiles  = false, ...
        LOGMODE                = 'BOTH');

    assert(inputSet.bc(1).POWER < 0, ...
        'Boundary condition POWER must be negative for condensation test.');
    assert(all(inputSet.bc(1).WPOWER(:) >= 0), ...
        'Boundary condition WPOWER must stay non-negative.');

    mixSolver = MixtureSolver(inputSet);
    mixSolver.solve();

    assert(mixSolver.STATE == Solvers.SolverState.INITIALSTEPCONVERGED, ...
        'Case %s did not converge to INITIALSTEPCONVERGED.', c.id);

    mix = mixSolver.mixture(1);
    q = mix.HFLUX;
    x = repmat(mix.X(),1,size(q,2));

    assert(any(q(:) < 0), ...
        'Case %s has no negative wall heat flux.', c.id);

    [twall,~] = mix.TWALL();
    assert(all(isfinite(twall(:))), ...
        'Case %s produced non-finite wall temperature.', c.id);
    assert(isreal(twall), ...
        'Case %s produced complex wall temperature.', c.id);

    idxCond = q < 0 & x > 0 & x < 1;
    hasTwoPhaseCondZone = any(idxCond(:));
    assert(hasTwoPhaseCondZone, ...
        'Case %s did not produce a two-phase condensation region (0 < x < 1 with q < 0).', c.id);

    [hwall,~] = mix.HWALL(twall);
    condResidual = abs(twall(idxCond) - (mixSolver.fluid.TSAT + q(idxCond)./hwall(idxCond)));
    maxCondResidual = max(condResidual,[],'omitnan');

    assert(maxCondResidual < 1E-2, ...
        'Case %s violates TWALL condensation relation. Max residual: %.3e K', c.id, maxCondResidual);

    results(end+1) = struct( ...
        'id', c.id, ...
        'state', string(mixSolver.STATE), ...
        'maxCondResidual', maxCondResidual, ...
        'hasTwoPhaseCondZone', hasTwoPhaseCondZone, ...
        'vfMean', mean(mix.VF(),'omitnan')); %#ok<SAGROW>

    % Store each case for post-processing and plotting.
    newCaseIdx = numel(plotCases) + 1;
    plotCases(newCaseIdx).id = c.id;
    plotCases(newCaseIdx).mix = mix;
    plotCases(newCaseIdx).fluid = mixSolver.fluid;
    plotCases(newCaseIdx).geom = inputSet.geometry;
    plotCases(newCaseIdx).twall = twall;
    plotCases(newCaseIdx).hwall = hwall;
end

% Explicitly exercise the Bestion angle effect (cos(ANGLE)) by comparing
% the same model in downflow vs upflow geometry.
bestionCaseId = 'PURECOND_NUSSELT_BESTION';

downSet = InputSet( ...
    modelFilePath          = tempModelFile, ...
    modelID                = bestionCaseId, ...
    optionsFilePath        = fullfile(inputsDir,'options.inp'), ...
    optionsID              = 'DEFAULT', ...
    geometryFilePath       = fullfile(inputsDir,'geom.inp'), ...
    geometryID             = 'PURECOND_DOWN', ...
    bcFilePath             = fullfile(inputsDir,'purecond_downflow.inp'), ...
    sessionDirName         = 'BESTION_ANGLECHECK_DOWN', ...
    sessionParentDir       = outputDir, ...
    overwriteSessionFiles  = false, ...
    LOGMODE                = 'BOTH');

upSet = InputSet( ...
    modelFilePath          = tempModelFile, ...
    modelID                = bestionCaseId, ...
    optionsFilePath        = fullfile(inputsDir,'options.inp'), ...
    optionsID              = 'DEFAULT', ...
    geometryFilePath       = fullfile(inputsDir,'geom.inp'), ...
    geometryID             = 'TUTORIAL1', ...
    bcFilePath             = fullfile(inputsDir,'purecond_downflow.inp'), ...
    sessionDirName         = 'BESTION_ANGLECHECK_UP', ...
    sessionParentDir       = outputDir, ...
    overwriteSessionFiles  = false, ...
    LOGMODE                = 'BOTH');

downSolver = MixtureSolver(downSet); downSolver.solve();
upSolver   = MixtureSolver(upSet);   upSolver.solve();

assert(downSolver.STATE == Solvers.SolverState.INITIALSTEPCONVERGED, ...
    'Bestion downflow comparison case did not converge.');
assert(upSolver.STATE == Solvers.SolverState.INITIALSTEPCONVERGED, ...
    'Bestion upflow comparison case did not converge.');

vfDown = downSolver.mixture(1).VF;
vfUp   = upSolver.mixture(1).VF;

assert(max(abs(vfDown-vfUp),[],'omitnan') > 1E-6, ...
    'Bestion upflow/downflow void fraction are unexpectedly identical; angle dependence may be inactive.');

figDir = fullfile(outputDir,'figures');
if ~isfolder(figDir)
    mkdir(figDir);
end

for iCase = 1:numel(plotCases)
    caseData = plotCases(iCase);
    caseId = caseData.id;
    mix = caseData.mix;
    fluid = caseData.fluid;
    geom = caseData.geom;
    twall = caseData.twall;

    z = mix.Z(:);                                                         % [m]
    q = mix.HFLUX;                                                        % [W/m^2]
    x = repmat(mix.X(),1,size(q,2));                                      % [-]
    mul = fluid.MUL(mix.liquid.H);                                        % [Pa.s]
    kl = fluid.KL(mix.liquid.H);                                          % [W/m/K]

    qprime = sum(geom.PERIM.*q,2);                                        % [W/m]
    mcondPerLength = -qprime./fluid.HFG;                                  % [kg/s/m] -q'/h_fg

    wWall = geom.RWALL.*mix.W;                                            % [kg/s] Wall-distributed flow
    ref = 4.*(1-x).*wWall./(max(1E-12,mul).*geom.PERIM);                  % [-] Film Reynolds number
    hcond = mix.HWALLCOND();                                              % [W/m^2/K] Condensation HTC
    deltaNusselt = kl./max(hcond,1E-12);                                  % [m] Nusselt film thickness estimate

    % Build wall-averaged traces for compact plotting.
    qMean = mean(q,2,'omitnan');
    refMean = mean(ref,2,'omitnan');
    deltaMean = mean(deltaNusselt,2,'omitnan');
    hcondMean = mean(hcond,2,'omitnan');
    twallMean = mean(twall,2,'omitnan');

    uLiquid = mix.liquid.U();                                             % [m/s]
    uVapor = mix.vapor.U();                                               % [m/s]
    uMix = mix.U();                                                       % [m/s]
    vf = mix.VF();                                                        % [-]

    caseFigDir = fullfile(figDir,caseId);
    if ~isfolder(caseFigDir)
        mkdir(caseFigDir);
    end

    saveLineFigure(z,mcondPerLength, 'z [m]', '-q''/h_{fg} [kg/s/m]', 'Condensation rate per length', fullfile(caseFigDir,'condensation_rate_per_length'));
    saveLineFigure(z,refMean,       'z [m]', 'Re_f [-]',             'Film Reynolds number',          fullfile(caseFigDir,'film_reynolds_number'));
    saveLineFigure(z,deltaMean,     'z [m]', '\delta [m]',          'Nusselt film thickness',        fullfile(caseFigDir,'nusselt_film_thickness'));
    saveLineFigure(z,hcondMean,     'z [m]', 'h_{cond} [W/m^2/K]',   'Condensation heat transfer coefficient', fullfile(caseFigDir,'h_cond'));
    saveLineFigure(z,twallMean,     'z [m]', 'T_w [K]',              'Wall temperature',              fullfile(caseFigDir,'t_wall'));
    saveMultiLineFigure( ...
        z, {uLiquid,uVapor,uMix}, {'Liquid','Gas','Mixture'}, ...
        'z [m]', 'Velocity [m/s]', 'Liquid, gas, and mixture velocity', ...
        fullfile(caseFigDir,'velocity_profiles'));
    saveLineFigure(z,vf,            'z [m]', 'v_f [-]',              'Void fraction',                 fullfile(caseFigDir,'vf_profile'));
    saveLineFigure(z,qMean,         'z [m]', 'HF [W/m^2]',           'Wall heat flux',                fullfile(caseFigDir,'hf_profile'));
end

disp('Downward pure-steam condensation integration test: PASSED');
disp(struct2table(results));
fprintf('Condensation figures saved in: %s\n',figDir);

clear warningCleanup

function saveLineFigure(x,y,xLabelText,yLabelText,plotTitle,fileStem)
%SAVELINEFIGURE Create, format, and save a single-line figure.

fig = figure('Visible','off');
plot(x,y,'LineWidth',1.6);
grid on;
xlabel(xLabelText);
ylabel(yLabelText);
title(plotTitle);

saveas(fig,[fileStem '.png']);
saveas(fig,[fileStem '.fig']);
close(fig);
end

function saveMultiLineFigure(x,ySeries,legendLabels,xLabelText,yLabelText,plotTitle,fileStem)
%SAVEMULTILINEFIGURE Create, format, and save a multi-line figure.

fig = figure('Visible','off');
hold on;
for i = 1:numel(ySeries)
    plot(x,ySeries{i},'LineWidth',1.6);
end
hold off;
grid on;
legend(legendLabels,'Location','best');
xlabel(xLabelText);
ylabel(yLabelText);
title(plotTitle);

saveas(fig,[fileStem '.png']);
saveas(fig,[fileStem '.fig']);
close(fig);
end