%% Downward pure-steam condensation integration test (single-case)
%
% Simplified integration test that runs only PURECOND_NUSSELT_ZIVI and
% generates/saves all post-processing figures.

clearvars

scriptPath = mfilename('fullpath');
scriptDir = fileparts(scriptPath);
openstreamRoot = fileparts(scriptDir);
addpath(openstreamRoot);

import Inputs.*
import Solvers.*
import Solvers.Mixture.*

inputsDir = fullfile(openstreamRoot,'inputs');
outputDir = fullfile(scriptDir,'outputs','integration_downward_condensation_v2');

if isfolder(outputDir)
    rmdir(outputDir,'s');
end
mkdir(outputDir);

tempModelFile = fullfile(outputDir,'models_downward_condensation_v2.inp');
fclose(fopen(tempModelFile,'w'));
tempBcFile = fullfile(outputDir,'purecond_downflow_twallin_v2.inp');

% Build a single model case: PURECOND_NUSSELT_ZIVI
caseId = 'PURECOND_NUSSELT_ZIVI';
Model.writeInputFile( ...
    tempModelFile, caseId, ...
    'NNODES', 100, ...
    'FLUID', 'WATER', ...
    'VOID', 'ZIVI', ...
    'CONDHTM', 'NUSSELT', ...
    'TPFM', 'SLIP', ...
    'MOMENTLIQUID', 'FULL', ...
    'MOMENTGAS', 'FULL', ...
    'MOMENTFILM', 'FULL', ...
    'MOMENTDROP', 'FULL');

warningState = warning;
warning('off','all');
warningCleanup = onCleanup(@() warning(warningState));

% Build a local BC file that explicitly exercises TWALLIN mode.
bcFid = fopen(tempBcFile,'w');
assert(bcFid > 0, 'Could not create temporary BC file: %s', tempBcFile);
fprintf(bcFid,'TIME          ! Time [s]                                                   > 0.0\n');
fprintf(bcFid,'PRESSURE      ! System pressure [Pa]                                       > 300000\n');
fprintf(bcFid,'HIN           ! Inlet enthalpy [J/kg]                                      > 2.72E6\n');
fprintf(bcFid,'MFLOW         ! Mass flow rate [kg/s]                                      > 0.05\n');
fprintf(bcFid,'BCTYPE        ! Wall thermal boundary-condition type [-]                   > TWALLIN\n');
fprintf(bcFid,'TWALL         ! Wall temperature distribution(s) [K]                       > 390 385\n');
fprintf(bcFid,'POWER         ! Total power [W]                                            > -31000\n');
fprintf(bcFid,'WMESH         ! Relative power node size distribution [m]                  > 1.5 1.5\n');
fprintf(bcFid,'WPOWER        ! Relative power distribution(s) [-]                         > 0.5 0.5\n');
fprintf(bcFid,'END\n');
fclose(bcFid);

% Solve case
inputSet = InputSet( ...
    modelFilePath          = tempModelFile, ...
    modelID                = caseId, ...
    optionsFilePath        = fullfile(inputsDir,'options.inp'), ...
    optionsID              = 'DEFAULT', ...
    geometryFilePath       = fullfile(inputsDir,'geom.inp'), ...
    geometryID             = 'PURECOND_DOWN', ...
    bcFilePath             = tempBcFile, ...
    sessionParentDir       = outputDir, ...
    overwriteSessionFiles  = false, ...
    LOGMODE                = 'BOTH');

assert(inputSet.bc(1).POWER < 0, ...
    'Boundary condition POWER must be negative for condensation test.');
assert(all(inputSet.bc(1).WPOWER(:) >= 0), ...
    'Boundary condition WPOWER must stay non-negative.');
assert(string(inputSet.bc(1).BCTYPE) == "TWALLIN", ...
    'Boundary condition BCTYPE must be TWALLIN for this test.');

mixSolver = MixtureSolver(inputSet);
mixSolver.solve();

assert(mixSolver.STATE == Solvers.SolverState.INITIALSTEPCONVERGED, ...
    'Case %s did not converge to INITIALSTEPCONVERGED.', caseId);

% For single-time inputs, steady-state converges in mixtureInit and
% transient solve is skipped. Use the converged state accordingly.
if numel(mixSolver.mixture) < 2
    mix = mixSolver.mixtureInit(end);
else
    mix = mixSolver.mixture(1);
end
fluid = mixSolver.fluid;
geom = inputSet.geometry;

q = mix.HFLUX;
x = repmat(mix.X(),1,size(q,2));
idxCond = q < 0 & x > 0 & x < 1;

assert(any(q(:) < 0), ...
    'Case %s has no negative wall heat flux.', caseId);
assert(any(idxCond(:)), ...
    'Case %s did not produce a two-phase condensation region (0 < x < 1 with q < 0).', caseId);

[twall,~] = mix.TWALL();
[hwall,~] = mix.HWALL(twall);
condResidual = abs(twall(idxCond) - (fluid.TSAT + q(idxCond)./hwall(idxCond)));
maxCondResidual = max(condResidual,[],'omitnan');

assert(all(isfinite(twall(:))) && isreal(twall), ...
    'Case %s produced non-finite or complex wall temperature.', caseId);
assert(maxCondResidual < 1E-2, ...
    'Case %s violates TWALL condensation relation. Max residual: %.3e K', caseId, maxCondResidual);

% TWALLIN-specific convergence check using dq inner-iteration metric.
dqMax = max(mix.ITR.DQ,[],'omitnan');
assert(dqMax < inputSet.options.ERRORQ, ...
    'Case %s did not satisfy dq convergence: max dq = %.3e W/m^2, target = %.3e W/m^2', ...
    caseId, dqMax, inputSet.options.ERRORQ);

% Post-processing metrics
z = mix.Z(:);                                                             % [m]
mul = fluid.MUL(mix.liquid.H);                                            % [Pa.s]
kl = fluid.KL(mix.liquid.H);                                              % [W/m/K]

qprime = sum(geom.PERIM.*q,2);                                            % [W/m]
mcondPerLength = -qprime./fluid.HFG;                                      % [kg/s/m] -q'/h_fg

wWall = geom.RWALL.*mix.W;                                                % [kg/s] Wall-distributed flow
ref = 4.*(1-x).*wWall./(max(1E-12,mul).*geom.PERIM);                      % [-] Film Reynolds number
hcond = mix.HWALLCOND();                                                  % [W/m^2/K] Condensation HTC
deltaNusselt = kl./max(hcond,1E-12);                                      % [m] Nusselt film thickness estimate
vf = mix.VF();                                                            % [-]
deltaFromvf = geom.PERIM/pi/2 * (1 - sqrt(vf));                           % [m] Film thickness estimate from void fraction


% Build wall-averaged traces for compact plotting.
qMean = mean(q,2,'omitnan');
refMean = mean(ref,2,'omitnan');
deltaMean = mean(deltaNusselt,2,'omitnan');
deltaFromvfMean = mean(deltaFromvf,2,'omitnan');
hcondMean = mean(hcond,2,'omitnan');
twallMean = mean(twall,2,'omitnan');

uLiquid = mix.liquid.U();                                                 % [m/s]
uVapor = mix.vapor.U();                                                   % [m/s]
uMix = mix.U();                                                           % [m/s]


% Save figures
figDir = fullfile(outputDir,'figures',caseId);
if ~isfolder(figDir)
    mkdir(figDir);
end

saveLineFigure(z,mcondPerLength, 'z [m]', '-q''/h_{fg} [kg/s/m]', 'Condensation rate per length', fullfile(figDir,'condensation_rate_per_length'));
saveLineFigure(z,refMean,        'z [m]', 'Re_f [-]',             'Film Reynolds number',          fullfile(figDir,'film_reynolds_number'));
saveLineFigure(z,deltaMean,      'z [m]', '\delta [m]',          'Nusselt film thickness',        fullfile(figDir,'nusselt_film_thickness'));
saveLineFigure(z,hcondMean,      'z [m]', 'h_{cond} [W/m^2/K]',   'Condensation heat transfer coefficient', fullfile(figDir,'h_cond'));
saveLineFigure(z,twallMean,      'z [m]', 'T_w [K]',              'Wall temperature',              fullfile(figDir,'t_wall'));
saveMultiLineFigure( ...
    z, {uLiquid,uVapor,uMix}, {'Liquid','Gas','Mixture'}, ...
    'z [m]', 'Velocity [m/s]', 'Liquid, gas, and mixture velocity', ...
    fullfile(figDir,'velocity_profiles'));
saveLineFigure(z,vf,             'z [m]', 'v_f [-]',              'Void fraction',                 fullfile(figDir,'vf_profile'));
saveLineFigure(z,deltaFromvf,    'z [m]', '\delta [m]',          'Film thickness from void fraction', fullfile(figDir,'delta_from_vf'));
saveLineFigure(z,qMean,          'z [m]', 'HF [W/m^2]',           'Wall heat flux',                fullfile(figDir,'hf_profile'));

disp('Downward pure-steam condensation integration test v2: PASSED');
fprintf('Case: %s\n', caseId);
fprintf('Max condensation TWALL residual: %.3e K\n', maxCondResidual);
fprintf('Max dq inner-iteration residual: %.3e W/m^2\n', dqMax);
fprintf('Condensation figures saved in: %s\n', figDir);

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
