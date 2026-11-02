function result = piconv_mixture_case(runId, opts)
%PICONV_MIXTURE_CASE Simulate one PICON-V pure-steam run with the mixture solver
%
% Runs a PICON-V downward condensation test (Ahn et al., 40 mm tube,
% 3 m cooled length, 10 measurement sections of 0.3 m) with the measured
% inner wall temperatures as TWALLIN boundary condition, for each
% condensation heat transfer model (CONDHTM) of the mixture solver, and
% stores the predicted local wall heat flux, quality and HTC next to the
% measured data.
%
% Inputs:
%
% - runId          — PICON-V run number (Data#)
% - opts.dataFile  — Experimental data workbook
% - opts.outputDir — Output directory for session files and the result file
% - opts.NNODES    — Number of axial nodes
% - opts.CONDHTM   — Condensation HTC models to simulate
% - opts.VOID      — Void fraction model
%
% Returns:
%
% - result — Struct with experimental data and model predictions, also saved to <outputDir>/PICONV_<runId>.mat

arguments
    runId (1,1) double {mustBeInteger,mustBePositive}
    opts.dataFile  (1,1) string = fullfile(fileparts(mfilename('fullpath')),'PICONV_data_Ahn.xlsx')
    opts.outputDir (1,1) string = fullfile(fileparts(mfilename('fullpath')),'outputs','validation_piconv_mixture')
    opts.NNODES    (1,1) double {mustBeInteger,mustBePositive} = 60
    opts.CONDHTM   (1,:) string = ["NUSSELT","SHAH","KIMMUDAWAR"]
    opts.VOID      (1,1) string = "ZIVI"
end

scriptDir = fileparts(mfilename('fullpath'));
openstreamRoot = fileparts(scriptDir);
addpath(openstreamRoot);

import Inputs.*
import Solvers.*

caseDir = fullfile(opts.outputDir,sprintf('PICONV_%03d',runId));
if ~isfolder(caseDir), mkdir(caseDir); end

warningState = warning;
warning('off','all');
warningCleanup = onCleanup(@() warning(warningState));

%% Experimental data

data = readPiconv(opts.dataFile, runId);

D     = data.D;                                                            % [m]
AREA  = pi*D^2/4;                                                          % [m^2]
PERIM = pi*D;                                                              % [m]
NSEC  = numel(data.L);                                                     % [-] Measurement sections
DSEC  = 2*data.L(1);                                                       % [m] Section length
LTOT  = NSEC*DSEC;                                                         % [m] Cooled length

coolprop = CoolPropWrapper.CoolPropWrapper('WATER');
HIN = coolprop.enthalpy('P',data.P,'Q',data.xin);                          % [J/kg] Inlet enthalpy from measured inlet quality

%% Input files

geomFile = fullfile(caseDir,'geom.inp');
fid = fopen(geomFile,'w');
fprintf(fid,'ID        ! Channel ID                     > PICONV\n');
fprintf(fid,'LENGTH    ! Axial length [m]               > %.6g\n',LTOT);
fprintf(fid,'AREA      ! Coolant area [m^2]             > %.8g\n',AREA);
fprintf(fid,'PERIM     ! Perimeters [m]                 > %.8g\n',PERIM);
fprintf(fid,'ANGLE     ! Angle [deg]                    > 180\n');
fprintf(fid,'END\n');
fclose(fid);

bcFile = fullfile(caseDir,'bc.inp');
fid = fopen(bcFile,'w');
fprintf(fid,'TIME          ! Time [s]                                                   > 0.0\n');
fprintf(fid,'PRESSURE      ! System pressure [Pa]                                       > %.8g\n',data.P);
fprintf(fid,'HIN           ! Inlet enthalpy [J/kg]                                      > %.8g\n',HIN);
fprintf(fid,'MFLOW         ! Mass flow rate [kg/s]                                      > %.8g\n',data.G*AREA);
fprintf(fid,'BCTYPE        ! Wall thermal boundary-condition type [-]                   > TWALLIN\n');
fprintf(fid,'TWALL         ! Measured inner wall temperature [K]                        > %s\n',sprintf('%.6g ',data.Twi+273.15));
fprintf(fid,'POWER         ! Total power [W] (unused with TWALLIN)                      > -1\n');
fprintf(fid,'WMESH         ! Measurement section lengths [m]                            > %s\n',sprintf('%.6g ',DSEC*ones(1,NSEC)));
fprintf(fid,'WPOWER        ! Relative power distribution(s) [-]                         > %s\n',sprintf('%g ',ones(1,NSEC)));
fprintf(fid,'END\n');
fclose(fid);

optionsFile = fullfile(caseDir,'options.inp');
fid = fopen(optionsFile,'w');
fprintf(fid,'ID            ! Options ID                                                 > PICONV\n');
fprintf(fid,'RELAXQ        ! Relaxation factor for TWALL-driven wall heat flux update [-] > 0.5\n');
fprintf(fid,'ERRORQ        ! Heat flux error target in inner iterations [W/m^2]        > 1E-1\n');
fprintf(fid,'SSTSTEP       ! Steady-state pseudo-time step [s] (slow condensate flow)   > 1E4\n');
fprintf(fid,'END\n');
fclose(fid);

modelFile = fullfile(caseDir,'models.inp');
fclose(fopen(modelFile,'w'));
common = {'NNODES',opts.NNODES,'FLUID','WATER','VOID',char(opts.VOID),'TPFM','SLIP', ...
    'MOMENTLIQUID','FULL','MOMENTGAS','FULL'};
for m = opts.CONDHTM
    Model.writeInputFile(modelFile,char(m),common{:},'CONDHTM',char(m));
end

%% Simulations

result = struct('runId',runId,'data',data,'z',[],'models',struct());

for m = opts.CONDHTM
    inputSet = InputSet( ...
        modelFilePath          = modelFile, ...
        modelID                = char(m), ...
        optionsFilePath        = optionsFile, ...
        optionsID              = 'PICONV', ...
        geometryFilePath       = geomFile, ...
        geometryID             = 'PICONV', ...
        bcFilePath             = bcFile, ...
        sessionParentDir       = caseDir, ...
        overwriteSessionFiles  = true, ...
        LOGMODE                = 'LOGTOFILEONLY');
    [result.models.(m), z] = runMixture(inputSet);
    if ~isempty(z), result.z = z; end
end

save(fullfile(opts.outputDir,sprintf('PICONV_%03d.mat',runId)),'result');

end

function [res, z] = runMixture(inputSet)
% Solve the mixture model and extract wall heat flux, quality and HTC

res = struct('ok',false,'state',"",'q',[],'x',[],'xeq',[],'h',[],'P',[],'dqmax',NaN,'time',NaN,'msg',"");
z = [];
try
    t0 = tic;
    mixSolver = Solvers.Mixture.MixtureSolver(inputSet);
    mixSolver.solve();
    if numel(mixSolver.mixture) < 2
        mix = mixSolver.mixtureInit(end);
    else
        mix = mixSolver.mixture(1);
    end
    res.time  = toc(t0);
    res.state = string(mixSolver.STATE);
    res.ok    = mixSolver.STATE == Solvers.SolverState.INITIALSTEPCONVERGED;
    res.q     = mix.HFLUX(:,1);
    res.x     = mix.X();
    res.xeq   = mix.XEQ();
    res.h     = mix.HWALLCOND();
    res.P     = mix.P(:);
    res.dqmax = max(mix.ITR.DQ,[],'omitnan');
    z         = mixSolver.Z;
catch ME
    res.msg = string(ME.message);
end
end

function data = readPiconv(dataFile, runId)
% Read the summary and local measurements of one PICON-V run

S = readcell(dataFile,'Sheet','PICONV_summary (118)');
S = S(4:end,:);
row = find(cellfun(@(v) isequal(v,runId), S(:,2)),1);
assert(~isempty(row),'PICON-V run %d not found in %s.',runId,dataFile);
assert(S{row,3} == 0,'PICON-V run %d is not a pure-steam test.',runId);

data.runId = runId;
data.D     = S{row,6};                                                     % [m]
data.P     = S{row,10}*1E3;                                                % [Pa] Measured pressure (kPa in workbook)
data.G     = S{row,11};                                                    % [kg/m^2/s]
data.xin   = S{row,13};                                                    % [-]
data.xout  = S{row,14};                                                    % [-]
data.Qw    = S{row,15};                                                    % [W] Total heat (labeled kW in workbook, values in W)

L = readcell(dataFile,'Sheet','Local data(1180)');
L = L(4:end,:);
rows = cellfun(@(v) isequal(v,runId), L(:,2));
loc  = cell2mat(L(rows,[5 7 8 13 15 16 17 20]));
loc  = sortrows(loc,1);

data.node     = loc(:,1);
data.L        = loc(:,2);                                                  % [m] Section center
data.valid    = loc(:,3) == 1;                                             % [-] Valid two-phase measurement
data.x        = loc(:,4);                                                  % [-] Local quality
data.Tsat     = loc(:,5);                                                  % [degC]
data.HF       = loc(:,6);                                                  % [W/m^2]
data.HTC      = loc(:,7);                                                  % [W/m^2/K]
data.Twi      = loc(:,8);                                                  % [degC] Inner wall temperature
end
