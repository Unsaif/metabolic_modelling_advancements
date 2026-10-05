% run_iem_reference_harvey_1_03d  Reference run of runIEM_HH.m (COBRA Toolbox) on Harvey 1.03d.
%
% Purpose: compare the MATLAB protocol with the Python port in this repository
% (gembench/wbm_constraints.py and scripts/run_wbm_iem.py, run v0.3).
% A copy of runIEM_HH.m is run with two changes and nothing else:
%   1. the model is loaded from ./Harvey_1_03d.mat (runIEM_HH loads the newest Harvey on the path);
%   2. line 9, "edit% This script ...", loses its stray "edit" command (it would open the Editor).
%
% Before running, in MATLAB:
%   initCobraToolbox(false)
%   changeCobraSolver('ibm_cplex', 'LP')   % or 'gurobi' / 'mosek': the LP solver the lab normally uses
%   cd <the folder holding this script, load_reference_harvey.m and Harvey_1_03d.mat>
% Then type at the command prompt (the script must run in the base workspace, because runIEM_HH reads
% its results back with evalin('base', ...)):
%   run_iem_reference_harvey_1_03d
%
% Outputs in ./out (please send back run_info.txt and the three matlab_*.mat files):
%   run_info.txt                            Toolbox commit, file checks, solver, MATLAB version, timings
%   matlab_setup_bounds_Harvey_1_03d.mat    rxns, lb, ub after the diet and physiological setup (minutes)
%   matlab_global_bounds_Harvey_1_03d.mat   rxns, lb, ub after the unified reaction constraints (modelO)
%   matlab_iem_results_Harvey_1_03d.mat     IEMSol_* of all IEMs, Table_IEM, Accuracy, ... (hours)
%   Results_IEM_Harvey_1_03.mat             the workspace runIEM_HH saves itself (large; not needed)

refKitDir = fileparts(mfilename('fullpath'));
refOutDir = fullfile(refKitDir, 'out');
if ~exist(refOutDir, 'dir')
    mkdir(refOutDir);
end
refModelFile = fullfile(refKitDir, 'Harvey_1_03d.mat');
refModelSha256 = '10e6cb6d02736fb88ae266e0c901e90766176bb80db2753e694e3a9330368c9c';
refLog = fopen(fullfile(refOutDir, 'run_info.txt'), 'w');
fprintf(refLog, 'Reference run of runIEM_HH.m on Harvey 1.03d\n');
fprintf(refLog, 'started %s\n', datestr(now, 'yyyy-mm-dd HH:MM:SS'));
fprintf(refLog, 'MATLAB %s on %s\n', version, computer);

% --- solver ---------------------------------------------------------------------------------------
global CBT_LP_SOLVER
if isempty(CBT_LP_SOLVER)
    fclose(refLog);
    error('No LP solver is set: run initCobraToolbox(false) and changeCobraSolver(''ibm_cplex'', ''LP'') first.');
end
fprintf(refLog, 'LP solver %s\n', CBT_LP_SOLVER);

% --- model file ------------------------------------------------------------------------------------
refFid = fopen(refModelFile, 'r');
if refFid < 0
    fclose(refLog);
    error('Harvey_1_03d.mat was not found next to this script.');
end
refBytes = fread(refFid, Inf, '*uint8');
fclose(refFid);
refDigest = java.security.MessageDigest.getInstance('SHA-256');
refDigest.update(refBytes);
refHash = lower(reshape(dec2hex(typecast(refDigest.digest(), 'uint8'), 2)', 1, []));
clear refBytes refDigest
fprintf(refLog, 'model %s sha256 %s\n', refModelFile, refHash);
if ~strcmp(refHash, refModelSha256)
    fclose(refLog);
    error('Harvey_1_03d.mat is not the expected file (sha256 %s).', refModelSha256);
end

% --- Toolbox version: files the protocol depends on, compared with commit 67c790d --------------------
refCbtDir = fileparts(which('initCobraToolbox'));
[refStatus, refOutput] = system(sprintf('git -C "%s" rev-parse HEAD', refCbtDir));
if refStatus == 0
    fprintf(refLog, 'COBRA Toolbox %s at commit %s\n', refCbtDir, strtrim(refOutput));
else
    fprintf(refLog, 'COBRA Toolbox %s (commit unknown)\n', refCbtDir);
end
refFiles = {
    'src/analysis/wholeBody/PSCMToolbox/runIEM_HH.m', '31d3d27cc297e4e59feefe2c9b2018d5e3688e58'
    'src/analysis/wholeBody/PSCMToolbox/checkIEM_WBM.m', 'edc042a677f753383f6fae6514209bb849c1b22e'
    'src/analysis/wholeBody/PSCMToolbox/optimizeWBModel.m', 'c625a1000d63c251e9e61344f4b48a9b5f30fb18'
    'src/analysis/wholeBody/PSCMToolbox/io/loadPSCMfile.m', 'c300ca139e67d32e5dc46572f796e21f225f68c1'
    'src/analysis/wholeBody/PSCMToolbox/io/OrganLists.m', '8cef0e8b6e4a31b0c0230456a058fa84f3b38c2d'
    'src/analysis/wholeBody/PSCMToolbox/setConstraints/physiologicalConstraintsHMDBbased.m', 'f31d50dce96dc2e56e29aa3f56ef70a3053622ca'
    'src/analysis/wholeBody/PSCMToolbox/setConstraints/setDietConstraints.m', 'bbbd039f6a43a0c655920a69f8f0459db2bf4433'
    'src/analysis/wholeBody/PSCMToolbox/setConstraints/standardPhysiolDefaultParameters.m', 'cebb9d51a87c475c2afd863985bd6534626e6c19'
    'src/analysis/wholeBody/PSCMToolbox/setConstraints/diets/EUAverageDietNew.m', 'ef3b29ac2e18027845ecc55a712f057ca8fa4a9b'
    'src/analysis/wholeBody/PSCMToolbox/setConstraints/organWeight/getOrganWeightFraction.m', '163c7731e3d9e1a2997b416f356499b152f264b3'
    'src/analysis/wholeBody/PSCMToolbox/setConstraints/inputData/NormalBloodConcExtractedHMDB.txt', '5973b202665a512fd061c79949c731f39f2451f9'
    'src/analysis/wholeBody/PSCMToolbox/setConstraints/inputData/NormalCSFConcExtractedHMDB.txt', 'fd5b8380ecf42633ab5d7b5b7340286971c4d3cf'
    'src/analysis/wholeBody/PSCMToolbox/setConstraints/inputData/NormalUrineConcExtractedHMDB.txt', 'd5d2f10656de98b41bcd7b33bfffeede779e5668'
    'src/analysis/wholeBody/PSCMToolbox/setConstraints/inputData/16_01_26_BloodFlowRatesPercentages.xlsx', 'd7473d92a43d74d1835d72b67c85a1776f0fcebb'
    'src/analysis/wholeBody/PSCMToolbox/setConstraints/inputData/Numbers_organWeightData_fromxls_16_01_29_OrganWeigths.mat', '953c2ab4bbf10a8845b148d076679036975547e6'
    'src/analysis/wholeBody/PSCMToolbox/hostMicrobeInteraction/AGORAEssentialMetabolites.m', '2910f3cbb532aa3a43495fe0b1e0a5e78767c823'
    'src/analysis/FBA/changeRxnBounds.m', 'e7d1a25697bdf227e16ddcc772ec52c8e2a01efe'
    'src/reconstruction/refinement/addDemandReaction.m', 'd4a4dfb15eb1e8aa35e2a7dd2f3ad1b21b94837f'
    'src/base/solvers/solveCobraLP.m', '9b1a2a2eb746d54374b9406c7985c6dd673794b2'
    };
refNDiffer = 0;
for refI = 1:size(refFiles, 1)
    refPath = fullfile(refCbtDir, refFiles{refI, 1});
    if ~exist(refPath, 'file')
        fprintf(refLog, 'MISSING   %s\n', refFiles{refI, 1});
        refNDiffer = refNDiffer + 1;
        continue
    end
    [refStatus, refOutput] = system(sprintf('git hash-object --no-filters "%s"', refPath));
    refOutput = strtrim(refOutput);
    if refStatus ~= 0
        fprintf(refLog, 'UNCHECKED %s (git is not available)\n', refFiles{refI, 1});
    elseif strcmp(refOutput, refFiles{refI, 2})
        fprintf(refLog, 'same      %s\n', refFiles{refI, 1});
    else
        fprintf(refLog, 'DIFFERENT %s\n', refFiles{refI, 1});
        refNDiffer = refNDiffer + 1;
    end
end
fprintf(refLog, '%d of %d checked Toolbox files differ from commit 67c790d\n', refNDiffer, size(refFiles, 1));
if refNDiffer > 0
    warning(['%d COBRA Toolbox files differ from commit 67c790d (see out/run_info.txt). ', ...
             'The run continues; the comparison is then with your Toolbox version.'], refNDiffer);
end

% --- step 1: model setup only, as at the top of runIEM_HH.m (Harvey branch) --------------------------
refStart = tic;
male = load_reference_harvey();
sex = male.sex;
standardPhysiolDefaultParameters;
male = physiologicalConstraintsHMDBbased(male, IndividualParameters);
EUAverageDietNew;
male = setDietConstraints(male, Diet);
rxns = male.rxns;
lb = male.lb;
ub = male.ub;
save(fullfile(refOutDir, 'matlab_setup_bounds_Harvey_1_03d.mat'), 'rxns', 'lb', 'ub', '-v7');
fprintf(refLog, 'step 1 (setup bounds) saved after %.0f s\n', toc(refStart));

% --- step 2: the full protocol from a minimally changed copy of runIEM_HH.m ------------------------------
refSource = fullfile(refCbtDir, 'src', 'analysis', 'wholeBody', 'PSCMToolbox', 'runIEM_HH.m');
refFid = fopen(refSource, 'r');
refText = char(fread(refFid, Inf, '*uint8')');
fclose(refFid);
refOld = [char([13 10]) '    male = loadPSCMfile(modelName);'];
refNew = [char([13 10]) '    male = load_reference_harvey();'];
if numel(strfind(refText, refOld)) ~= 1
    refOld = [char(10) '    male = loadPSCMfile(modelName);'];
    refNew = [char(10) '    male = load_reference_harvey();'];
end
if numel(strfind(refText, refOld)) ~= 1
    fclose(refLog);
    error('Could not find the Harvey load line in runIEM_HH.m exactly once.');
end
refText = strrep(refText, refOld, refNew);
refText = strrep(refText, [char(10) 'edit% This script'], [char(10) '% This script']);
refCopy = fullfile(refOutDir, 'runIEM_HH_ref.m');
refFid = fopen(refCopy, 'w');
fwrite(refFid, uint8(refText), 'uint8');
fclose(refFid);
fprintf(refLog, 'runIEM_HH copy written to %s\n', refCopy);

clearvars -except ref*
addpath(refKitDir);
modelName = 'Harvey';
resultsPath = [refOutDir filesep];
run(refCopy);

fprintf(refLog, 'step 2 (runIEM_HH) finished after %.0f s\n', toc(refStart));
refVars = who('-regexp', '^IEMSol_|^Table_IEM$|^Table_sum_results$|^Accuracy$|^Precision$|^FalseDiscoveryRate$|^NumDiseases$|^NumBiomarkers$');
save(fullfile(refOutDir, 'matlab_iem_results_Harvey_1_03d.mat'), refVars{:}, '-v7');
if exist('modelO', 'var')
    rxns = modelO.rxns;
    lb = modelO.lb;
    ub = modelO.ub;
    save(fullfile(refOutDir, 'matlab_global_bounds_Harvey_1_03d.mat'), 'rxns', 'lb', 'ub', '-v7');
end
if exist('Accuracy', 'var')
    fprintf(refLog, 'Accuracy as defined at the end of runIEM_HH.m: %.4f\n', Accuracy);
end
fprintf(refLog, 'finished %s\n', datestr(now, 'yyyy-mm-dd HH:MM:SS'));
fclose(refLog);
fprintf('Done. Please send back the files in %s (run_info.txt and matlab_*.mat).\n', refOutDir);
