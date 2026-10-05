% iem_reference_step2_runiem  Runner job 2: runIEM_HH.m on Harvey 1.03d (all 57 IEMs). Takes hours.
%
% Runs a copy of the Toolbox's runIEM_HH.m with two changes and nothing else: the model is loaded from
% /inputs/iem_reference/Harvey_1_03d.mat (runIEM_HH would load the newest Harvey on the path), and the
% stray "edit" command on line 9 is removed. Outputs in output_dir:
%   run_info_step2.txt, runIEM_HH_ref.m (the executed copy),
%   matlab_iem_results_Harvey_1_03d.mat (IEMSol_*, Table_IEM, Accuracy, ...),
%   matlab_global_bounds_Harvey_1_03d.mat (bounds of modelO after the unified reaction constraints),
%   Results_IEM_Harvey_1_03.mat (the workspace runIEM_HH saves itself; large).

refKitDir = '/inputs/iem_reference';
addpath(refKitDir);
refLog = fopen(fullfile(output_dir, 'run_info_step2.txt'), 'w');
fprintf(refLog, 'IEM reference run, step 2 (runIEM_HH)\nstarted %s\n', datestr(now, 'yyyy-mm-dd HH:MM:SS'));
global CBT_LP_SOLVER
fprintf(refLog, 'LP solver %s\n', CBT_LP_SOLVER);
refCbtDir = fileparts(which('initCobraToolbox'));
refSource = fullfile(refCbtDir, 'src', 'analysis', 'wholeBody', 'PSCMToolbox', 'runIEM_HH.m');
refFid = fopen(refSource, 'r');
assert(refFid >= 0, 'runIEM_HH.m not found in the Toolbox');
refText = char(fread(refFid, Inf, '*uint8')');
fclose(refFid);
refOld = [char([13 10]) '    male = loadPSCMfile(modelName);'];
refNew = [char([13 10]) '    male = load_reference_harvey();'];
if numel(strfind(refText, refOld)) ~= 1
    refOld = [char(10) '    male = loadPSCMfile(modelName);'];
    refNew = [char(10) '    male = load_reference_harvey();'];
end
assert(numel(strfind(refText, refOld)) == 1, 'Could not find the Harvey load line in runIEM_HH.m exactly once');
refText = strrep(refText, refOld, refNew);
refText = strrep(refText, [char(10) 'edit% This script'], [char(10) '% This script']);
refCopy = fullfile(output_dir, 'runIEM_HH_ref.m');
refFid = fopen(refCopy, 'w');
fwrite(refFid, uint8(refText), 'uint8');
fclose(refFid);
clear refText refOld refNew refFid refSource
fprintf(refLog, 'executing %s\n', refCopy);

refStart = tic;
modelName = 'Harvey';
resultsPath = [output_dir filesep];
run(refCopy);

fprintf(refLog, 'runIEM_HH finished after %.0f s\n', toc(refStart));
refVars = who('-regexp', '^IEMSol_|^Table_IEM$|^Table_sum_results$|^Accuracy$|^Precision$|^FalseDiscoveryRate$|^NumDiseases$|^NumBiomarkers$');
save(fullfile(output_dir, 'matlab_iem_results_Harvey_1_03d.mat'), refVars{:}, '-v7');
if exist('modelO', 'var')
    rxns = modelO.rxns;
    lb = modelO.lb;
    ub = modelO.ub;
    save(fullfile(output_dir, 'matlab_global_bounds_Harvey_1_03d.mat'), 'rxns', 'lb', 'ub', '-v7');
end
if exist('Accuracy', 'var')
    fprintf(refLog, 'Accuracy as defined at the end of runIEM_HH.m: %.4f\n', Accuracy);
end
fprintf(refLog, 'finished %s\n', datestr(now, 'yyyy-mm-dd HH:MM:SS'));
fclose(refLog);
disp(jsonencode(struct('status', 'step2_done', 'n_iem_results', numel(refVars), 'seconds', toc(refStart))));
