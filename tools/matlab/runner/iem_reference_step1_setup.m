% iem_reference_step1_setup  Runner job 1: environment checks and the runIEM_HH model setup on Harvey 1.03d.
%
% For the matlab-agent-runner (MATLAB R2024b, COBRA Toolbox, Gurobi). Inputs are read from
% /inputs/iem_reference (Harvey_1_03d.mat, load_reference_harvey.m); outputs go to output_dir.
% Writes run_info_step1.txt and matlab_setup_bounds_Harvey_1_03d.mat (rxns, lb, ub after the
% diet and physiological setup of runIEM_HH.m). Takes minutes.

refKitDir = '/inputs/iem_reference';
addpath(refKitDir);
refLog = fopen(fullfile(output_dir, 'run_info_step1.txt'), 'w');
fprintf(refLog, 'IEM reference run, step 1 (checks and setup bounds)\n');
fprintf(refLog, 'started %s\n', datestr(now, 'yyyy-mm-dd HH:MM:SS'));
fprintf(refLog, 'MATLAB %s on %s\n', version, computer);
global CBT_LP_SOLVER
fprintf(refLog, 'LP solver %s\n', CBT_LP_SOLVER);

% model file
refFid = fopen(fullfile(refKitDir, 'Harvey_1_03d.mat'), 'r');
assert(refFid >= 0, 'Harvey_1_03d.mat is missing from /inputs/iem_reference');
refBytes = fread(refFid, Inf, '*uint8');
fclose(refFid);
refDigest = java.security.MessageDigest.getInstance('SHA-256');
refDigest.update(refBytes);
refHash = lower(reshape(dec2hex(typecast(refDigest.digest(), 'uint8'), 2)', 1, []));
clear refBytes refDigest
fprintf(refLog, 'model sha256 %s\n', refHash);
assert(strcmp(refHash, '10e6cb6d02736fb88ae266e0c901e90766176bb80db2753e694e3a9330368c9c'), ...
    'Harvey_1_03d.mat is not the expected file');

% Toolbox files the protocol depends on: git blob SHA-1 computed here (no git needed), compared with commit 67c790d
refCbtDir = fileparts(which('initCobraToolbox'));
fprintf(refLog, 'COBRA Toolbox at %s\n', refCbtDir);
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
    refFid = fopen(refPath, 'r');
    if refFid < 0
        fprintf(refLog, 'MISSING   %s\n', refFiles{refI, 1});
        refNDiffer = refNDiffer + 1;
        continue
    end
    refBytes = fread(refFid, Inf, '*uint8');
    fclose(refFid);
    refDigest = java.security.MessageDigest.getInstance('SHA-1');
    refDigest.update([uint8(sprintf('blob %d', numel(refBytes)))'; uint8(0); refBytes(:)]);
    refBlob = lower(reshape(dec2hex(typecast(refDigest.digest(), 'uint8'), 2)', 1, []));
    if strcmp(refBlob, refFiles{refI, 2})
        fprintf(refLog, 'same      %s\n', refFiles{refI, 1});
    else
        fprintf(refLog, 'DIFFERENT %s (%s)\n', refFiles{refI, 1}, refBlob);
        refNDiffer = refNDiffer + 1;
    end
end
clear refBytes refDigest
fprintf(refLog, '%d of %d checked Toolbox files differ from commit 67c790d\n', refNDiffer, size(refFiles, 1));

% model setup only, as at the top of runIEM_HH.m (Harvey branch)
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
save(fullfile(output_dir, 'matlab_setup_bounds_Harvey_1_03d.mat'), 'rxns', 'lb', 'ub', '-v7');
fprintf(refLog, 'setup bounds saved after %.1f s\n', toc(refStart));
fprintf(refLog, 'finished %s\n', datestr(now, 'yyyy-mm-dd HH:MM:SS'));
fclose(refLog);
disp(jsonencode(struct('status', 'step1_done', 'toolbox_files_differing', refNDiffer, ...
    'n_rxns', numel(rxns), 'setup_seconds', toc(refStart))));
clear male sex IndividualParameters Diet rxns lb ub
