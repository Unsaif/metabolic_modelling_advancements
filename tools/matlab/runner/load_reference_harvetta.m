function model = load_reference_harvetta()
% Load ./Harvetta_1_03d.mat as loadPSCMfile + correctWBMfields do, without searching the MATLAB path.
%
% USAGE:
%    model = load_reference_harvetta()
%
% Used by the Harvetta runner jobs in place of loadPSCMfile('Harvetta'), which loads the newest
% Harvey version found in 2020_WholeBodyModelling/Data rather than version 1.03d.
kitDir = fileparts(mfilename('fullpath'));
S = load(fullfile(kitDir, 'Harvetta_1_03d.mat'));
names = fieldnames(S);
model = S.(names{1});
% correctWBMfields(model, 'female'), COBRA Toolbox commit 67c790d
if isfield(model, 'gender')
    model.sex = model.gender;
    model = rmfield(model, 'gender');
else
    model.sex = 'female';
end
if isfield(model, 'rxnGeneMat')
    model = rmfield(model, 'rxnGeneMat');
end
end
