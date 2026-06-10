% run_cmcs_validation.m
%
% Validate cMCS reaction knockouts against associated gene knockouts for
% the glutamine target described in Supplementary Data 4, section 1.2.
%
% Required files:
%   ../CellNetAnalyzer/Scripts/results.mat
%   ./CellNetAnalyzer/Scripts/iJN1462D.mat
%   ../other_info.mat, containing iJN1462
%
% COBRA Toolbox must be initialized before running this script.

changeCobraSolver('gurobi', 'all');

load('../CellNetAnalyzer/Scripts/results.mat', 'mcs');
load('./CellNetAnalyzer/Scripts/iJN1462D.mat', 'del_exchanges', 'reac_off');
load('../other_info.mat', 'iJN1462');

model = setupPutidaKT2440Model(iJN1462, del_exchanges, reac_off);
FBA = optimizeCbModel(model);
WT_biomass = FBA.f;

met_num = 1120; % glutamine in the study workflow
sizeofMCS = calculateSizeOfMCS(mcs);

[Bio, Bio_gene, KO_rxn_names, KO_geneList] = testSetcMCS_yield(met_num, sizeofMCS, mcs, model);

save('cmcs_validation_results.mat', ...
    'Bio', 'Bio_gene', 'KO_rxn_names', 'KO_geneList', 'WT_biomass', 'met_num', 'sizeofMCS');
