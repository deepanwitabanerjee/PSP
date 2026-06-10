% run_indigoidine_substrate_screen.m
%
% Test which substrates can be growth-coupled using the indigoidine cMCS
% strategy from Supplementary Data 4, section 2.3.
%
% Required file:
%   iJN1462.mat containing iJN1462 and del_exchanges
%
% COBRA Toolbox must be initialized before running this script.

changeCobraSolver('gurobi', 'all');

load('iJN1462.mat', 'iJN1462', 'del_exchanges');

model = removeRxns(iJN1462, iJN1462.rxns(530));
model = changeObjective(model, model.rxns(529));
model = changeRxnBounds(model, putidaAromaticRxns(), 0, 'b');
model = changeRxnBounds(model, model.rxns(del_exchanges), 0, 'u');

FBA = optimizeCbModel(model);

KOs = {'ANHMK', 'D_LACt2pp', 'GLCDpp', 'H2CO3D', 'HCO3E', 'MDH', ...
    'MDH2', 'ME2', 'NACODA', 'ORNDC', 'PHAPC40', 'PPS', 'PROD2', 'TALA'};
geneKO = {'PP_0100', 'PP_0434', 'PP_0654', 'PP_0751', 'PP_0864', ...
    'PP_1251', 'PP_1444', 'PP_2082', 'PP_2168', 'PP_2925', 'PP_4735', ...
    'PP_4947', 'PP_5003', 'PP_5005', 'PP_5085', 'PP_5186'};

model_INDG = addIndigoidineRxn(model, 'gln__L_c', true);
FBA_indg = optimizeCbModel(model_INDG);
[minINDG, maxINDG] = fluxVariability(model_INDG, 100, 'max', {'EX_ind_e'});

model_R = changeRxnBounds(model_INDG, KOs, 0, 'b');
model_G = deleteModelGenes(model_INDG, geneKO);

Substrate = model_R.rxns(del_exchanges);
Bio_R = zeros(length(del_exchanges), 1);
min_R = zeros(length(del_exchanges), 1);
max_R = zeros(length(del_exchanges), 1);
Bio_G = zeros(length(del_exchanges), 1);
min_G = zeros(length(del_exchanges), 1);
max_G = zeros(length(del_exchanges), 1);

for i = 1:length(del_exchanges)
    tmp = changeRxnBounds(model_R, model_R.rxns(del_exchanges(i)), -10, 'l');
    FBA_R = optimizeCbModel(tmp);
    Bio_R(i) = FBA_R.f;
    if Bio_R(i) ~= 0
        [min_R(i), max_R(i)] = fluxVariability(tmp, 100, 'max', {'EX_ind_e'});
    end
end

for i = 1:length(del_exchanges)
    tmp = changeRxnBounds(model_G, model_G.rxns(del_exchanges(i)), -10, 'l');
    FBA_G = optimizeCbModel(tmp);
    Bio_G(i) = FBA_G.f;
    if Bio_G(i) ~= 0
        [min_G(i), max_G(i)] = fluxVariability(tmp, 100, 'max', {'EX_ind_e'});
    end
end

SubstrateResults = table(Substrate(:), Bio_R, min_R, max_R, Bio_G, min_G, max_G, ...
    'VariableNames', {'Substrate', 'BiomassReactionKO', 'MinIndigoidineReactionKO', ...
    'MaxIndigoidineReactionKO', 'BiomassGeneKO', 'MinIndigoidineGeneKO', ...
    'MaxIndigoidineGeneKO'});

save('indigoidine_substrate_screen_results.mat', ...
    'FBA', 'FBA_indg', 'minINDG', 'maxINDG', 'SubstrateResults', ...
    'Bio_R', 'min_R', 'max_R', 'Bio_G', 'min_G', 'max_G');
