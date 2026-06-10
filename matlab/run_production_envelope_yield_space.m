% run_production_envelope_yield_space.m
%
% Calculate indigoidine production envelopes and normalized yield space from
% Supplementary Data 4, section 2.4.
%
% Required files:
%   ../CellNetAnalyzer/Scripts/results.mat
%   ./CellNetAnalyzer/Scripts/iJN1462D.mat
%   ../other_info.mat, containing iJN1462
%
% COBRA Toolbox must be initialized before running this script.

changeCobraSolver('gurobi', 'all');

load('../CellNetAnalyzer/Scripts/results.mat');
load('./CellNetAnalyzer/Scripts/iJN1462D.mat', 'del_exchanges', 'reac_off');
load('../other_info.mat', 'iJN1462');

model = setupPutidaKT2440Model(iJN1462, del_exchanges, reac_off);
FBA = optimizeCbModel(model);
model_INDG = addIndigoidineRxn(model, 'gln__L_c', true);

geneKO = {'PP_0100', 'PP_0434', 'PP_0654', 'PP_0751', 'PP_0864', ...
    'PP_1251', 'PP_1444', 'PP_2082', 'PP_2168', 'PP_2925', 'PP_4735', ...
    'PP_4947', 'PP_5003', 'PP_5005', 'PP_5085', 'PP_5186'};
geneKO_tar = {'PP_0434', 'PP_0654', 'PP_0864', 'PP_1251', 'PP_1444', ...
    'PP_2082', 'PP_2168', 'PP_2925', 'PP_4735', 'PP_4947', 'PP_5003', ...
    'PP_5005', 'PP_5085', 'PP_5186'};

[biomassValuesWT, targetValuesWT, lineHandleWT] = ...
    productionEnvelope(model_INDG, {}, 'k', 'EX_ind_e', 'BiomassKT2440_Core2');
[biomassValues16, targetValues16, lineHandle16] = ...
    productionEnvelope(model_INDG, geneKO, 'r', 'EX_ind_e', 'BiomassKT2440_Core2', true);
[biomassValues14, targetValues14, lineHandle14] = ...
    productionEnvelope(model_INDG, geneKO_tar, 'b', 'EX_ind_e', 'BiomassKT2440_Core2', true);

nPoints = min([20, size(biomassValuesWT, 1), size(biomassValues16, 1), size(biomassValues14, 1)]);
Yield = zeros(nPoints, 9);

for i = 1:nPoints
    Yield(i, 1) = biomassValuesWT(i, 1) / max(biomassValuesWT);
    Yield(i, 2) = targetValuesWT(i, 1) / max(targetValuesWT(:, 2));
    Yield(i, 3) = targetValuesWT(i, 2) / max(targetValuesWT(:, 2));
    Yield(i, 4) = biomassValues16(i, 1) / max(biomassValuesWT);
    Yield(i, 5) = targetValues16(i, 1) / max(targetValuesWT(:, 2));
    Yield(i, 6) = targetValues16(i, 2) / max(targetValuesWT(:, 2));
    Yield(i, 7) = biomassValues14(i, 1) / max(biomassValuesWT);
    Yield(i, 8) = targetValues14(i, 1) / max(targetValuesWT(:, 2));
    Yield(i, 9) = targetValues14(i, 2) / max(targetValuesWT(:, 2));
end

save('production_envelope_yield_space_results.mat', ...
    'FBA', 'biomassValuesWT', 'targetValuesWT', 'lineHandleWT', ...
    'biomassValues16', 'targetValues16', 'lineHandle16', ...
    'biomassValues14', 'targetValues14', 'lineHandle14', 'Yield');
