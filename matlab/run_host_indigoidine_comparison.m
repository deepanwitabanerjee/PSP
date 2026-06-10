% run_host_indigoidine_comparison.m
%
% Compare biomass, glutamine demand, glucose uptake, ATP maintenance, and
% indigoidine export across the five host models from Supplementary Data 4,
% section 2.2.
%
% Required file:
%   iCW773.mat containing iCW773, iJN1462, iMM904, iRhto1108C, iML1515
%
% COBRA Toolbox must be initialized before running this script.

changeCobraSolver('gurobi', 'all');

load('iCW773.mat', 'iCW773', 'iJN1462', 'iMM904', 'iRhto1108C', 'iML1515');

hostNames = {'iCW773', 'iMM904', 'iRhto1108C', 'iJN1462', 'iML1515'};

iCW773 = changeObjective(iCW773, iCW773.rxns(971));
iMM904 = changeObjective(iMM904, iMM904.rxns(1521));
iRhto1108C = changeObjective(iRhto1108C, iRhto1108C.rxns(2115));
iJN1462 = changeObjective(iJN1462, iJN1462.rxns(529));
iML1515 = changeObjective(iML1515, iML1515.rxns(2669));

model_gln_773 = changeObjective(addDemandReaction(iCW773, 'gln-L_c'), 'DM_gln-L_c');
model_gln_904 = changeObjective(addDemandReaction(iMM904, 'gln__L_c'), 'DM_gln__L_c');
model_gln_1108 = changeObjective(addDemandReaction(iRhto1108C, 'gln__L_c'), 'DM_gln__L_c');
model_gln_1462 = changeObjective(addDemandReaction(iJN1462, 'gln__L_c'), 'DM_gln__L_c');
model_gln_1515 = changeObjective(addDemandReaction(iML1515, 'gln__L_c'), 'DM_gln__L_c');

FBA_773 = optimizeCbModel(iCW773);
FBA_904 = optimizeCbModel(iMM904);
FBA_1108 = optimizeCbModel(iRhto1108C);
FBA_1462 = optimizeCbModel(iJN1462);
FBA_1515 = optimizeCbModel(iML1515);

FBA_773_g = optimizeCbModel(model_gln_773);
FBA_904_g = optimizeCbModel(model_gln_904);
FBA_1108_g = optimizeCbModel(model_gln_1108);
FBA_1462_g = optimizeCbModel(model_gln_1462);
FBA_1515_g = optimizeCbModel(model_gln_1515);

Results_wt = zeros(4, 5);
Results_wt(1, :) = [FBA_773.f, FBA_904.f, FBA_1108.f, FBA_1462.f, FBA_1515.f];
Results_wt(2, :) = [FBA_773_g.f, FBA_904_g.f, FBA_1108_g.f, FBA_1462_g.f, FBA_1515_g.f];
Results_wt(3, :) = [ ...
    FBA_773.x(strcmp(iCW773.rxns, 'EX_glc(e)')), ...
    FBA_904.x(strcmp(iMM904.rxns, 'EX_glc__D_e')), ...
    FBA_1108.x(strcmp(iRhto1108C.rxns, 'EX_glc__D_e')), ...
    FBA_1462.x(strcmp(iJN1462.rxns, 'EX_glc__D_e')), ...
    FBA_1515.x(strcmp(iML1515.rxns, 'EX_glc__D_e'))];
Results_wt(4, :) = [ ...
    FBA_773.x(strcmp(iCW773.rxns, 'ATPM')), ...
    FBA_904.x(strcmp(iMM904.rxns, 'ATPM')), ...
    FBA_1108.x(strcmp(iRhto1108C.rxns, 'ATPM_c')), ...
    FBA_1462.x(strcmp(iJN1462.rxns, 'ATPM')), ...
    FBA_1515.x(strcmp(iML1515.rxns, 'ATPM'))];

model_INDG_773 = addIndigoidineRxn(iCW773, 'gln-L_c', false);
model_INDG_904 = addIndigoidineRxn(iMM904, 'gln__L_c', false);
model_INDG_1108 = addIndigoidineRxn(iRhto1108C, 'gln__L_c', false);
model_INDG_1462 = addIndigoidineRxn(iJN1462, 'gln__L_c', true);
model_INDG_1515 = addIndigoidineRxn(iML1515, 'gln__L_c', true);

FBA_773_b = optimizeCbModel(model_INDG_773);
model_INDG_773 = changeObjective(model_INDG_773, 'EX_ind_e');
FBA_773_i = optimizeCbModel(model_INDG_773);
model_GLN_773 = changeObjective(addDemandReaction(model_INDG_773, 'gln-L_c'), 'DM_gln-L_c');
FBA_773_g = optimizeCbModel(model_GLN_773);

FBA_904_b = optimizeCbModel(model_INDG_904);
model_INDG_904 = changeObjective(model_INDG_904, 'EX_ind_e');
FBA_904_i = optimizeCbModel(model_INDG_904);
model_GLN_904 = changeObjective(addDemandReaction(model_INDG_904, 'gln__L_c'), 'DM_gln__L_c');
FBA_904_g = optimizeCbModel(model_GLN_904);

FBA_1108_b = optimizeCbModel(model_INDG_1108);
model_INDG_1108 = changeObjective(model_INDG_1108, 'EX_ind_e');
FBA_1108_i = optimizeCbModel(model_INDG_1108);
model_GLN_1108 = changeObjective(addDemandReaction(model_INDG_1108, 'gln__L_c'), 'DM_gln__L_c');
FBA_1108_g = optimizeCbModel(model_GLN_1108);

FBA_1462_b = optimizeCbModel(model_INDG_1462);
model_INDG_1462 = changeObjective(model_INDG_1462, 'EX_ind_e');
FBA_1462_i = optimizeCbModel(model_INDG_1462);
model_GLN_1462 = changeObjective(addDemandReaction(model_INDG_1462, 'gln__L_c'), 'DM_gln__L_c');
FBA_1462_g = optimizeCbModel(model_GLN_1462);

FBA_1515_b = optimizeCbModel(model_INDG_1515);
model_INDG_1515 = changeObjective(model_INDG_1515, 'EX_ind_e');
FBA_1515_i = optimizeCbModel(model_INDG_1515);
model_GLN_1515 = changeObjective(addDemandReaction(model_INDG_1515, 'gln__L_c'), 'DM_gln__L_c');
FBA_1515_g = optimizeCbModel(model_GLN_1515);

Results_indg = zeros(5, 5);
Results_indg(1, :) = [FBA_773_b.f, FBA_904_b.f, FBA_1108_b.f, FBA_1462_b.f, FBA_1515_b.f];
Results_indg(2, :) = [FBA_773_g.f, FBA_904_g.f, FBA_1108_g.f, FBA_1462_g.f, FBA_1515_g.f];
Results_indg(3, :) = Results_wt(3, :);
Results_indg(4, :) = Results_wt(4, :);
Results_indg(5, :) = [FBA_773_i.f, FBA_904_i.f, FBA_1108_i.f, FBA_1462_i.f, FBA_1515_i.f];

rowNamesWT = {'biomass', 'glutamineDemand', 'glucoseUptake', 'atpm'};
rowNamesINDG = {'biomass', 'glutamineDemand', 'glucoseUptake', 'atpm', 'indigoidineExport'};
Results_wt_table = array2table(Results_wt, 'VariableNames', hostNames, 'RowNames', rowNamesWT);
Results_indg_table = array2table(Results_indg, 'VariableNames', hostNames, 'RowNames', rowNamesINDG);

save('host_indigoidine_comparison_results.mat', ...
    'Results_wt', 'Results_indg', 'Results_wt_table', 'Results_indg_table');
