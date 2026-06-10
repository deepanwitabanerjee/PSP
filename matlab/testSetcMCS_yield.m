function [Bio, Bio_gene, KO_rxn_names, KO_geneList] = testSetcMCS_yield(met_num, sizeofMCS, mcs, model)
%TESTSETCMCS_YIELD Test reaction and gene KO phenotypes for one cMCS target.
%
% [Bio, Bio_gene, KO_rxn_names, KO_geneList] = testSetcMCS_yield(met_num, sizeofMCS, mcs, model)
%
% Requires COBRA Toolbox. Bio and Bio_gene columns are:
%   1 biomass before demand reaction
%   2 biomass after demand reaction
%   3 minimum target demand flux at 100% biomass optimum
%   4 maximum target demand flux at 100% biomass optimum

if sizeofMCS(met_num) == 0
    warning('testSetcMCS_yield:noCutset', 'This metabolite does not have a cut set.');
    Bio = cell(0, 0);
    Bio_gene = cell(0, 0);
    KO_rxn_names = cell(0, 0);
    KO_geneList = cell(0, 0);
    return;
end

Modes_met = mcs{1, met_num};
KO_rxn_names = cell(size(Modes_met, 2), size(Modes_met, 1));
KO_geneList = cell(size(Modes_met, 2), size(Modes_met, 1));
Bio = zeros(size(Modes_met, 1), 4);
Bio_gene = zeros(size(Modes_met, 1), 4);

for j = 1:size(Modes_met, 1)
    KO_rxn_index = find(Modes_met(j, :));
    KO_rxn_names(1:length(KO_rxn_index), j) = model.rxns(KO_rxn_index);
    genes = findGenesFromRxns(model, KO_rxn_names(1:length(KO_rxn_index), j));
    KO_geneList(1:length(genes), j) = genes;
end

for j = 1:size(Modes_met, 1)
    rxnsToDelete = KO_rxn_names(~cellfun('isempty', KO_rxn_names(:, j)), j);
    tmp_RXN = changeRxnBounds(model, rxnsToDelete, 0, 'b');
    FBA_prev = optimizeCbModel(tmp_RXN);
    tmp_RXN = addDemandReaction(tmp_RXN, tmp_RXN.mets(met_num));
    demandRxn = tmp_RXN.rxns{end};
    FBA_post = optimizeCbModel(tmp_RXN);

    Bio(j, 1) = FBA_prev.f;
    Bio(j, 2) = FBA_post.f;
    if FBA_post.f ~= 0
        [Bio(j, 3), Bio(j, 4)] = fluxVariability(tmp_RXN, 100, 'max', {demandRxn});
    end
end

for j = 1:size(Modes_met, 1)
    tmp_GENE = model;
    genesToDelete = KO_geneList(~cellfun('isempty', KO_geneList(:, j)), j);
    for i = 1:length(genesToDelete)
        tmp_GENE = deleteModelGenes(tmp_GENE, genesToDelete{i});
    end

    FBA_prev_gene = optimizeCbModel(tmp_GENE);
    tmp_GENE = addDemandReaction(tmp_GENE, tmp_GENE.mets(met_num));
    demandRxn = tmp_GENE.rxns{end};
    FBA_post_gene = optimizeCbModel(tmp_GENE);

    Bio_gene(j, 1) = FBA_prev_gene.f;
    Bio_gene(j, 2) = FBA_post_gene.f;
    if FBA_post_gene.f ~= 0
        [Bio_gene(j, 3), Bio_gene(j, 4)] = fluxVariability(tmp_GENE, 100, 'max', {demandRxn});
    end
end
end
