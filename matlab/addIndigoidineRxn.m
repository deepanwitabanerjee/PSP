function model = addIndigoidineRxn(model, glnMetId, includeFmn)
%ADDINDIGOIDINERXN Add indigoidine production and export reactions.
%
% model = addIndigoidineRxn(model)
% model = addIndigoidineRxn(model, glnMetId)
% model = addIndigoidineRxn(model, glnMetId, includeFmn)
%
% COBRA Toolbox is required. The default glutamine metabolite id follows
% BiGG-style models such as iJN1462, iMM904, and iML1515.

if nargin < 2 || isempty(glnMetId)
    glnMetId = 'gln__L_c';
end

if nargin < 3 || isempty(includeFmn)
    includeFmn = true;
end

if includeFmn
    indgFormula = sprintf(['2 %s + 2 atp_c + 2 coa_c + 2 fmn_c + 2.5 o2_c -> ' ...
        '2 pap_c + 2 fmnh2_c + 2 ppi_c + 2 amp_c + 2 ptth_c + h2o_c + 2 pi_c + ind_c'], glnMetId);
else
    indgFormula = sprintf(['2 %s + 2 atp_c + 2 coa_c + 2.5 o2_c -> ' ...
        '2 pap_c + 2 ppi_c + 2 amp_c + 2 ptth_c + h2o_c + 2 pi_c + ind_c'], glnMetId);
end

model = addReaction(model, 'INDG', indgFormula);
model = addReaction(model, 'APNPT_db', 'atp_c + ptth_c -> adp_c + h_c + pan4p_c');
model = addReaction(model, 'DF_INDG', 'ind_c -> ind_e');
model = addExchangeRxn(model, 'ind_e');
model = changeRxnBounds(model, 'EX_ind_e', 0, 'l');
end
