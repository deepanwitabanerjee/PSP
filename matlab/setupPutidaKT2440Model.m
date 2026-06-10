function model = setupPutidaKT2440Model(iJN1462, del_exchanges, reac_off)
%SETUPPUTIDAKT2440MODEL Apply constraints used for the KT2440 simulations.
%
% model = setupPutidaKT2440Model(iJN1462, del_exchanges, reac_off)

model = removeRxns(iJN1462, iJN1462.rxns(530));
model = changeObjective(model, model.rxns(529));

if nargin >= 3 && ~isempty(reac_off)
    model = changeRxnBounds(model, model.rxns(reac_off), 0, 'b');
end

model = changeRxnBounds(model, putidaAromaticRxns(), 0, 'b');
model = changeRxnBounds(model, model.rxns(del_exchanges), 0, 'u');
model = changeRxnBounds(model, {'sink_PHAg', 'sink_pqqA_kt_c_'}, 0, 'b');
end
