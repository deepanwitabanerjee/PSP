% calculate_coupling_mcs_ppu_template.m
%
% Template for the CellNetAnalyzer constrained minimal cut set calculation
% described in Supplementary Data 4, section 1.1. Run from the
% CellNetAnalyzer main directory after `startcna`, using MATLAB R2017b,
% CNA 2018.2, and CPLEX 12.8.
%
% This file records the study-specific inputs and parameters. The final MCS
% calculation call depends on the CNA coupling-cMCS script distribution
% (`calculate_coupling_mcs_coli.m`) and should be adapted there.

load('Scripts/iJN1462D.mat', ...
    'atpm', 'del_exchanges', 'exchange_ind', 'glucose', 'mue', 'noKO', ...
    'o2in', 'iJN1462Bs', 'iJN1411Bs_reacNames', 'reac_off');

glucose_uptake_limit = 6.3;
min_yield = 0.8;
min_mue = 0.01;
min_atpm = 0.92;
screen_metabolites = 1120; % glutamine

% Expected output variables, saved by the adapted CNA cMCS workflow:
%   gluc_up
%   mcs
%   max_mue
%   max_prod
%
% save('Scripts/results.mat', 'gluc_up', 'mcs', 'max_mue', 'max_prod');
