function sizeofMCS = calculateSizeOfMCS(mcs)
%CALCULATESIZEOFMCS Return the number of reaction cuts in each MCS entry.
%
% sizeofMCS = calculateSizeOfMCS(mcs)
%
% mcs is the CellNetAnalyzer output cell array described in Supplementary
% Data 4. The size is calculated from the first cut set for each metabolite,
% matching the original workflow.

sizeofMCS = zeros(length(mcs), 1);

for j = 1:length(mcs)
    modes = mcs{1, j};
    if ~isempty(modes)
        sizeofMCS(j, 1) = sum(modes(1, :) ~= 0);
    end
end
end
