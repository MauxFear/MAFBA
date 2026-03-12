function [growthAtAcetate, thresholdIndex] = detect_acetate_threshold(biomassValues, acetateValues, respirationValues)
    dataTable = table(biomassValues(:), acetateValues(:), respirationValues(:), ...
        'VariableNames', {'Growth', 'Acetate', 'Respiration'});
    dataTable = sortrows(dataTable, 'Growth');

    slopeResp = diff(dataTable.Respiration) ./ diff(dataTable.Growth);
    finiteSlope = slopeResp(isfinite(slopeResp));

    thresholdIndex = [];
    if ~isempty(finiteSlope)
        maxSlope = max(finiteSlope);
        slopeThresholds = [0.01, 0.05, 0.20];
        acetateThresholds = [0.01, 0, 0];

        for idx = 1:numel(slopeThresholds)
            candidateRows = find(slopeResp <= maxSlope * slopeThresholds(idx));
            if isempty(candidateRows)
                continue;
            end

            acetateIdx = find(dataTable.Acetate(candidateRows) >= acetateThresholds(idx), 1, 'first');
            if ~isempty(acetateIdx)
                thresholdIndex = candidateRows(acetateIdx);
                break;
            end
        end
    end

    if isempty(thresholdIndex)
        acetateIdx = find(dataTable.Acetate >= 0.01, 1, 'first');
        if isempty(acetateIdx)
            thresholdIndex = 1;
        else
            thresholdIndex = acetateIdx;
        end
    end

    growthAtAcetate = dataTable.Growth(thresholdIndex);
end
