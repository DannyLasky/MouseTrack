function allTrials = AverageMouse(allTrials)
    if ~isempty(allTrials)
        uniqueMouse = unique(allTrials(:, 2));
        tempVector = nan(size(uniqueMouse));
        for n = 1:length(uniqueMouse)
            tempVector(n) = mean(allTrials(allTrials(:,2) == uniqueMouse(n), 1));
        end
        allTrials = tempVector;
    end
end