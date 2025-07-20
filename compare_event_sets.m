function isEqual = compare_event_sets(OmsA, OmsB, mk)
    % Number of events in each set
    nA = size(OmsA, 1) / mk;
    nB = size(OmsB, 1) / mk;

    % Sanity check
    if mod(nA,1) ~= 0 || mod(nB,1) ~= 0
        error('Each matrix must be a multiple of mk rows.');
    end
    if size(OmsA, 2) ~= size(OmsB, 2)
        error('Event matrices must have the same number of columns.');
    end

    % Extract blocks
    blocksA = cell(nA, 1);
    for i = 1:nA
        blocksA{i} = OmsA((i-1)*mk+1:i*mk, :);
    end

    blocksB = cell(nB, 1);
    for i = 1:nB
        blocksB{i} = OmsB((i-1)*mk+1:i*mk, :);
    end

    % Check A ⊆ B
    unmatched = true(1, nA);
    for i = 1:nA
        for j = 1:nB
            if isequal(blocksA{i}, blocksB{j})
                unmatched(i) = false;
                break;
            end
        end
    end

    % If all events in A are found in B, and same size, sets are equal
    isEqual = ~(any(unmatched) || nA ~= nB);
end