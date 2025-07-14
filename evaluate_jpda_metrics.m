function tracking_metrics = evaluate_jpda_metrics(xhv, xv, T, nt, association_threshold, scan_start, scan_eval, scan_final)

% 🔒 Перетворення індексів у цілі числа
scan_start = round(scan_start);
scan_eval = round(scan_eval);
scan_final = round(scan_final);

% Визначення вкладеної структури
nested = iscell(xhv{1});  % Якщо xhv{1} — cell, то це вкладена структура: xhv{t}{i}

get_vec = @(X, t, i) nested .* X{t}{i} + ~nested .* X{t, i};


nCases = 0;
track_indices = [];

% --- nCases ---
for i = 1:nt
    dx = norm(xhv{scan_start, i}([1,3]) - xv{scan_start, i}([1,3]));
    if dx < association_threshold
        nCases = nCases + 1;
        track_indices(end+1) = i;
    end
end

% --- nOK ---
nOK = 0;
for i = track_indices
    dx = norm(xhv{scan_eval, i}([1,3]) - xv{scan_eval, i}([1,3]));
    if dx < association_threshold
        nOK = nOK + 1;
    end
end

% --- nSwitched ---
nSwitched = 0;
for i = track_indices
    dx_to_true = norm(xhv{scan_eval, i}([1,3]) - xv{scan_eval, i}([1,3]));
    switched = false;
    for j = 1:nt
        if j ~= i
            dx_other = norm(xhv{scan_eval, i}([1,3]) - xv{scan_eval, j}([1,3]));
            if dx_other < association_threshold && dx_to_true > association_threshold
                switched = true;
                break;
            end
        end
    end
    if switched
        nSwitched = nSwitched + 1;
    end
end

% --- nLost ---
nLost = 0;
for i = track_indices
    matched = false;
    for j = 1:nt
        dx = norm(xhv{scan_eval, i}([1,3]) - xv{scan_eval, j}([1,3]));
        if dx < association_threshold
            matched = true;
            break;
        end
    end
    if ~matched
        nLost = nLost + 1;
    end
end

% --- nMerged ---
target_counts = zeros(1, nt);
for i = 1:nt
    for j = 1:nt
        % if i ~= j 
            dx = norm(xhv{scan_eval, j}([1,3]) - xv{scan_eval, i}([1,3]));
            if dx < association_threshold
                target_counts(i) = target_counts(i) + 1;
            end
        % end
    end
end
nMerged = sum(target_counts > 1);

% --- nResult ---
nResult = 0;
for i = 1:nt
    dx = norm(xhv{scan_final, i}([1,3]) - xv{scan_final, i}([1,3]));
    if dx < association_threshold
        nResult = nResult + 1;
    end
end

% --- CFT ---
CFT = 0;
for i = 1:nt
    ever_matched = false;
    for k = 1:T
        for j = 1:nt
            dx = norm(xhv{k, i}([1,3]) - xv{k, j}([1,3]));
            if dx < association_threshold
                ever_matched = true;
                break;
            end
        end
        if ever_matched
            break;
        end
    end
    if ~ever_matched
        CFT = CFT + 1;
    end
end

% --- Return structure ---
tracking_metrics = struct('nCases', nCases, ...
                          'nOK', nOK, ...
                          'nSwitched', nSwitched, ...
                          'nMerged', nMerged, ...
                          'nLost', nLost, ...
                          'nResult', nResult, ...
                          'CFT', CFT);
end