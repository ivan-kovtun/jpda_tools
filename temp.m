% Example input
ind = {
    [1 2],        % Level 1
    [1 2 3],      % Level 2
    [1 3 4]         % Level 3
};

numLevels = numel(ind);
sizes = cellfun(@numel, ind);
totalPaths = prod(sizes);

% Preallocate with max possible size (optional)
validPaths = zeros(totalPaths, numLevels);  
count = 0;

% Create index ranges
indices = arrayfun(@(n) 1:n, sizes, 'UniformOutput', false);

% Loop over all index combinations with ndgrid-style indexing
[subs{1:numLevels}] = ndgrid(indices{:});
for i = 1:numel(subs{1})
    path = zeros(1, numLevels);
    for j = 1:numLevels
        path(j) = ind{j}(subs{j}(i));
    end

    % Validate: all elements unique except 1
    pathWithoutOnes = path(path ~= 1);
    if numel(pathWithoutOnes) == numel(unique(pathWithoutOnes))
        count = count + 1;
        validPaths(count, :) = path;
    end
end

% Trim unused preallocated rows
validPaths = validPaths(1:count, :);

% Display
disp(validPaths);