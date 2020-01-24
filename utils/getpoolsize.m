function n = getpoolsize()

% getpoolsize finds the number of workers in the current parallel pool
% and returns them in the output variable n
n = 0;
if ~isempty(ver('distcomp')) % Check if parallel toolbox installed
    p = gcp('nocreate');
    if ~isempty(p)
        n = p.NumWorkers;
    end
end
