function open_pool(numWorkers)
% OPEN_POOL Opens the current Matlab pool
% OPEN_POOL Opens the current Matlab pool in a Matlab-version independent
%   way.
% 
%   Copyright: Laboratory of Acoustics and Studio Technologies, 2015.

v = version('-release');
year = sscanf(v, '%d');
if year >= 2015
    poolobj = parpool(numWorkers); %#ok<NASGU>
else
    matlabpool(numWorkers); %#ok<DPOOL>
end

end % of open_pool
