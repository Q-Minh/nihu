function clear_pool
% CLEAR_POOL Closes the current Matlab pool
% CLEAR_POOL Closes the current Matlab pool in a Matlab-version independent
%   way.
% 
%   Copyright: Laboratory of Acoustics and Studio Technologies, 2015.

v = version('-release');
year = sscanf(v, '%d');
if year >= 2015
    % Delete the pool
    delete(gcp('nocreate'));
else
    try
        matlabpool('close');    %#ok<DPOOL>
    catch
    end
end

end % of clear_pool
