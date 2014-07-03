%INSTALL Install NiHu under Matlab
%   INSTALL includes the necessary directories of toolbox NiHu to the
%   Matlab path

% get the directory of this m file
name = mfilename('fullpath');
% extract the path only
path = fileparts(name);

% these directories will be added
dirs = {'analytic', 'compatibility', 'meshing', 'nihudemos'};
for i = 1 : length(dirs)
    addpath(fullfile(path, dirs{i}));
end
% save the path if you have admin rights
savepath();

% run the demos
rehash toolbox
demo toolbox NiHu