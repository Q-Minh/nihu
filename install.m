%% NiHu INSTALL
% This script installs the toolbox NiHu to the MATLAB system.
% Make sure that the archive 'nihu.zip' is in the current directory.
% The installation process consists of:
% 
% * Verifying MATLAB compatibility
% * Finding an install directory
% * Extracting the NiHu files to the install directory
% * Trying to find precompiled MEX files suitable for the MATLAB release
% * Compiling NiHu MEX files if previous step failed
% * Adding NiHu to the Matlab path
% * Trying to run precompiled MEX files
% * Running the tutorials in the web browser
% 
clear all;
close all;
clc;

%% Check if Matlab release is supported by NIHU
verstring = version('-release');
if str2double(verstring(1:end-1)) < 2006
    errordlg(sprintf('Matlab version r%s is not supported by NIHU. Update your Matlab to r2006a or newer.',...
        verstring), 'Unsupported Matlab version');
    return;
end

%% Get installation directory
installdir = uigetdir(fullfile(matlabroot, 'toolbox'),...
    ['Select installation directory for toolbox NIHU.'...
	'If you select ../dir, the toolbox will be installed in ../dir/nihu/']);
if installdir == 0
    disp('Installation cancelled by the User.');
    return;
end
installdir = fullfile(installdir, 'nihu');

%% Check whether install dir is empty or not
d = dir(installdir);
if isempty(d)
    mkdir(installdir);
    fprintf(1, 'Directory %s has been created\n', installdir);
else
    if length(d) > 2
        button = questdlg(...
            sprintf('The directory \n%s\nis not empty. Would you like to overwrite previous installation files?', installdir),...
            'Directory not empty','Overwrite','Abort','Overwrite');
        pause(0.5);
        switch button
            case 'Abort'
                disp('Installation cancelled by the User.');
                return;
        end
        % clear all functions from memory so that they can be deleted
        clear functions;
    end
end

%% Copy NIHU files to install dir
d = dir('nihu.zip');
if isempty(d)
    errordlg('Unable to find installation archive ''nihu.zip''');
    return;
else
    disp('Extracting NIHU files to install directory...');
    unzip('nihu.zip', installdir);
    % change the file attributes to write enable
    fileattrib(installdir,'+w','','s');
end

%% Copy appropriate precompiled MEX files or compile MEX scripts
disp('Locating precompiled MEX files...');
bemdir = fullfile(installdir, 'bem');
precfiles = fullfile(bemdir, 'mex_source',...
    ['precomp_r' verstring], ['*.' mexext]);
if ~isempty(dir(precfiles))
    copyfile(precfiles, bemdir);
else
    disp('Compiling MEX scripts...');
    try
        cd(fullfile(bemdir, 'mex_source'));
        make;
        movefile(sprintf('*.%s', mexext), '../');
        delete *.m;
    catch %#ok<CTCH>
        errordlg('Could not find precompiled MEX files for your system or compile C scripts.');
    end
end

try %#ok<TRYNC>
    rmdir(fullfile(bemdir, 'mex_source', 'precomp_r*'), 's');
end

%% Add NIHU to the Matlab path
disp('Adding NIHU to the Matlab path...');
path(fullfile(installdir, 'bem'), path);
path(fullfile(installdir, 'fem'), path);
path(fullfile(installdir, 'iem'), path);
path(fullfile(installdir, 'meshing'), path);
path(fullfile(installdir, 'misc'), path);
path(fullfile(installdir, 'nihudemos'), path);
savepath;

%% Test if MEX files are working correctly
cd(fullfile(installdir, 'bem'));
try
    if mex_test(1) ~= 2
        errordlg('Precompiled MEX functions are not running on Your system. Try to recompile them.');
    end
catch %#ok<CTCH>
    errordlg('Precompiled MEX functions are not running on Your system. Try to recompile them.');
end
cd(installdir);

%% Confirm installation
disp('Installation finished successfully.');

%% Run demos
demo toolbox nihu