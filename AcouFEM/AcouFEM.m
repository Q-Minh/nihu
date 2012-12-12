function AcouFEM
%AcouFEM Acoustic Finite Element software GUI

%% Figure initialization
fh = figure('MenuBar', 'none',...
    'Name', 'AcouFEM',...
    'NumberTitle', 'off',...
    'Resize', 'on', 'ResizeFcn', @ResizeCallback, ...
    'Toolbar', 'figure',...
    'Visible', 'on',...
    'Units', 'pixels');

%% Axes initialization
axes('Parent', fh, ...
    'Tag', 'Axes', ...
    'Units', 'pixels');
xlabel('x');
ylabel('y');
zlabel('z');
axis equal;
colorbar;
grid;

%% StatusBar Initialization
uicontrol(fh, 'Style', 'Edit', ...
    'Tag', 'StatusBar', ...
    'HorizontalAlignment', 'left', ...
    'Min', 1, 'Max', 3);

%% Menubar initialization
% Model Menu
mhMod = uimenu(fh, 'Label', 'Model');
ehMod1 = uimenu(mhMod,'Label','New Model');
uimenu(ehMod1,'Label','3D Cartesian',...
    'Callback', {@NewModel, '3DCart'});
uimenu(ehMod1,'Label','2D Cartesian',...
    'Callback', {@NewModel, '2DCart'});
sehMod13 = uimenu(ehMod1,'Label','2.5D Cartesian');
uimenu(sehMod13,'Label','Inifinite Geometry',...
    'Callback', {@NewModel, '2.5DCartInf'});
uimenu(sehMod13,'Label','Finite Geometry',...
    'Callback', {@NewModel, '2.5DCartFin'});
uimenu(ehMod1,'Label','2D Cylindrical',...
    'Callback', {@NewModel, '2DCyl'});
sehMod15 = uimenu(ehMod1,'Label','2.5D Cylindrical');
uimenu(sehMod15,'Label','Full geometry',...
    'Callback', {@NewModel, '2.5DCylFull'});
uimenu(sehMod15,'Label','Finite Geometry',...
    'Callback', {@NewModel, '2.5DCylFin'});
uimenu(mhMod, 'Label', 'Load Model Data', 'Callback', @LoadModelData);
uimenu(mhMod, 'Label', 'Save Model Data', 'Tag', 'mSaveModel', ...
    'Callback', @SaveModelData, 'Separator','on');
uimenu(mhMod, 'Label', 'Export Data to Worksapce', 'Tag', 'mExportData', ...
    'Callback', @ExportData, 'Separator','on');
uimenu(mhMod, 'Label', 'Exit', 'Tag', 'mExit', ...
    'Callback', @Exit, 'Separator','on');
% Mesh Menu
mhMes = uimenu(fh,'Label','Mesh', 'Tag', 'mMesh');
uimenu(mhMes,'Label','Load Mesh','CallBack', {@LoadMesh});
uimenu(mhMes,'Label','Import Mesh','Callback',{@ImportMesh});
uimenu(mhMes,'Label','Save Mesh','Tag','mSaveMesh',...
    'Callback',{@SaveMesh}, 'Separator','on');
% Matrices Menu
mhMat = uimenu(fh,'Label','System Matrices', 'Tag', 'mSystemMatrices');
uimenu(mhMat,'Label', 'Compute Mass and Stiffness Matrices', ...
    'Callback', {@ComputeMK});
% Constraints Menu
uimenu(fh, 'Label', 'Constraints', 'Tag', 'mConstraints');
% Modal Analysis menu
mhMode = uimenu(fh, 'Label', 'Modal Analysis', 'Tag', 'mModalAnalysis');
uimenu(mhMode, 'Label', 'Compute Modes', 'Callback', {@ComputeModes});
uimenu(mhMode, 'Label', 'Display Mode', 'Tag', 'mDisplayModes', ...
    'Callback', {@DisplayMode}, 'Separator','on');
% Excitation menu
mhExc = uimenu(fh, 'Label', 'Excitation', 'Tag', 'mExcitation');
uimenu(mhExc, 'Label', 'Frequency Range', 'Tag', 'mFrequencyRange', ...
    'Callback', {@FrequencyRange});
uimenu(mhExc, 'Label', 'Surface Normal velocity Excitation', ...
    'Tag', 'mNormalVelocityExcitation', ...
    'Callback', {@NormalVelocityExcitation}, 'Separator','on');
uimenu(mhExc, 'Label', 'Display Excitation', 'Tag', 'mDisplayExcitation', ...
    'Callback', {@DisplayExcitation});
uimenu(mhExc, 'Label', 'Transform Excitation to Wavenumber Domain', ...
    'Tag', 'mTransformExcitation', 'Callback', {@TransformExcitation});
uimenu(mhExc, 'Label', 'Compute Excitation Matrix', ...
    'Tag', 'mComputeExcitationMatrix', 'Callback', {@ComputeExcitationMatrix});
% Solution menu
mhSol = uimenu(fh, 'Label', 'Solution', 'Tag', 'mSolution');
uimenu(mhSol, 'Label', 'Direct Solution', 'Tag', 'mDirectSolution', ...
    'Callback', {@DirectSolution});
uimenu(mhSol, 'Label', 'Modal Solution', 'Tag', 'mModalSolution', ...
    'Callback', {@ModalSolution});
uimenu(mhSol, 'Label', 'Transform Solution', ...
    'Tag', 'mTransformSolution', 'Callback', {@TransformSolution}, 'Separator','on');
uimenu(mhSol, 'Label', 'Display Solution', 'Tag', 'mDisplaySolution', ...
    'Callback', {@DisplaySolution}, 'Separator','on');
% Help Menu
mhHel = uimenu(fh, 'Label', 'Help');
uimenu(mhHel, 'Label', 'AcouFEM documentation', ...
    'Callback', {@Documentation});
uimenu(mhHel, 'Label', 'About', 'Callback', {@About}, 'Separator','on');

%% data and LogFile initialisation
data.LogString = {};
AcouFEMLog(sprintf('%%%% %s - AcouFEM started', datestr(now)));
data = SInitData(data);
AcouFEMLog('data = SInitData;');
guidata(fh, data);

UpdateMenu(fh);
UpdateStatusBar(fh, 'AcouFEM Started');


%% Window Resize Callback
function ResizeCallback(hObject, eventdata)

% define borders
sbBorder = 10;
axBorder = 50;
sbHeight = 80;
% get figure position and size
fig = gcbo;
figpos = get(fig,'Position');
FigWidth = figpos(3);
FigHeight = figpos(4);
% set StatusBar position and size
h = guihandles(hObject);
sb = h.StatusBar;
sbPos = [sbBorder sbBorder FigWidth-2*sbBorder, sbHeight];
set(sb, 'Position', sbPos);
% set Figure Axes position and size
ax = h.Axes;
axPos = [sbBorder sbBorder+sbHeight+axBorder FigWidth-2*sbBorder FigHeight-2*axBorder-sbHeight-sbBorder];
set(ax, 'Position', axPos);


%% Menu updating
function UpdateMenu(hObject)

h = guihandles(hObject);
set(h.mSaveModel, 'Enable', 'off');
set(h.mMesh, 'Enable', 'off');
set(h.mSaveMesh, 'Enable', 'off');
set(h.mSystemMatrices, 'Enable', 'off');
set(h.mConstraints, 'Enable', 'off');
set(h.mModalAnalysis, 'Enable', 'off');
set(h.mDisplayModes, 'Enable', 'off');
set(h.mExcitation, 'Enable', 'off');
set(h.mDisplayExcitation, 'Enable', 'off');
set(h.mTransformExcitation, 'Enable', 'off');
set(h.mComputeExcitationMatrix, 'Enable', 'off');
set(h.mSolution, 'Enable', 'off');
set(h.mDirectSolution, 'Enable', 'off');
set(h.mModalSolution, 'Enable', 'off');
set(h.mTransformSolution, 'Enable', 'off');
set(h.mDisplaySolution, 'Enable', 'off');

data = guidata(hObject);
if ~isempty(data.ModelType)
    set(h.mSaveModel, 'Enable', 'on');
    set(h.mMesh, 'Enable', 'on');
end

if exists(data.Model, 'Domain')
    set(h.mSaveMesh, 'Enable', 'on');
    set(h.mSystemMatrices, 'Enable', 'on');
    set(h.mExcitation, 'Enable', 'on');
end

if exists(data.Matrices, 'M')
    set(h.mModalAnalysis, 'Enable', 'on');
end

if exists(data.Modes, 'Phi')
    set(h.mDisplayModes, 'Enable', 'on');
end

if exists(data.Excitation, 'NormalVelocity')
    set(h.mDisplayExcitation, 'Enable', 'on');
    set(h.mComputeExcitationMatrix, 'Enable', 'on');
    if strcmp(data.ModelType, '2.5DCartInf') || strcmp(data.ModelType, '2.5DCylFull')
        set(h.mTransformExcitation, 'Enable', 'on');
    end
end

if exists(data.Excitation, 'NormalVelocity')
    NV = data.Excitation.NormalVelocity(1);
    set(h.mSolution, 'Enable', 'on');
    if ((isfield(NV, 'vhat') && ...
            (strcmp(data.ModelType, '2DCart') || strcmp(data.ModelType, '3DCart') || strcmp(data.ModelType, '2DCyl'))) || ...
            (isfield(NV, 'vtil') && ...
            (strcmp(data.ModelType, '2.5DCartInf') || strcmp(data.ModelType, '2.5DCartFin') || strcmp(data.ModelType, '2.5DCylFin') || strcmp(data.ModelType, '2.5DCylFull')))) && ...
            (exists(NV, 'A'))
        set(h.mDirectSolution, 'Enable', 'on');
        if exists(data.Modes, 'Phi')
            set(h.mModalSolution, 'Enable', 'on');
        end
    end
end

if exists(data.Response, 'ptil')
    set(h.mTransformSolution, 'Enable', 'on');
end

if exists(data.Response, 'phat')
    set(h.mDisplaySolution, 'Enable', 'on');
end


%% StatusBar Updating
function UpdateStatusBar(hObject, String)

h = guihandles(hObject);
sb = h.StatusBar;
data = guidata(hObject);
data.LogString{end+1} = sprintf('%s\t%s', datestr(now), String);
nString = length(data.LogString);
str = '';
for iString = nString : -1 : 1
    str = [str sprintf('%3d.  %s\n', iString, data.LogString{iString})];
end
set(sb, 'String', str(1:end-1));
guidata(hObject ,data);


%% New Model Menu CallBack
function NewModel(hObject, eventdata, type)

data = guidata(hObject);

switch type
    case '2.5DCartFin'
        answers = paramdlg('Longitudinal (z) Axis', ...
            {'Zmax', '1', 'isreal(str2double(x)) & isfinite(str2double(x)) & str2double(x) > 0', 'Zmax should be a positive real value.';
            'nZ', '10', 'floor(str2double(x)) == str2double(x) & floor(str2double(x)) > 0', 'nZ should be a positive integer.'});
        if ~isempty(answers)
            Zmax = str2double(answers{1});
            Nz = str2double(answers{2});
            data = SNewModel(data, type, Zmax, Nz);
            AcouFEMLog(sprintf('data = SNewModel(data, ''%s'', %f, %d);', type, Zmax, Nz));
        end
    case '2.5DCartInf'
        answers = paramdlg('Longitudinal (z) Axis', ...
            {'Zmin', '0', 'isreal(str2double(x)) & isfinite(str2double(x))', 'Zmin should be a finite real value.';
            'Zmax', '1', 'isreal(str2double(x)) & isfinite(str2double(x))', 'Zmax should be a finite real value.';
            'nZ', '10', 'floor(str2double(x)) == str2double(x) & floor(str2double(x)) > 0', 'nZ should be a positive integer.'});
        if ~isempty(answers)
            Zmin = str2double(answers{1});
            Zmax = str2double(answers{2});
            Nz = str2double(answers{3});
            data = SNewModel(data, type, Zmin, Zmax, Nz);
            AcouFEMLog(sprintf('data = SNewModel(data, ''%s'', %f, %f, %d);', type, Zmin, Zmax, Nz));
        end
    case '2.5DCylFin'
        answers = paramdlg('Angular (phi) Axis', ...
            {'PhiMax', 'pi', 'isreal(str2double(x)) & isfinite(str2double(x)) & str2double(x) > 0 & str2double(x) < 2*pi', 'Zmax should be a real value between 0 and 2pi.';
            'nPhi', '10', 'floor(str2double(x)) == str2double(x) & floor(str2double(x)) > 0', 'nPhi should be a positive integer.'});
        if ~isempty(answers)
            Phimax = str2double(answers{1});
            nPhi = str2double(answers{2});
            data = SNewModel(data, type, Phimax, nPhi);
            AcouFEMLog(sprintf('data = SNewModel(data, ''%s'', %f, %d);', type, Phimax, nPhi));
        end
    case '2.5DCylFull'
        answers = paramdlg('Angular (phi) Axis', ...
            {'nPhi', '10', 'floor(str2double(x)) == str2double(x) & floor(str2double(x)) > 0', 'nPhi should be a positive integer.'});
        if ~isempty(answers)
            nPhi = str2double(answers{1});
            data = SNewModel(data, type, nPhi);
            AcouFEMLog(sprintf('data = SNewModel(data, ''%s'', %d);', type, nPhi));
        end
    otherwise
        data = SNewModel(data, type);
        AcouFEMLog(sprintf('data = SNewModel(data, ''%s'');', type));
end
guidata(hObject, data);
cla;
UpdateMenu(hObject);
UpdateStatusBar(hObject, sprintf('New Model Cerated. Model type: %s', type));


%% Load Model Menu CallBack
function LoadModelData(hObject, eventdata)

[filename, pathname] = uigetfile({'*.mat', 'Matlab MAT file'}, 'Load Model Data');
if filename ~= 0
    load([pathname filename], 'data');
    guidata(hObject, data);
end
UpdateMenu(hObject);


%% Save Model Menu CallBack
function SaveModelData(hObject, eventdata)

data = guidata(hObject);
[filename, pathname] = uiputfile({'*.mat', 'Matlab MAT file'}, 'Save Model Data');
if filename ~= 0
    save([pathname filename], 'data');
end
UpdateMenu(hObject);


%% ExportModelData
function ExportData(hObject, eventdata)

data = guidata(hObject);

checkLabels = {'Model named:' ...
    'Matrices named:' ...
    'Modes named:' ...
    'Excitation named:' ...
    'Response named:'};
varNames = {'Model','Matrices', 'Modes', 'Excitation', 'Response'};
items = {data.Model, data.Matrices, data.Modes, data.Excitation, data.Response};
export2wsdlg(checkLabels,varNames,items,'Export data to Workspace');


%% Exit Function
function Exit(hObject, eventdata)

delete(gcf);


%% Load Mesh Menu CallBack
function LoadMesh(hObject, eventdata)

data = guidata(hObject);
[filename, pathname] = uigetfile({'*.mat', 'Matlab MAT file'}, 'Load Mesh');
if filename ~= 0
    data = SLoadMesh(data, [pathname filename]);
    AcouFEMLog(sprintf('data = SLoadMesh(data, ''%s'');', [pathname filename]));
    guidata(hObject, data);
    PlotModel(data.Model);
    UpdateMenu(hObject);
    [nNodes, nElements] = mesh_statistics(data.Model.Domain);
    UpdateStatusBar(hObject, sprintf('New Mesh Loaded. Number of nodes: %d. Number of elements: %d', nNodes, nElements));
end


%% Import Mesh Menu CallBack
function ImportMesh(hObject, eventdata)

data = guidata(hObject);
[filename, pathname] = uigetfile({'*.bdf', 'Nastran Bulk format'}, 'Import Mesh');
if filename ~= 0
    data = SImportMesh(data, [pathname filename]);
    AcouFEMLog(sprintf('data = SImportMesh(data, ''%s'');', [pathname filename]));
    guidata(hObject, data);
    PlotModel(data.Model);
    UpdateMenu(hObject);
    [nNodes, nElements] = mesh_statistics(data.Model.Domain);
    UpdateStatusBar(hObject, sprintf('New Mesh Imported. Number of nodes: %d. Number of elements: %d', nNodes, nElements));
end


%% Plot Model
function PlotModel(Model)

cla;
plot_mesh(Model.PlotDomain);


%% Save Mesh Menu CallBack
function SaveMesh(hObject, eventdata)

data = guidata(hObject);
Domain = data.Model.Domain;
[filename, pathname] = uiputfile({'*.mat', 'Matlab MAT file'}, 'Save mesh');
if filename ~= 0
    save([pathname filename], 'Domain');
    UpdateStatusBar(hObject, sprintf('Mesh Saved to file %s', filename));
end


%% Compute MK Menu CallBack
function ComputeMK(hObject, eventdata)

data = guidata(hObject);
tic;
data = SComputeMK(data);
AcouFEMLog('data = SComputeMK(data);');
tElapsed = toc;
guidata(hObject, data);
UpdateMenu(hObject);
UpdateStatusBar(hObject, sprintf('System matrices computed. Time elapsed: %.2f seconds', tElapsed));


%% Compute Modes Menu CallBack
function ComputeModes(hObject, eventdata)

data = guidata(hObject);

Values = paramdlg('Modal Analysis Options', {'Number of modes', '100', 'str2double(x) == floor(str2double(x)) & str2double(x) > 0', 'Number of modes must be a positive integer.';
    'Maximal number of iterations', '300', 'str2double(x) == floor(str2double(x)) & str2double(x) > 0', 'Maximal number of iterations bust be a positive integer'});
if ~isempty(Values)
    tic;
    data = SComputeModes(data, str2double(Values{1}), str2double(Values{2}));
    AcouFEMLog(sprintf('data = SComputeModes(data, %d, %d);', str2double(Values{1}), str2double(Values{2})));
    tElapsed = toc;
    guidata(hObject, data);
    UpdateMenu(hObject);
    UpdateStatusBar(hObject, sprintf('%d modes determined. Time elapsed: %.2f seconds', size(data.Modes.Om,1), tElapsed));
end


%% Display Modes Menu CallBack
function DisplayMode(hObject, eventdata)

data = guidata(hObject);
nModes = size(data.Modes.Phi,2);

ListString = cell(nModes,1);
for iModes = 1 : nModes
    ListString{iModes} = sprintf('Mode #%d at %.2f Hz', iModes, data.Modes.Om(iModes)/2/pi);
end
[iMode, ok] = listdlg('ListString', ListString, 'SelectionMode', 'single', 'Name', 'Mode Selection');
if ok
    cla;
    plot_mesh(data.Model.Domain, 'node', data.Modes.DOF,...
        data.Modes.Phi(:,iMode)./max(abs(data.Modes.Phi(:,iMode))));
    set(gca, 'Clim', [-1, 1]);
    title(sprintf('Mode #%d at frequency %.2f Hz', iMode, data.Modes.Om(iMode)/2/pi));
end


%% Frequency Range Menu CallBack
function FrequencyRange(hObject, eventdata)

data = guidata(hObject);
Values = paramdlg('Frequency Range', {'Scale Type (lin / log)', 'lin', '1', '';
    'Start Value [Hz]', '0', '1', '';
    'Stop Value [Hz]', '100', '1', '';
    'Number of frequency steps', '100', '1', ''});
if ~isempty(Values)
    type = Values{1};
    start = str2double(Values{2});
    stop = str2double(Values{3});
    steps = str2double(Values{4});
    data = SFrequencyRange(data, type, start, stop, steps);
    AcouFEMLog(sprintf('data = SFrequencyRange(data, ''%s'', %d, %d, %d);', type, start, stop, steps));
    guidata(hObject, data);
    UpdateMenu(hObject);
    UpdateStatusBar(hObject, sprintf('Frequency vector defined. Type: %s, %d steps between %.2f Hz and %.2f Hz ', type, steps, start, stop));
end


%% Normal Velocity Excitation Menu CallBack
function NormalVelocityExcitation(hObject, eventdata)

data = guidata(hObject);

switch data.ModelType
    case {'3DCart', '2DCart', '2DCyl'}
        boundary = data.Model.Boundary;
    case {'2.5DCartFin', '2.5DCylFin', '2.5DCartInf', '2.5DCylFull'}
        boundary = data.Model.PlotBoundary;
end

cla;
plot_mesh(boundary);

Values = paramdlg('Normal Velocity Excitation', {'Element Selector', '', '1', '';
    'Velocity function [m/s/Hz]', '', '1', ''});
if ~isempty(Values)
    selector = Values{1};
    expression = Values{2};
    data = SNormalVelocityExcitation(data, selector, expression);
    AcouFEMLog(sprintf('data = SNormalVelocityExcitation(data, ''%s'', ''%s'');', selector, expression));
    guidata(hObject, data);
    UpdateMenu(hObject);
    nV = length(data.Excitation.NormalVelocity);
    UpdateStatusBar(hObject, sprintf('Normal velocity excitation #%d added.', nV));
end


%% Display Excitation Menu CallBack
function DisplayExcitation(hObject, eventdata)

data = guidata(hObject);

NV = data.Excitation.NormalVelocity;
nV = length(NV);
ListString = cell(nV,1);
for iV = 1 : nV
    ListString{iV} = sprintf('Excitation #%d: %s', iV, NV(iV).ElemSelector);
end
[iV, ok] = listdlg('ListString', ListString, 'SelectionMode', 'single', 'Name', 'Excitation Selection');
if ok
    cla;
    plot_mesh(data.Model.PlotBoundary, 'node', data.Model.PlotBoundary.Nodes(:,1), NV(iV).vhat);
    mn = min(NV(iV).vhat);
    mx = max(NV(iV).vhat);
    if mn == mx
        mn = mn-1;
        mx = mx+1;
    end
    set(gca, 'clim', [mn mx]);
    colorbar;
    title('Excitation vector');
end


%% Transform Excitation Vector Menu CallBack
function TransformExcitation(hObject, eventdata)

data = guidata(hObject);
tic;
data = STransformExcitation(data);
AcouFEMLog('data = STransformExcitation(data);');
tElapsed = toc;
guidata(hObject, data);
UpdateMenu(hObject);
UpdateStatusBar(hObject, sprintf('Excitation transformed into wavenumber domain. Time elapsed: %.2f seconds', tElapsed));


%% Compute Excitation Matrix Menu CallBack
function ComputeExcitationMatrix(hObject, eventdata)

data = guidata(hObject);
tic;
data = SComputeExcitationMatrix(data);
AcouFEMLog('data = SComputeExcitationMatrix(data);');
tElapsed = toc;
guidata(hObject, data);
UpdateMenu(hObject);
UpdateStatusBar(hObject, sprintf('Excitation matrix computed. Time elapsed: %.2f seconds', tElapsed));


%% Direct Solution Menu CallBack
function DirectSolution(hObject, eventdata)

data = guidata(hObject);
tic;
data = SDirectSolution(data);
AcouFEMLog('data = SDirectSolution(data);');
tElapsed = toc;
guidata(hObject, data);
UpdateMenu(hObject);
UpdateStatusBar(hObject, sprintf('Direct solution performed. Time elapsed: %.2f seconds', tElapsed));


%% Modal Solution Menu CallBack
function ModalSolution(hObject, eventdata)

data = guidata(hObject);
tic;
data = SModalSolution(data);
AcouFEMLog('data = SModalSolution(data);');
tElapsed = toc;
guidata(hObject, data);
UpdateMenu(hObject);
UpdateStatusBar(hObject, sprintf('Modal solution performed. Time elapsed: %.2f seconds', tElapsed));


%% Transform Solution Menu CallBack
function TransformSolution(hObject, eventdata)

data = guidata(hObject);
tic;
data = STransformSolution(data);
AcouFEMLog('data = STransformSolution(data);');
tElapsed = toc;
guidata(hObject, data);
UpdateMenu(hObject);
UpdateStatusBar(hObject, sprintf('Solution Transformed into spatial domain. Time elapsed: %.2f seconds', tElapsed));


%% Display Solution Menu CallBack
function DisplaySolution(hObject, eventdata)

data = guidata(hObject);
Values = paramdlg('Solution Frequency [Hz]', {'Solution Frequency', '', '1', ''});
if ~isempty(Values)
    freq = str2double(Values{1});
    [thrash, ind] = min(abs(data.Excitation.FrequencyVector-freq));
    phat = data.Response.phat(:,ind);
    freq = data.Excitation.FrequencyVector(ind);
    cla;
    phat = imag(phat);
    DOF = data.Model.PlotDomain.Nodes(:,1);
    plot_mesh(data.Model.PlotDomain, 'node', DOF, phat);
    colorbar;
    set(gca, 'clim', [min(phat), max(phat)]);
    title(sprintf('Sound pressure at f = %.2f Hz', freq));
end


%% Documentation Menu CallBack
function Documentation(hObject, eventdata)

web http://vibac.hit.bme.hu/download/fiala/AcouFEM.pdf


%% About Menu CallBack
function About(hObject, eventdata)

fh = gcf;
rgb = imread('logoU_small.jpg');
s = [size(rgb,1), size(rgb,2)];
p = get(fh, 'Position');
p = p(1:2) + p(3:4)/2 - s/2;
h = figure('MenuBar', 'none', 'Name', 'About AcouFEM',...
    'NumberTitle', 'off', 'Position', [p s],...
    'Resize', 'off', 'Toolbar', 'none', 'Visible', 'on',...
    'Units', 'pixels');
axes('Parent', h, 'Position', [0 0 1 1]);
image(rgb);
axis off;