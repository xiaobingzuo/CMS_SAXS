function varargout = agbh_findpeaks_12IDB(varargin)
% AGBH_FINDPEAKS_12IDB MATLAB code for agbh_findpeaks_12IDB.fig
%      AGBH_FINDPEAKS_12IDB, by itself, creates a new AGBH_FINDPEAKS_12IDB or raises the existing
%      singleton*.
%
%      H = AGBH_FINDPEAKS_12IDB returns the handle to a new AGBH_FINDPEAKS_12IDB or the handle to
%      the existing singleton*.
%
%      AGBH_FINDPEAKS_12IDB('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in AGBH_FINDPEAKS_12IDB.M with the given input arguments.
%
%      AGBH_FINDPEAKS_12IDB('Property','Value',...) creates a new AGBH_FINDPEAKS_12IDB or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before agbh_findpeaks_12IDB_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to agbh_findpeaks_12IDB_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help agbh_findpeaks_12IDB

% Last Modified by GUIDE v2.5 07-Aug-2020 12:03:58

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @agbh_findpeaks_12IDB_OpeningFcn, ...
                   'gui_OutputFcn',  @agbh_findpeaks_12IDB_OutputFcn, ...
                   'gui_LayoutFcn',  [] , ...
                   'gui_Callback',   []);
if nargin && ischar(varargin{1})
    gui_State.gui_Callback = str2func(varargin{1});
end

if nargout
    [varargout{1:nargout}] = gui_mainfcn(gui_State, varargin{:});
else
    gui_mainfcn(gui_State, varargin{:});
end
% End initialization code - DO NOT EDIT


% --- Executes just before agbh_findpeaks_12IDB is made visible.
function agbh_findpeaks_12IDB_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to agbh_findpeaks_12IDB (see VARARGIN)

% Choose default command line output for agbh_findpeaks_12IDB
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes agbh_findpeaks_12IDB wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = agbh_findpeaks_12IDB_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on button press in pb_load.
function pb_load_Callback(hObject, eventdata, handles)
% hObject    handle to pb_load (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
[filename, pathname, filterindex] = uigetfile( ...
   {'*.h5','HDF5 data (*.h5)'; ...
    '*.tif','TIFF image (*.tif)'; ...
    '*.bmp','Bitmap (*.bmp)'; ...
    '*.*',  'All Files (*.*)'}, ...
    'Pick a file');

%This code checks if the user pressed cancel on the dialog.

%[filename, pathname] = uigetfile('*.m', 'Pick a MATLAB code file');
if isequal(filename,0) || isequal(pathname,0)
   disp('User pressed cancel')
else
   disp(['User selected ', fullfile(pathname, filename)])
end

%a = imread(fullfile(pathname, filename));
[~, fn9, ext9] = fileparts(filename);

if strcmpi(ext9, '.tif')      % tif file
    a = imread(fullfile(pathname, filename));
elseif strcmpi(ext9, '.h5')   % hdf5 data
    if contains(fn9, 'master')
        fprintf('%s is a master file, will not be proecssed!\n',filename);       
    else
        dat9 = readHdf5(fullfile(pathname, filename), '12-ID-B');
        a = dat9.data';
    end
else
    fprintf('%s is not supported and proecssed!\n',cfile);
end  

clim1=min(a(:));
clim2=max(a(:));

a = double(a);%a = flipud(a);
figh = figure;
imagesc(a);
axis image; %axis xy;
set(handles.ed_clim1, 'string', num2str(clim1));
set(handles.ed_clim2, 'string', num2str(clim2));
hold on;
setappdata(gcbf, 'imagehandle', figh);
setappdata(gcbf, 'image', a);

hMenubar = findall(figh,'tag','figMenuFile');
set(findall(hMenubar), 'visible', 'off')
hMenubar = findall(figh,'tag','figMenuView');
set(findall(hMenubar), 'visible', 'off')
hMenubar = findall(figh,'tag','figMenuEdit');
set(findall(hMenubar), 'visible', 'off')
hMenubar = findall(figh,'tag','figMenuInsert');
set(findall(hMenubar), 'visible', 'off')
hMenubar = findall(figh,'tag','figMenuTools');
set(findall(hMenubar), 'visible', 'off')
hMenubar = findall(figh,'tag','figMenuDesktop');
set(findall(hMenubar), 'visible', 'off')
hMenubar = findall(figh,'tag','figMenuWindow');
set(findall(hMenubar), 'visible', 'off')
hMenubar = findall(figh,'tag','figMenuHelp');
set(findall(hMenubar), 'visible', 'off')
%get(findall(hMenubar),'tag')
hToolbar = findall(figh,'tag','FigureToolBar');
hStandardTools = findall(hToolbar,'-regexp', 'tag', 'Standard');
set(hStandardTools, 'visible', 'off');
hMenuFile = uimenu(figh,...
    'Label','Scale',...
    'Position',1,...
    'HandleVisibility','callback',...
    'Tag','SSL_MenuScale');
hMenuScale1 = uimenu(hMenuFile,...
    'Label','50%...',...
    'Position',1,...
    'HandleVisibility','callback',...
    'Tag','SSL_MenuFileOpen',...
    'callback',{@setcolor, 0.5});
hMenuScale2 = uimenu(hMenuFile,...
    'Label','10%...',...
    'Position',2,...
    'HandleVisibility','callback',...
    'Tag','SSL_MenuFileOpen',...
    'callback',{@setcolor, 0.1});
hMenuScale3 = uimenu(hMenuFile,...
    'Label','5%...',...
    'Position',3,...
    'HandleVisibility','callback',...
    'Tag','SSL_MenuFileOpen',...
    'callback',{@setcolor, 0.05});
hMenuScale4 = uimenu(hMenuFile,...
    'Label','1%...',...
    'Position',4,...
    'HandleVisibility','callback',...
    'Tag','SSL_MenuFileOpen',...
    'callback',{@setcolor, 0.01});
hMenuScale5 = uimenu(hMenuFile,...
    'Label','0.1%...',...
    'Position',5,...
    'HandleVisibility','callback',...
    'Tag','SSL_MenuFileOpen',...
    'callback',{@setcolor, 0.001});
%get(findall(hToolbar),'tag')
function setcolor(varargin)
    f = varargin{3};
    ax = findobj(gcbf, 'type', 'axes');
    img = get(findobj(ax, 'type', 'image'), 'CData');
    mv = max(img(:));
    set(ax, 'CLim', [0, f*mv]);
    

% --- Executes on button press in pb_findpeaks.
function pb_findpeaks_Callback(hObject, eventdata, handles)
% hObject    handle to pb_findpeaks (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
h = getappdata(gcbf, 'imagehandle');
hax = findobj(h, 'type', 'axes');
%a = getappdata(gcbf, 'image');
a = get(findobj(hax, 'type', 'image'), 'CData');
a = double(a);
if size(a, 1) < 1000
    p = fc_vertical(a, str2double(get(handles.ed_step, 'string')));
    %p = fc_vertical2(a, s.center);
else
    s = getgihandle;
    p = fc_vertical2(a, s.center);
    t = cellfun(@numel, p);
    p(t<30) = [];
    for i=1:numel(p)
        p{1} = p{1}(:,1:2);
    end
end
%    str2double(get(handles.ed_threshhold, 'string')));
setappdata(gcbf, 'allpeaks', p);
disp('Peaks are found')

function ed_step_Callback(hObject, eventdata, handles)
% hObject    handle to ed_step (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of ed_step as text
%        str2double(get(hObject,'String')) returns contents of ed_step as a double


% --- Executes during object creation, after setting all properties.
function ed_step_CreateFcn(hObject, eventdata, handles)
% hObject    handle to ed_step (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function ed_threshhold_Callback(hObject, eventdata, handles)
% hObject    handle to ed_threshhold (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of ed_threshhold as text
%        str2double(get(hObject,'String')) returns contents of ed_threshhold as a double


% --- Executes during object creation, after setting all properties.
function ed_threshhold_CreateFcn(hObject, eventdata, handles)
% hObject    handle to ed_threshhold (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in pb_showandhide.
function pb_showandhide_Callback(hObject, eventdata, handles)
% hObject    handle to pb_showandhide (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
p = getappdata(gcbf, 'allpeaks');
h = getappdata(gcbf, 'imagehandle');
hax = findobj(h, 'type', 'axes');
if strcmp(get(hObject, 'string'), 'Show All')
    set(hObject, 'string', 'Hide All');
    k = [];
    for i=1:numel(p)
        k = [k;plot(p{i}(:,1), p{i}(:,2), 'yo', 'parent', hax)];
    end
    setappdata(gcbf, 'allpeakhandle', k);
else
    set(hObject, 'string', 'Show All');
    k = getappdata(gcbf, 'allpeakhandle');
    delete(k);
    setappdata(gcbf, 'allpeakhandle', []);
end

% --- Executes on button press in pb_OK.
function pb_OK_Callback(hObject, handles)
% hObject    handle to pb_OK (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
k = getappdata(gcbf, 'currentP');
p = getappdata(gcbf, 'allpeaks');
t = getappdata(gcbf,'currentPhandle');
peak = getappdata(gcbf, 'peak');
if isempty(peak)
    peak = {[],[],[],[],[],[],[],[], []};
end
%h = handles.uipanel_select.SelectedObject;
%h = get(handles.uipanel_select, 'SelectedObject');
h = hObject;

%h = getappdata(handles.uipanel_selected, 'selected');
peakN = get(h, 'string');peakN = str2double(peakN(5));
if isfield(peak{peakN}, 'groupN')
    gN = unique([peak{peakN}.groupN, k]);
else
    gN = k;
end
%peak{peakN}.position = [peak{peakN}; p{k}];
peak{peakN}.groupN = gN;
peak{peakN}.position = cat(1, p{gN});
peak{peakN}.q = str2double(get(handles.(['edit_q', num2str(peakN)]), 'string'));
%set(t, 'ForegroundColor', get(h, 'ForegroundColor'));
set(t, 'color', get(h, 'ForegroundColor'));
setappdata(gcbf, 'peak', peak);
assignin('base', 'peak', peak);
t = 1;data = [];
for i=1:numel(peak)
    if ~isempty(peak{i})
        data(t).q = peak{i}.q;
        data(t).xy = peak{i}.position';
        data(t).xy(2,:) = 619 - data(t).xy(2,:);
        %data(t).xy(2,:) = data(t).xy(2,:);
        t = t+1;
    end
end
try
    qCalData = evalin('base', 'qCalData');
catch
    qCalData = [];
end
%if ~isempty(qCalData);
    qCalData.calibrant.wavelength = eng2wl(14);
    qCalData.calibrant.data = data;
    assignin('base', 'qCalData', qCalData);
%else
%    assignin('base', 'data', data);
%end


% --- Executes on button press in pb_showforassign.
function pb_showforassign_Callback(hObject, eventdata, handles)
% hObject    handle to pb_showforassign (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
k = getappdata(gcbf, 'currentP');
p = getappdata(gcbf, 'allpeaks');
h = getappdata(gcbf, 'imagehandle');
hax = findobj(h, 'type', 'axes');

if isempty(k)
    k = 1;
else
    k = k + 1;
end
if k>numel(p)
    k=0;
    disp('All done')
    return
end
t = plot(p{k}(:,1), p{k}(:,2), 'wo', 'parent', hax);
setappdata(gcbf,'currentP', k);
setappdata(gcbf,'currentPhandle', t);


function radiocallback(hObject, handles)
t = getappdata(gcbf,'currentPhandle');
h = hObject;
set(t, 'color', get(h, 'ForegroundColor'));
pb_OK_Callback(hObject, handles)


% --- Executes when selected object is changed in uipanel_select.
function uipanel_select_SelectionChangeFcn(hObject, eventdata, handles)
% hObject    handle to the selected object in uipanel_select 
% eventdata  structure with the following fields (see UIBUTTONGROUP)
%	EventName: string 'SelectionChanged' (read only)
%	OldValue: handle of the previously selected object or empty if none was selected
%	NewValue: handle of the currently selected object
% handles    structure with handles and user data (see GUIDATA)
radiocallback(hObject, handles)
%setappdata(handles.uipanel_selected, 'selected', hObject);
%eventdata.NewValue
%hObject
%get(handles.uipanel_select, 'SelectedObject')


% --- Executes on button press in radiobutton2.
function radiobutton2_Callback(hObject, eventdata, handles)
% hObject    handle to radiobutton2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of radiobutton2
radiocallback(hObject, handles)

% --- Executes on button press in pb_clear.
function pb_clear_Callback(hObject, eventdata, handles)
% hObject    handle to pb_clear (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
setappdata(gcbf, 'currentP', 0);
p = getappdata(gcbf, 'allpeaks');
h = getappdata(gcbf, 'imagehandle');
setappdata(gcbf, 'peak', []);
hax = findobj(h, 'type', 'axes');
k = findobj(hax, 'type', 'line');
delete(k);
assignin('base', 'peak', []);



% --- Executes on key press with focus on figure1 and none of its controls.
function figure1_KeyPressFcn(hObject, eventdata, handles)
% hObject    handle to figure1 (see GCBO)
% eventdata  structure with the following fields (see FIGURE)
%	Key: name of the key that was pressed, in lower case
%	Character: character interpretation of the key(s) that was pressed
%	Modifier: name(s) of the modifier key(s) (i.e., control, shift) pressed
% handles    structure with handles and user data (see GUIDATA)
eventdata.Key
eventdata.Character
eventdata.Modifier
if strcmp(eventdata.Modifier, 'alt') & strcmp(eventdata.Character, 'a') 
    pb_showforassign_Callback;
end
if strcmp(eventdata.Modifier, 'alt') & strcmp(eventdata.Character, 's') 
    pb_OK_Callback(hObject, handles)
end


% --- Executes on button press in pb_plotagbh.
function pb_plotagbh_Callback(hObject, eventdata, handles)
% hObject    handle to pb_plotagbh (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
contents = cellstr(get(handles.pm_reference,'String'));
refsample = contents{get(handles.pm_reference,'Value')};
switch refsample
    case 'Ag Behenate'
%        Refpeak = [0.9684, 1.076, 1.1836, 1.372, 1.391, 1.640, 1.828, 2.372];
%        Refpeak = [0.9684, 1.076, 1.1836, 1.3686, 1.3839, 1.6346, 1.8119, 2.3676];
        % d-spacing = 58.380A; q=0.107626
        Refpeak = [0.8610 0.9686, 1.0763, 1.1838, 1.3686, 1.3839, 1.6346, 1.8119, 2.3676];
    case 'Al2O3/Si'
        Refpeak = [1.80495626128487,2.46244953614232,2.63966109615577,...
            3.01200160395548,3.60991252256974,3.92239500066561,4.14650316655710,...
            4.47205535779290,5.04295232788956,5.26452308463249,5.66778911854108,...
            5.82536698657524,6.01159818008358,6.29601567766764,6.54459914659945,...
            6.84450561956098,7.31708428864174,7.58651446143875];
end
% assign peaks in the edit boxes.
for i=1:numel(Refpeak)
    if i>9
        break
    end
    set(handles.(['edit_q', num2str(i)]), 'string', Refpeak(i));
end

agbhW = reference(refsample);
% Load agbh profile and plot
% Agbh data
%load agbhWAXS.mat;
figure;
scalef = 1/max(agbhW(:,2))*0.085;
agbhW(:,2) = agbhW(:,2)*scalef;
plot(agbhW(:,1), agbhW(:,2)/max(agbhW(:,2))*0.08, 'b'); 
set(gca, 'xlim', [0.9, max(agbhW(:,1))*0.9], 'ylim', [0.00, 0.1]);
xlabel(sprintf('q (%c^{-1})', char(197)), 'fontsize', 12);
title('WAXS of Silver Behenate')
for i=1:numel(Refpeak)
    x = [Refpeak(i), Refpeak(i)];
    temp = abs(agbhW(:,1)-Refpeak(i));
    [~, indmin] = min(temp);
    y = [0, agbhW(indmin(1), 2)];
    t = line(x, y, 'color' ,'r');set(t, 'LineStyle', ':');
    text(x(2), y(2)+0.001, sprintf('%0.4f', Refpeak(i)), 'rotation', 90)
end


function edit_q1_Callback(hObject, eventdata, handles)
% hObject    handle to edit_q1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_q1 as text
%        str2double(get(hObject,'String')) returns contents of edit_q1 as a double


% --- Executes during object creation, after setting all properties.
function edit_q1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_q1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit_q2_Callback(hObject, eventdata, handles)
% hObject    handle to edit_q2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_q2 as text
%        str2double(get(hObject,'String')) returns contents of edit_q2 as a double


% --- Executes during object creation, after setting all properties.
function edit_q2_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_q2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit_q3_Callback(hObject, eventdata, handles)
% hObject    handle to edit_q3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_q3 as text
%        str2double(get(hObject,'String')) returns contents of edit_q3 as a double


% --- Executes during object creation, after setting all properties.
function edit_q3_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_q3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit_q4_Callback(hObject, eventdata, handles)
% hObject    handle to edit_q4 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_q4 as text
%        str2double(get(hObject,'String')) returns contents of edit_q4 as a double


% --- Executes during object creation, after setting all properties.
function edit_q4_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_q4 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit_q5_Callback(hObject, eventdata, handles)
% hObject    handle to edit_q5 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_q5 as text
%        str2double(get(hObject,'String')) returns contents of edit_q5 as a double


% --- Executes during object creation, after setting all properties.
function edit_q5_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_q5 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit_q6_Callback(hObject, eventdata, handles)
% hObject    handle to edit_q6 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_q6 as text
%        str2double(get(hObject,'String')) returns contents of edit_q6 as a double


% --- Executes during object creation, after setting all properties.
function edit_q6_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_q6 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit_q8_Callback(hObject, eventdata, handles)
% hObject    handle to edit_q8 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_q8 as text
%        str2double(get(hObject,'String')) returns contents of edit_q8 as a double


% --- Executes during object creation, after setting all properties.
function edit_q8_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_q8 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in radiobutton1.
function radiobutton1_Callback(hObject, eventdata, handles)
% hObject    handle to radiobutton1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
radiocallback(hObject, handles)


% --- Executes on button press in radiobutton3.
function radiobutton3_Callback(hObject, eventdata, handles)
% hObject    handle to radiobutton3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
radiocallback(hObject, handles)


% --- Executes on button press in radiobutton4.
function radiobutton4_Callback(hObject, eventdata, handles)
% hObject    handle to radiobutton4 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
radiocallback(hObject, handles)


% --- Executes on button press in radiobutton5.
function radiobutton5_Callback(hObject, eventdata, handles)
% hObject    handle to radiobutton5 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
radiocallback(hObject, handles)


% --- Executes on button press in radiobutton6.
function radiobutton6_Callback(hObject, eventdata, handles)
% hObject    handle to radiobutton6 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
radiocallback(hObject, handles)

function radiobutton7_Callback(hObject, eventdata, handles)
% hObject    handle to radiobutton6 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
radiocallback(hObject, handles)


% --- Executes on button press in radiobutton8.
function radiobutton8_Callback(hObject, eventdata, handles)
% hObject    handle to radiobutton8 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
radiocallback(hObject, handles)

% --- Executes on button press in radiobutton9.
function radiobutton9_Callback(hObject, eventdata, handles)
% hObject    handle to radiobutton9 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
radiocallback(hObject, handles)


% --- Executes during object creation, after setting all properties.
function figure1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to figure1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called


% --- Executes on button press in pbfindDetSetup.
function pbfindDetSetup_Callback(hObject, eventdata, handles)
% hObject    handle to pbfindDetSetup (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


peak = getappdata(gcbf, 'peak');
if isempty(peak)
    peak = evalin('base', 'peak');
end

if isempty(peak)
    peak = {[],[],[],[],[],[],[],[]};
end
Np = cellfun(@isempty, peak);
qv = zeros(1, sum(Np));
XY = cell(1, sum(Np));
k = 1;
for i=1:numel(peak)
    if ~Np(i)
        qv(k) = peak{i}.q;
        XY{k} = peak{i}.position;
        k = k+1;
    end
end

% other parameters not to fit
psize = str2double(get(handles.edit_psize, 'string'));
otherFixedParameters.pixelsize = psize;
xeng = str2double(get(handles.edit_energy, 'string'));
wl = eng2wl(xeng);
otherFixedParameters.waveln = wl;
% initial parameters to fit
p = zeros(1, 6);LB = p; UB = p;
for i=1:6
    p(i) = str2double(get(handles.(sprintf('edit_p%i', i)), 'string'));
    LB(i) = str2double(get(handles.(sprintf('LB_p%i', i)), 'string'));
    UB(i) = str2double(get(handles.(sprintf('UB_p%i', i)), 'string'));
end
%p = [centerX, centerY, sdd, Pitch, Yaw, Roll];
%LB = [XL, YL, sddL, PitchL, YawL, RollL];
%UB = [XU, YU, sddU, PitchU, YawU, RollU];
options = optimset('fminsearch');
%options = optimset(options, 'TolX',0.1E-8);
%    options = optimset(options, 'PlotFcns',@optimplotx);
%    options = optimset(options, 'OutputFcn',@BLFit_outfun);
options = optimset(options, 'MaxIter',1000);
options = optimset(options, 'MaxFunEvals', 2000);

[INLP, val] = fminsearchcon(@(x) fitDetSetup(x, XY, qv, otherFixedParameters), p,LB,UB, [], [], [], options);
fprintf('Fitting done with val = %0.5e\n', val);
for i=1:6
    set(handles.(sprintf('edit_p%i', i)), 'string', INLP(i));
end

h = getappdata(gcbf, 'imagehandle');
hax = findobj(h, 'type', 'axes');
xl = get(hax, 'xlim');
yl = get(hax, 'ylim');
if isempty(xl)
    return
end
k = findobj(hax, 'type', 'line');
delete(k);
drawnow;
% waveln = varargin{2};
% center = varargin{3};
% pixelsize =varargin{4};
% sdd =varargin{5};
for i=1:numel(XY)
%    qfit{i} = pixel2q(XY{i}, INLP(1:2), INLP(3), psize, INLP(4:6), wl);
    center = INLP(1:2);
    sdd = INLP(3);
    detangle = INLP(4:6);
    if isempty(XY{i})
        continue;
    end
    if size(XY{i}, 2) == 4
        XY{i}(:,3:4) = [];
    end
    pix = bsxfun(@minus,XY{i},center);
    qv0 = pixel2qv0(pix, wl, sdd, psize, detangle);
    qfit{i}  = sqrt(sum(qv0.*qv0, 2));
end

assignin('base', 'qfit', qfit)
assignin('base', 'qv', qv)
h = zeros(size(qv));
for i=1:numel(qv)
%    pos = q2pixel(qv(i), wl, INLP(1:2), psize, INLP(3), INLP(4:6));
    pos = q_powder2pixel(qv(i), wl, INLP(3), psize, INLP(4:6));
    pos = pos + repmat(INLP(1:2), numel(pos(:,1)), 1);
    x = pos(:, 1);
    y = pos(:, 2);
    t = (x>xl(1))&(x<xl(2));
    t = t & (y>yl(1))&(y<yl(2));
    if sum(t) > 0
        line(XY{i}(:,1), XY{i}(:,2),'marker', 'o', 'color', 'w', 'parent', hax)
        h(i) = line(x(t), y(t), 'color', 'r', 'parent', hax);
    end
end


function edit_p1_Callback(hObject, eventdata, handles)
% hObject    handle to edit_p1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_p1 as text
%        str2double(get(hObject,'String')) returns contents of edit_p1 as a double


% --- Executes during object creation, after setting all properties.
function edit_p1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_p1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit_p2_Callback(hObject, eventdata, handles)
% hObject    handle to edit_p2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_p2 as text
%        str2double(get(hObject,'String')) returns contents of edit_p2 as a double


% --- Executes during object creation, after setting all properties.
function edit_p2_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_p2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit_p3_Callback(hObject, eventdata, handles)
% hObject    handle to edit_p3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_p3 as text
%        str2double(get(hObject,'String')) returns contents of edit_p3 as a double


% --- Executes during object creation, after setting all properties.
function edit_p3_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_p3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit_p4_Callback(hObject, eventdata, handles)
% hObject    handle to edit_p4 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_p4 as text
%        str2double(get(hObject,'String')) returns contents of edit_p4 as a double


% --- Executes during object creation, after setting all properties.
function edit_p4_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_p4 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit_p5_Callback(hObject, eventdata, handles)
% hObject    handle to edit_p5 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_p5 as text
%        str2double(get(hObject,'String')) returns contents of edit_p5 as a double


% --- Executes during object creation, after setting all properties.
function edit_p5_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_p5 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit_p6_Callback(hObject, eventdata, handles)
% hObject    handle to edit_p6 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_p6 as text
%        str2double(get(hObject,'String')) returns contents of edit_p6 as a double


% --- Executes during object creation, after setting all properties.
function edit_p6_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_p6 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function LB_p1_Callback(hObject, eventdata, handles)
% hObject    handle to LB_p1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of LB_p1 as text
%        str2double(get(hObject,'String')) returns contents of LB_p1 as a double


% --- Executes during object creation, after setting all properties.
function LB_p1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to LB_p1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function LB_p2_Callback(hObject, eventdata, handles)
% hObject    handle to LB_p2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of LB_p2 as text
%        str2double(get(hObject,'String')) returns contents of LB_p2 as a double


% --- Executes during object creation, after setting all properties.
function LB_p2_CreateFcn(hObject, eventdata, handles)
% hObject    handle to LB_p2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function LB_p3_Callback(hObject, eventdata, handles)
% hObject    handle to LB_p3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of LB_p3 as text
%        str2double(get(hObject,'String')) returns contents of LB_p3 as a double


% --- Executes during object creation, after setting all properties.
function LB_p3_CreateFcn(hObject, eventdata, handles)
% hObject    handle to LB_p3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function LB_p4_Callback(hObject, eventdata, handles)
% hObject    handle to LB_p4 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of LB_p4 as text
%        str2double(get(hObject,'String')) returns contents of LB_p4 as a double


% --- Executes during object creation, after setting all properties.
function LB_p4_CreateFcn(hObject, eventdata, handles)
% hObject    handle to LB_p4 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function LB_p5_Callback(hObject, eventdata, handles)
% hObject    handle to LB_p5 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of LB_p5 as text
%        str2double(get(hObject,'String')) returns contents of LB_p5 as a double


% --- Executes during object creation, after setting all properties.
function LB_p5_CreateFcn(hObject, eventdata, handles)
% hObject    handle to LB_p5 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function LB_p6_Callback(hObject, eventdata, handles)
% hObject    handle to LB_p6 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of LB_p6 as text
%        str2double(get(hObject,'String')) returns contents of LB_p6 as a double


% --- Executes during object creation, after setting all properties.
function LB_p6_CreateFcn(hObject, eventdata, handles)
% hObject    handle to LB_p6 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function UB_p1_Callback(hObject, eventdata, handles)
% hObject    handle to UB_p1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of UB_p1 as text
%        str2double(get(hObject,'String')) returns contents of UB_p1 as a double


% --- Executes during object creation, after setting all properties.
function UB_p1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to UB_p1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function UB_p2_Callback(hObject, eventdata, handles)
% hObject    handle to UB_p2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of UB_p2 as text
%        str2double(get(hObject,'String')) returns contents of UB_p2 as a double


% --- Executes during object creation, after setting all properties.
function UB_p2_CreateFcn(hObject, eventdata, handles)
% hObject    handle to UB_p2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function UB_p3_Callback(hObject, eventdata, handles)
% hObject    handle to UB_p3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of UB_p3 as text
%        str2double(get(hObject,'String')) returns contents of UB_p3 as a double


% --- Executes during object creation, after setting all properties.
function UB_p3_CreateFcn(hObject, eventdata, handles)
% hObject    handle to UB_p3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function UB_p4_Callback(hObject, eventdata, handles)
% hObject    handle to UB_p4 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of UB_p4 as text
%        str2double(get(hObject,'String')) returns contents of UB_p4 as a double


% --- Executes during object creation, after setting all properties.
function UB_p4_CreateFcn(hObject, eventdata, handles)
% hObject    handle to UB_p4 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function UB_p5_Callback(hObject, eventdata, handles)
% hObject    handle to UB_p5 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of UB_p5 as text
%        str2double(get(hObject,'String')) returns contents of UB_p5 as a double


% --- Executes during object creation, after setting all properties.
function UB_p5_CreateFcn(hObject, eventdata, handles)
% hObject    handle to UB_p5 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function UB_p6_Callback(hObject, eventdata, handles)
% hObject    handle to UB_p6 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of UB_p6 as text
%        str2double(get(hObject,'String')) returns contents of UB_p6 as a double


% --- Executes during object creation, after setting all properties.
function UB_p6_CreateFcn(hObject, eventdata, handles)
% hObject    handle to UB_p6 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit_psize_Callback(hObject, eventdata, handles)
% hObject    handle to edit_psize (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_psize as text
%        str2double(get(hObject,'String')) returns contents of edit_psize as a double


% --- Executes during object creation, after setting all properties.
function edit_psize_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_psize (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit_energy_Callback(hObject, eventdata, handles)
% hObject    handle to edit_energy (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_energy as text
%        str2double(get(hObject,'String')) returns contents of edit_energy as a double


% --- Executes during object creation, after setting all properties.
function edit_energy_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_energy (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in pb_transfer.
function pb_transfer_Callback(hObject, eventdata, handles)
% hObject    handle to pb_transfer (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
psize = str2double(get(handles.edit_psize, 'string'));
xeng = str2double(get(handles.edit_energy, 'string'));
wl = eng2wl(xeng);
% initial parameters to fit
p = zeros(1, 6);
for i=1:6
    p(i) = str2double(get(handles.(sprintf('edit_p%i', i)), 'string'));
end

s.SDD = p(3);
s.center = p(1:2);
s.tiltangle = p(4:6);
s.ai = 0;
s.tthi = 0;
s.edensity = 0;
s.beta = 0;
s.xeng = xeng;
s.psize = psize;
s.waveln = wl;
s.Q = 0.10763;
px = SDDcal(s.SDD, s.Q, s.psize, s.waveln, 'sdd');
s.px = px;

setgihandle(s);
gisaxsleenew('SetSAXSValues')

function ref = reference(sample)
switch sample
    case 'Ag Behenate'
        ref = load('agbh.dat');
    case 'Al2O3/Si'
        ref = load('al2o3simixture.dat');
        ref(:,3) = [];
end
           
            


% --- Executes on selection change in pm_reference.
function pm_reference_Callback(hObject, eventdata, handles)
% hObject    handle to pm_reference (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns pm_reference contents as cell array
%        contents{get(hObject,'Value')} returns selected item from pm_reference


% --- Executes during object creation, after setting all properties.
function pm_reference_CreateFcn(hObject, eventdata, handles)
% hObject    handle to pm_reference (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in radiobutton10.
function radiobutton10_Callback(hObject, eventdata, handles)
% hObject    handle to radiobutton10 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of radiobutton10
figh = getappdata(gcbf, 'imagehandle');
a = getappdata(gcbf, 'image');
if get(hObject, 'value')
    a = flipud(a);
end
figure(figh);clf;
imagesc(a);axis image; 
axis xy;
resetImageColorScale(handles);
hold on;



function ed_clim1_Callback(hObject, eventdata, handles)
% hObject    handle to ed_clim1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
resetImageColorScale(handles)
% Hints: get(hObject,'String') returns contents of ed_clim1 as text
%        str2double(get(hObject,'String')) returns contents of ed_clim1 as a double


% --- Executes during object creation, after setting all properties.
function ed_clim1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to ed_clim1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function ed_clim2_Callback(hObject, eventdata, handles)
% hObject    handle to ed_clim2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
resetImageColorScale(handles)
% Hints: get(hObject,'String') returns contents of ed_clim2 as text
%        str2double(get(hObject,'String')) returns contents of ed_clim2 as a double


% --- Executes during object creation, after setting all properties.
function ed_clim2_CreateFcn(hObject, eventdata, handles)
% hObject    handle to ed_clim2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in pb_resetColorScales.
function pb_resetColorScales_Callback(hObject, eventdata, handles)
% hObject    handle to pb_resetColorScales (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
resetImageColorScale(handles)

function resetImageColorScale(handles)
clim1=str2num(get(handles.ed_clim1, 'string'));
clim2=str2num(get(handles.ed_clim2, 'string'));
figh = getappdata(gcbf, 'imagehandle');
hax = findobj(figh, 'type', 'axes');
set(hax, 'Clim', [clim1 clim2]);


% --- Executes on button press in radiobutton10.
function pushbutton13_Callback(hObject, eventdata, handles)
% hObject    handle to radiobutton10 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
radiocallback(hObject, handles)


function edit_q9_Callback(hObject, eventdata, handles)
% hObject    handle to edit_q9 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_q9 as text
%        str2double(get(hObject,'String')) returns contents of edit_q9 as a double


% --- Executes during object creation, after setting all properties.
function edit_q9_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_q9 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
