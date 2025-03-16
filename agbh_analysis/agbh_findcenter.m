function varargout = agbh_findcenter(varargin)
% AGBH_FINDCENTER MATLAB code for agbh_findcenter.fig
%      AGBH_FINDCENTER, by itself, creates a new AGBH_FINDCENTER or raises the existing
%      singleton*.
%
%      H = AGBH_FINDCENTER returns the handle to a new AGBH_FINDCENTER or the handle to
%      the existing singleton*.
%
%      AGBH_FINDCENTER('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in AGBH_FINDCENTER.M with the given input arguments.
%
%      AGBH_FINDCENTER('Property','Value',...) creates a new AGBH_FINDCENTER or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before agbh_findcenter_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to agbh_findcenter_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help agbh_findcenter

% Last Modified by GUIDE v2.5 09-Jul-2019 02:56:20

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @agbh_findcenter_OpeningFcn, ...
                   'gui_OutputFcn',  @agbh_findcenter_OutputFcn, ...
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


% --- Executes just before agbh_findcenter is made visible.
function agbh_findcenter_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to agbh_findcenter (see VARARGIN)

% Choose default command line output for agbh_findcenter
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes agbh_findcenter wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = agbh_findcenter_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on button press in pb_AutoAnalysis.
function pb_AutoAnalysis_Callback(hObject, eventdata, handles)
% hObject    handle to pb_AutoAnalysis (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
    hSAXSimageview = evalin('base', 'SAXSimagehandle');
%    cut.handle = line(cut.X, cut.Y, 'parent', hSAXSimageview);
    saxs = getgihandle;
    sf = str2double(get(handles.ed_ScaleFactor, 'string'));
    [center, SDDpeak] = agbhSAXS(saxs.image/sf, 0);
    
    % draw a circle...
    %[xc,yc,Re,a] = circfit(handles.selectedX, handles.selectedY);
    th = linspace(0, 2*pi, 100);
    xc = center(1); yc = center(2);
    Re = SDDpeak;
    xe = Re*cos(th)+yc; ye = Re*sin(th)+xc;
    tmph = line(xe,ye, 'parent', hSAXSimageview, 'color', 'g');
        %tmph = line([xe;xe(1)],[ye;ye(1)],'-', 'parent', handles.ImageAxes);

    try
        sampleinfo = evalin('base', 'sampleinfo');
        saxs.ai = sampleinfo.th;
        saxs.xeng = sampleinfo.Energy;
        saxs.waveln = eng2wl(saxs.xeng);
    catch
        fprintf('sampleinfo not available\n');
    end    
    try
        saxs.center = [yc, xc];
        saxs.px = Re;
        setgihandle(saxs);
        gisaxsleenew('SetSAXSValues')
    catch
        error('Error in putting the setup info to gisaxsleenew');
    end

function ed_ScaleFactor_Callback(hObject, eventdata, handles)
% hObject    handle to ed_ScaleFactor (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of ed_ScaleFactor as text
%        str2double(get(hObject,'String')) returns contents of ed_ScaleFactor as a double


% --- Executes during object creation, after setting all properties.
function ed_ScaleFactor_CreateFcn(hObject, eventdata, handles)
% hObject    handle to ed_ScaleFactor (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in pb_centerpicking.
function pb_centerpicking_Callback(hObject, eventdata, handles)
% hObject    handle to pb_centerpicking (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
    hSAXSimageview = evalin('base', 'SAXSimagehandle');
%    cut.handle = line(cut.X, cut.Y, 'parent', hSAXSimageview);
    saxs = getgihandle;
    sf = str2double(get(handles.ed_ScaleFactor, 'string'));
    [center, SDDpeak] = agbhSAXS(saxs.image/sf, 2);
    
    % draw a circle...
    %[xc,yc,Re,a] = circfit(handles.selectedX, handles.selectedY);
    th = linspace(0, 2*pi, 100);
    xc = center(1); yc = center(2);
    Re = SDDpeak;
    xe = Re*cos(th)+yc; ye = Re*sin(th)+xc;
    tmph = line(xe,ye, 'parent', hSAXSimageview, 'color', 'g');
        %tmph = line([xe;xe(1)],[ye;ye(1)],'-', 'parent', handles.ImageAxes);

    try
        sampleinfo = evalin('base', 'sampleinfo');
        saxs.ai = sampleinfo.th;
        saxs.xeng = sampleinfo.Energy;
        saxs.waveln = eng2wl(saxs.xeng);
    catch
        fprintf('sampleinfo not available\n');
    end    
    try
        saxs.center = [yc, xc];
        saxs.px = Re;
        setgihandle(saxs);
        gisaxsleenew('SetSAXSValues')
    catch
        error('Error in putting the setup info to gisaxsleenew');
    end


% --- Executes on button press in pb_startfindcenter.
function pb_startfindcenter_Callback(hObject, eventdata, handles)
% hObject    handle to pb_startfindcenter (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
APS_findcenter


% --- Executes on button press in pushbutton5.
function pushbutton5_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton5 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
findcenter_graphical
