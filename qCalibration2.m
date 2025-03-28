

%%
% "qCalibration2" performs detector q value calibration and 2-D image data
% processing for isotropic x-ray scattering images.
%
%#Data and parameters' storage and passing:
%  All internal data are stored in "qCalData" that is also written to the base 
%   workspace in the same name. When the data/parameters are needed to changed, 
%   the program first reads "qCalData" from base workspace, then modify through 
%   the interface, and write to the base workspace after modification done.
%
%#Data Structure:
%  "qCalData" has THREE fields: "userData", "Scatt", and "calibrant".  
%  "userData" stores all parameters at the GUI and has the following fileds: 
%     xeng, wavelength, pSize, detRow, detCol, isPilatus2m, offset, limits, 
%     beamx, beamxLb, beamxUb, beamy, beamyLb, beamyUb, SDD, SDDLb, SDDUb, 
%     yaw, yawLb, yawUb,  # parameters to fit
%     qNum, qMin, qMax, dq1, q2,dq2, q3, dq3, combQs, qChoice   # q setup
%
%   "isPilatus2m": 1 for pilatus2m, 0 for pilatus300k, 2 for PE detector, 3 for Mar300,  -1 for others
%
%  "Scatt" stores data, parameters and info necessary for processing scattering images. 
%   It has following fields: 
%     eng, wavelength, BeamXY (= [centerx centery]) , yaw, SDD, SDDrel, 
%     offset, limits, isPilatus, row, col, pSize, qNum, qMin, qMax, combQs, 
%     qChoice, qMode, qArray (1D), qMap (2D), qCMap (2D), qRMap, mask.
%     #"Scatt" will be saved in xxxsetup.mat file, which can be loaded later for data processing. 
%
%  "calibrant" stores standards' ring [x y]_j & q_j data, which will be used to
%   calibrate the beam center position (i.e. BeamXY) and sample-to-detector distance (SDD). 
%   It has two fields: calibrant.wavelength, and calibrant.data.  
%   calibrant.data is a struct array, and each row has two fields, xy, and q. 
%   Here is an example for manupilating calibrant.data: 
%      cdata =[]; cdata.xy=[x;y]; # x,y are coordinates of some points on a ring       
%      cdata.q = q;               # and q is the q value of this diffraction ring.
%      calibrant.data = [calibrant.data cdata]
%
%  Some other functions can be carried out with the current program: 
%    (1) customized qArray: users can make their own qArray in setup file,
%    and load into the program, and re-"Calculate Q Index Map" to make it
%    into effect.
%
% Known issues:
%   fixed. 1. Comb Linear & comb Log have no dynamic setting, need set parameters
%    before select!
%
% Written by Xiaobing Zuo, APS, Argonne National Laboratory  2011-
%%

function varargout = qCalibration2(varargin)
% QCALIBRATION2 MATLAB code for qCalibration2.fig
%      QCALIBRATION2, by itself, creates a new QCALIBRATION2 or raises the existing
%      singleton*.
%
%      H = QCALIBRATION2 returns the handle to a new QCALIBRATION2 or the handle to
%      the existing singleton*.
%
%      QCALIBRATION2('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in QCALIBRATION2.M with the given input arguments.
%
%      QCALIBRATION2('Property','Value',...) creates a new QCALIBRATION2 or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before qCalibration2_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to qCalibration2_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help qCalibration2

% Last Modified by GUIDE v2.5 10-Mar-2025 10:20:34

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @qCalibration2_OpeningFcn, ...
                   'gui_OutputFcn',  @qCalibration2_OutputFcn, ...
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

global detQCalibApp ;




% --- Executes just before qCalibration2 is made visible.
function qCalibration2_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to qCalibration2 (see VARARGIN)

% Choose default command line output for qCalibration2
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);



%fig= gcbf;
%setfield(fig, 'scatt', []);
%handles.calibrant.data=[];
%handles.calibrant.wavelength=0.;
%handles.userData=[];
%handles.scatt=[];
%h=findObj();




% UIWAIT makes qCalibration2 wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = qCalibration2_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;

qCalData=[];
qCalData.scatt=[];
%qCalData.calibrant=[];
qCalData.calibrant.wavelength=0.;
qCalData.calibrant.data=[];
qCalData.userData=[];
%handles.output.userData=qCalData;
assignin('base','qCalData', qCalData);
updateAll(handles);



function edit2_Callback(hObject, eventdata, handles)
% hObject    handle to edit2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit2 as text
%        str2double(get(hObject,'String')) returns contents of edit2 as a double


% --- Executes during object creation, after setting all properties.



function edit4_Callback(hObject, eventdata, handles)
% hObject    handle to edit4 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit4 as text
%        str2double(get(hObject,'String')) returns contents of edit4 as a double


function xc_init_Callback(hObject, eventdata, handles)
% hObject    handle to xc_init (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of xc_init as text
%        str2double(get(hObject,'String')) returns contents of xc_init as a double


% --- Executes during object creation, after setting all properties.

function edit8_Callback(hObject, eventdata, handles)
% hObject    handle to edit8 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit8 as text
%        str2double(get(hObject,'String')) returns contents of edit8 as a double





function yaw_init_Callback(hObject, eventdata, handles)
% hObject    handle to yaw_init (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of yaw_init as text
%        str2double(get(hObject,'String')) returns contents of yaw_init as a double




function edit10_Callback(hObject, eventdata, handles)
% hObject    handle to edit10 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit10 as text
%        str2double(get(hObject,'String')) returns contents of edit10 as a double



function edit11_Callback(hObject, eventdata, handles)
% hObject    handle to edit11 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit11 as text
%        str2double(get(hObject,'String')) returns contents of edit11 as a double




% --- Executes on selection change in popupmenu3.
function popupmenu3_Callback(hObject, eventdata, handles)
% hObject    handle to popupmenu3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns popupmenu3 contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popupmenu3



% --- Executes on button press in pbLoadSetup.
function pbLoadSetup_Callback(hObject, eventdata, handles)
% hObject    handle to pbLoadSetup (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% maybe need more work on combined linear array!

[filename, pathname] = uigetfile( ...
       {'*.m;*.fig;*.mat;*.mdl', 'All Setup Files (*.mat)'}, ...
        'Open as');
setupFile = [pathname, filename]; 
if setupFile==0 return;end

vars = whos('-file',setupFile);   
tmp=load(setupFile);

if ismember('saxs',{vars.name})
    cScatt = tmp.saxs;
elseif ismember('waxs',{vars.name})
    cScatt = tmp.waxs;
elseif ismember('laxs',{vars.name})
    cScatt = tmp.laxs;    
else
    disp('Should not be here, in pbLoadSetup_Callback()!');
end

syncScatt(handles, cScatt);
% 
% isPilatus2m = get(handles.rbPilatus2m, 'Value');
% 
% if isPilatus2m==1
%     if ismember('saxs',{vars.name})
%         %setupFile
%         tmp=load(setupFile);
%         cScatt = tmp.saxs;
%         syncScatt(handles, cScatt);
%     else
%         disp('This is NOT SAXS setup file!!!')
%     end
% elseif isPilatus2m==0
%    if ismember('waxs',{vars.name})
%         tmp=load(setupFile);
%         cScatt = tmp.waxs;
%         syncScatt(handles, cScatt);
%     else
%         disp('This is NOT WAXS setup file!!!')
%    end
% else
%     disp('Place holder for other detector!');
% end

function syncScatt(handles, cScatt)
%
[qCalData, userData, scatt, calibrant]= getData;
scatt=cScatt;

setDetectorInitials(handles, scatt.isPilatus2m);

% label for Pilatus 2M, 300, or PE
if scatt.isPilatus2m ==1        % Pilatus2m
    set(handles.rbPilatus2m, 'Value', 1);
    scatt.detector = 'Pilatus2M';
elseif scatt.isPilatus2m ==0    % Pilatus300
    set(handles.rbPilatus300, 'Value', 1);
    scatt.detector = 'Pilatus300K';
elseif scatt.isPilatus2m ==2    % PE
    set(handles.rbPE, 'Value', 1); 
    scatt.detector = 'PEXRpad';
elseif scatt.isPilatus2m ==3    % Pilatus900K
    set(handles.rbPilatus900k, 'Value', 1); 
    %scatt.detector = 'Mar300';
    scatt.detector = 'Pilatus900K';
elseif scatt.isPilatus2m ==4    % Eiger9M
    set(handles.rbEiger9m, 'Value', 1); 
    scatt.detector = 'Eiger9M';
elseif scatt.isPilatus2m ==-1  
    set(handles.rbOther, 'Value', 1);
    set(handles.edRow, 'String', scatt.row);
    set(handles.edCol, 'String', scatt.col);
    scatt.detector = 'Undefined';
else
    disp('should not be here, in syncScatt()!')
end

%setDetectorInitials(handles, scatt.isPilatus2m)

set(handles.edEnergy, 'String', scatt.eng);
set(handles.edPixelsize, 'String', scatt.pSize);
set(handles.edBeamX, 'String', scatt.BeamXY(1));
set(handles.edBeamY, 'String', scatt.BeamXY(2));
set(handles.edSDD, 'String', scatt.SDD);
set(handles.edYaw, 'String', scatt.yaw);
set(handles.edOffset, 'String', scatt.offset);
set(handles.edLimitLb, 'String', scatt.limits(1));
set(handles.edLimitUb, 'String', scatt.limits(2));
%if exist(handles.edAbsCoeff) & isfield(scatt, 'absIntCoeff')
if isfield(scatt, 'absIntCoeff')
    set(handles.edAbsCoeff, 'String', scatt.absIntCoeff);
end
set(handles.edPtsNum, 'String', scatt.qNum);
set(handles.edQMin, 'String', scatt.qMin);
set(handles.edQMax, 'String', scatt.qMax);
if isfield(scatt, 'qMode')
    switch(scatt.qMode)
        case 'Linear'
            set(handles.rbLinearQ, 'Value', 1);
        case 'LogQ'
            set(handles.rbLogQ, 'Value', 1);
        case 'Pixel'
            set(handles.rbPixelBasedQ, 'Value', 1);
        case 'CombLinear'
            set(handles.rbCombLinearQ, 'Value', 1);
            dq1= 0;
            q2 = 0.;
            dq2 =0.;
            q3 = 0.;
            dq3 = 0.;
            if isfield(scatt,'combQs')
                dq1 = scatt.combQs(1);
                q2  = scatt.combQs(2);
                dq2 = scatt.combQs(3);
                q3  = scatt.combQs(4);
                dq3 = scatt.combQs(5);
            end
            set(handles.edDQ1, 'String', num2str(dq1));
            set(handles.edQ2, 'String', num2str(q2));
            set(handles.edDQ2, 'String', num2str(dq2));
            set(handles.edQ3, 'String', num2str(q3));
            set(handles.edDQ3, 'String', num2str(dq3));
        case 'CombLogLinear'
            set(handles.combLogLinear, 'Value', 1);
            logQ2= 0;
            logPts = 0.;
            logDQ2 =0.;
            logQ3 = 0.;
            logDQ3 = 0.;
            if isfield(scatt,'combQs')
                logQ2 = scatt.combQs(1);
                logPts  = scatt.combQs(2);
                logDQ2 = scatt.combQs(3);
                logQ3  = scatt.combQs(4);
                logDQ3 = scatt.combQs(5);
            end
            set(handles.cLogQ2,  'String', num2str(logQ2));
            set(handles.cLogPts, 'String', num2str(logPts));
            set(handles.cLogDQ2, 'String', num2str(logDQ2));
            set(handles.cLogQ3,  'String', num2str(logQ3));
            set(handles.cLogDQ3, 'String', num2str(logDQ3));        
    end
end
if isfield(scatt, 'absIntStandardIndex')
    set(handles.absIntStandard, 'Value', scatt.absIntStandardIndex);
else
    set(handles.absIntStandard, 'Value', 1);
    scatt.absIntStandardIndex = 1;
    scatt.absIntStandard = 'None';
    userData.absIntCalibrantIndex = 1;
    userData.absIntCalibrant = 'None';    
end

if isfield(scatt, 'qStandardIndex')
    set(handles.qStandard, 'Value', scatt.qStandardIndex);
else
    set(handles.qStandard, 'Value', 1);
    scatt.qStandardIndex = 1;
    scatt.qStandard = 'None';
end

if isfield(scatt, 'BeamlineIndex')
    set(handles.Beamline, 'Value', scatt.BeamlineIndex);
else
    set(handles.Beamline, 'Value', 1);
    scatt.BeamlineIndex = 1;
    scatt.Beamline = '12-ID-B, APS';
end
setData(qCalData, userData, scatt, calibrant);
updateUserdata(handles);
        


% --- Executes on button press in pbSaveSetup.
function pbSaveSetup_Callback(hObject, eventdata, handles)
% hObject    handle to pbSaveSetup (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
qCalData = evalin('base','qCalData');

if qCalData.scatt.isPilatus2m == 1   % saxs, Pilatus2m
    saxs = qCalData.scatt;
    [filename, pathname] = uiputfile('saxssetup.mat', 'Save file name(*.mat)');
    if filename==0
        return;
    else
        %save([pathname, filename], 'saxs');
        save([pathname, filename], 'saxs', '-v7.3', '-nocompression');
        disp('setup file is saved!')
    end
elseif qCalData.scatt.isPilatus2m == 4   % saxs, Eiger 9M
    saxs = qCalData.scatt;
    [filename, pathname] = uiputfile('saxssetup.mat', 'Save file name(*.mat)');
    if filename==0
        return;
    else
        %save([pathname, filename], 'saxs');
        %save([pathname, filename], 'saxs', '-v7.3', '-nocompression');
        save([pathname, filename], 'saxs', '-v7.3');
        disp('setup file is saved!')
    end    
elseif qCalData.scatt.isPilatus2m == 0   % waxs, Pilatus300
    waxs = qCalData.scatt;
    [filename, pathname] = uiputfile( 'waxssetup.mat', 'Save file name(*.mat)');
    if filename==0 return;
    else
        %save([pathname, filename], 'waxs');
        save([pathname, filename], 'waxs', '-v7.3', '-nocompression');
        disp('setup file is saved!')
    end       
elseif qCalData.scatt.isPilatus2m == 2   % laxs, PE detector
    laxs = qCalData.scatt;
    [filename, pathname] = uiputfile( 'PEsetup.mat', 'Save file name(*.mat)');
    if filename==0 return;
    else
        %save([pathname, filename], 'laxs');
        %save([pathname, filename], 'laxs', '-v7.3', '-nocompression');
        save([pathname, filename], 'laxs', '-v7.3');
        disp('setup file is saved!')
    end       

elseif qCalData.scatt.isPilatus2m == 3 % for Pilatus900k detector
    waxs = qCalData.scatt;
    [filename, pathname] = uiputfile('waxssetup.mat', 'Save file name(*.mat)');
    if filename==0
        return;
    else
        %save([pathname, filename], 'saxs');
        save([pathname, filename], 'waxs', '-v7.3', '-nocompression');
        disp('setup file is saved!')
    end    
    
elseif qCalData.scatt.isPilatus2m == -1 % for other detector
    saxs = qCalData.scatt;
    [filename, pathname] = uiputfile('saxssetup.mat', 'Save file name(*.mat)');
    if filename==0
        return;
    else
        save([pathname, filename], 'saxs');
        disp('setup file is saved!')
    end    
else
    disp('why I am here? please check in pbSaveSetup_Callback....')
end

%disp(pathname);
%disp(filename);
% save seperate setup parameters in text and mask file
setupFolder = fullfile(pathname, 'SetupFiles');
setupText = [filename(1:end-4) '_param.txt'];

if ~exist(setupFolder, 'dir')
    mkdir(setupFolder);
end
myStr=evalc('disp(qCalData.scatt)');
fid = fopen(fullfile(setupFolder, setupText), 'w');
fprintf(fid, '%s', myStr);
fclose(fid);

if isfield(qCalData.scatt, 'mask')
    setupMask = fullfile(setupFolder, [filename(1:end-4) '_mask.bmp']);
    imwrite(logical(qCalData.scatt.mask), setupMask);
end


% --- Executes on button press in pbInitilize.
function pbInitilize_Callback(hObject, eventdata, handles)
% hObject    handle to pbInitilize (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
qCalData = evalin('base', 'qCalData');
calibrant = qCalData.calibrant;
calibrant.data=[];
calibrant.wavelength=qCalData.userData.wavelength;
qCalData.calibrant = calibrant;
assignin('base','qCalData', qCalData);




% --- Executes on button press in pbGetptsQ.
function pbGetptsQ_Callback(hObject, eventdata, handles)
% hObject    handle to pbGetptsQ (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

qCalData=evalin('base','qCalData');
calibrant = qCalData.calibrant;
try 
    mousePts = evalin('base', 'mousepnt');
    %cImgName = evalin('base', 'cImgName');
    %cImg = imgUpsideDn(double(imread(cImgName)));
    cImgName2 = getgihandle;
    cImg = cImgName2.image;
    [row, col] = size(cImg);
    bnds =[col row];
    halfW = 30;
    pks=[];
    rsqVal = 0.95;
    %for kk = 1:length(mousePts(:,1))
    for kk = 1:length(mousePts.x)
        x=round(mousePts.x(kk));
        y=round(mousePts.y(kk));
        %x=round(mousePts(kk,1));
        %y=round(mousePts(kk,2));        
        cpt = [x y];
        pkV = findVertPeak(cpt, cImg, halfW, bnds, rsqVal);
        pkH = findHorzPeak(cpt, cImg, halfW, bnds, rsqVal);
        newPKs = [pkV pkH];
        if isempty(newPKs) 
            newPKs = [x;y];
        end
        pks = [pks newPKs];
    end
    
    cdat.xy = sortrows(pks')';
   
    %cdat.xy = [mousePts.x; mousePts.y];
    
    cq = input('Please enter the q-value for the selected ring:');
    cdat.q=cq;

    calibrant.data = [calibrant.data cdat]
    qCalData.calibrant = calibrant;
    assignin('base','qCalData', qCalData);
catch
    disp('Pick calibrant pixels in SAXSimageviewer window ....');
    calibrant.data = [];
    calibrant.wavelength = 0.;
    qCalData.calibrant.wavelength = qCalData.userData.wavelength;
    assignin('base','qCalData', qCalData);
end




function edEnergy_Callback(hObject, eventdata, handles)
% hObject    handle to edEnergy (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edEnergy as text
%        str2double(get(hObject,'String')) returns contents of edEnergy as a double
updateAll(handles);
qCalData=evalin('base', 'qCalData');
qCalData.calibrant.wavelength = qCalData.userData.wavelength;
assignin('base','qCalData', qCalData);




function edPixelsize_Callback(hObject, eventdata, handles)
% hObject    handle to edPixelsize (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edPixelsize as text
%        str2double(get(hObject,'String')) returns contents of edPixelsize as a double
updateAll(handles);
%handles.pSize = str2double(get(handles.edPixelsize, 'string'));



function edBeamX_Callback(hObject, eventdata, handles)
% hObject    handle to edBeamX (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edBeamX as text
%        str2double(get(hObject,'String')) returns contents of edBeamX as a double
%beamx = str2double(get(handles.edBeamX, 'string'));
updateAll(handles);

function edBeamY_Callback(hObject, eventdata, handles)
% hObject    handle to edBeamY (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edBeamY as text
%        str2double(get(hObject,'String')) returns contents of edBeamY as a double
%beamy  =str2double(get(handles.edBeamY, 'string'));
updateAll(handles);


function edSDD_Callback(hObject, eventdata, handles)
% hObject    handle to edSDD (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edSDD as text
%        str2double(get(hObject,'String')) returns contents of edSDD as a double
%SDD = str2double(get(handles.edSSD, 'string'));
updateAll(handles);


function edYaw_Callback(hObject, eventdata, handles)
% hObject    handle to edYaw (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edYaw as text
%        str2double(get(hObject,'String')) returns contents of edYaw as a double
updateAll(handles);
%yaw = str2double(get(handles.edYaw, 'string'));




function edBeamxLb_Callback(hObject, eventdata, handles)
% hObject    handle to edBeamxLb (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edBeamxLb as text
%        str2double(get(hObject,'String')) returns contents of edBeamxLb as a double
updateAll(handles);



function edBeamyLb_Callback(hObject, eventdata, handles)
% hObject    handle to edBeamyLb (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edBeamyLb as text
%        str2double(get(hObject,'String')) returns contents of edBeamyLb as a double
updateAll(handles);




function edSDDLb_Callback(hObject, eventdata, handles)
% hObject    handle to edSDDLb (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edSDDLb as text
%        str2double(get(hObject,'String')) returns contents of edSDDLb as a double
updateAll(handles);


function edYawLb_Callback(hObject, eventdata, handles)
% hObject    handle to edYawLb (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edYawLb as text
%        str2double(get(hObject,'String')) returns contents of edYawLb as a double
updateAll(handles);


function edBeamxUb_Callback(hObject, eventdata, handles)
% hObject    handle to edBeamxUb (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edBeamxUb as text
%        str2double(get(hObject,'String')) returns contents of edBeamxUb as a double
updateAll(handles);


function edBeamyUb_Callback(hObject, eventdata, handles)
% hObject    handle to edBeamyUb (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edBeamyUb as text
%        str2double(get(hObject,'String')) returns contents of edBeamyUb as a double
updateAll(handles);


function edSDDUb_Callback(hObject, eventdata, handles)
% hObject    handle to edSDDUb (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edSDDUb as text
%        str2double(get(hObject,'String')) returns contents of edSDDUb as a double
updateAll(handles);




function edYawUb_Callback(hObject, eventdata, handles)
% hObject    handle to edYawUb (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edYawUb as text
%        str2double(get(hObject,'String')) returns contents of edYawUb as a double
updateAll(handles);



% --- Executes on button press in pbFit.
function pbFit_Callback(hObject, eventdata, handles)
% hObject    handle to pbFit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


qCalData=evalin('base','qCalData');
userData = qCalData.userData;
calibrant = qCalData.calibrant;
scatt = qCalData.scatt;

userData.xeng = str2double(get(handles.edEnergy, 'string'));
userData.wavelength = 12.398418746/userData.xeng;
calibrant.wavelength = userData.wavelength;

x0 = [userData.beamx userData.beamy userData.SDD/userData.pSize userData.yaw];
lb = [userData.beamxLb userData.beamyLb userData.SDDLb/userData.pSize userData.yawLb];
ub = [userData.beamxUb userData.beamyUb userData.SDDUb/userData.pSize userData.yawUb];

if ~isempty(calibrant.data)
    [xval, fval] = fminsearchcon(@(x) tiltcircle2(x,calibrant), x0, lb, ub);
    plotTiltCircle3(xval,calibrant);
    disp(sprintf("Current penalty: %.5e", fval));

    beamx = xval(1);
    beamy = xval(2);
    sdd   = xval(3) * userData.pSize;
    yaw   = xval(4);

    set(handles.edBeamX, 'String', beamx);
    set(handles.edBeamY, 'String', beamy);
    set(handles.edSDD, 'String', sdd);
    set(handles.edYaw, 'String', yaw);

    userData.beamx = beamx;
    userData.beamy = beamy;
    userData.SDD = sdd;
    userData.yaw = yaw;

    scatt = [];
    scatt.eng = userData.xeng;
    scatt.wavelength = userData.wavelength;
    scatt.BeamXY = [beamx beamy];
    scatt.yaw = yaw;
    scatt.SDD = sdd;
    scatt.SDDrel = xval(3);
    scatt.limits = userData.limits;
    scatt.offset = userData.offset;
    qCalData.scatt = scatt;
    qCalData.userData = userData;
    assignin('base','qCalData', qCalData);
else
    disp('Please pick the pionts first ....')
end


% --- Executes on button press in pbCalQCorrMaps.
function pbCalQCorrMaps_Callback(hObject, eventdata, handles)
% hObject    handle to pbCalQCorrMaps (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
updateAll(handles);

[qCalData, userData, scatt, calibrant]= getData;
% if ~isfield(scatt, 'wavelength')
%     scatt.wavelength = 12.398418746/scatt.xeng
% end
%change 8/7/2014
%[scatt.qMap, scatt.qCMap, qMax, qMin]=calQCorrMap2(scatt.BeamXY, scatt.SDDrel, [scatt.yaw 0 0], scatt.wavelength, [scatt.row scatt.col]);
[scatt.qMap, scatt.qCMap, qMax, qMin]=calQCorrMap2(scatt.BeamXY-1, scatt.SDDrel, [scatt.yaw 0 0], scatt.wavelength, [scatt.row scatt.col]);

setData(qCalData, userData, scatt, calibrant);

set(handles.edQMin, 'String', num2str(qMin, '%.5f'));
set(handles.edQMax, 'String', num2str(qMax, '%.5f'));
updateAll(handles);


function edPtsNum_Callback(hObject, eventdata, handles)
% hObject    handle to edPtsNum (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edPtsNum as text
%        str2double(get(hObject,'String')) returns contents of edPtsNum as a double
updateAll(handles);

%here here
% --- Executes during object creation, after setting all properties.
%%function edit39_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edPtsNum (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
%%if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
%%    set(hObject,'BackgroundColor','white');
%%end




% --- Executes when selected object is changed in bgDetector.
function bgDetector_SelectionChangeFcn(hObject, eventdata, handles)
% hObject    handle to the selected object in bgDetector 
% eventdata  structure with the following fields (see UIBUTTONGROUP)
%	EventName: string 'SelectionChanged' (read only)
%	OldValue: handle of the previously selected object or empty if none was selected
%	NewValue: handle of the currently selected object
% handles    structure with handles and user data (see GUIDATA)

switch get(eventdata.NewValue, 'Tag')
    
    case 'rbPilatus2m'
        setPilatus2mInitials(handles);
        updateAll(handles);
        %scatt=struct( 'row',1679, 'col', 1475);
        
    case 'rbPilatus300'
        setPilatus300Initials(handles);
        updateAll(handles);
        %scatt=struct( 'row',619, 'col', 487);
        
    case 'rbEiger9m'
        setEiger9mInitials(handles);
        updateAll(handles);
        %scatt=struct( 'row',619, 'col', 487);
        
    case 'rbPE'
        setPEInitials(handles);
        updateAll(handles);
        %scatt=struct( 'row',4318, 'col', 4320);
        
    case 'rbPilatus900k'
        %setMar300Initials(handles);
        setPilatus900kInitials(handles);
        updateAll(handles);
        %scatt=struct( 'row',2048, 'col', 2048);
        
    case 'rbOther'
        setOtherDetectorInitials(handles);
        updateAll(handles);
        %scatt=struct( 'row',619, 'col', 487);    
        
    otherwise    
        disp('Should not be here Detector Choice!');
end


function setPilatus300Initials(handles)
set(handles.edBeamX,   'String', 200.0);
set(handles.edBeamxLb, 'String', 50.0);
set(handles.edBeamxUb, 'String', 400.0);

set(handles.edBeamY,   'String', 900.0);
set(handles.edBeamyLb, 'String', 600.0);
set(handles.edBeamyUb, 'String', 2000.0);

set(handles.edSDD,   'String', 500.0);
set(handles.edSDDLb, 'String', 300.0);
set(handles.edSDDUb, 'String', 600.0);

set(handles.edYaw,   'String', 16.0);
set(handles.edYawLb, 'String', 14.0);
set(handles.edYawUb, 'String', 17.0);

set(handles.edPixelsize, 'String', 0.172);
set(handles.edRow, 'String', 619);
set(handles.edCol, 'String', 487);

set(handles.edOffset,   'String', 0);
set(handles.edLimitLb,   'String', -1);
set(handles.edLimitUb,   'String', 1000000);


function setPilatus2mInitials(handles)
set(handles.edBeamX,   'String', 750.0);
set(handles.edBeamxLb, 'String', 500.0);
set(handles.edBeamxUb, 'String', 2000.0);

set(handles.edBeamY,   'String', 180.0);
set(handles.edBeamyLb, 'String', 90.0);
set(handles.edBeamyUb, 'String', 200.0);

set(handles.edSDD,   'String', 2000.0);
set(handles.edSDDLb, 'String', 800.0);
set(handles.edSDDUb, 'String', 4000.0);

set(handles.edYaw,   'String', 0.0);
set(handles.edYawLb, 'String', 0.0);
set(handles.edYawUb, 'String', 0.0);

set(handles.edPixelsize, 'String', 0.172);
set(handles.edRow, 'String', 1679);
set(handles.edCol, 'String', 1475);

set(handles.edOffset,   'String', 0);
set(handles.edLimitLb,   'String', -1);
set(handles.edLimitUb,   'String', 1000000);

function setEiger9mInitials(handles)
set(handles.edBeamX,   'String', 1440.0);
set(handles.edBeamxLb, 'String', 500.0);
set(handles.edBeamxUb, 'String', 2000.0);

set(handles.edBeamY,   'String', 206.0);
set(handles.edBeamyLb, 'String', 90.0);
set(handles.edBeamyUb, 'String', 300.0);

set(handles.edSDD,   'String', 2000.0);
set(handles.edSDDLb, 'String', 800.0);
set(handles.edSDDUb, 'String', 4000.0);

set(handles.edYaw,   'String', 0.0);
set(handles.edYawLb, 'String', 0.0);
set(handles.edYawUb, 'String', 0.0);

set(handles.edPixelsize, 'String', 0.075);
set(handles.edRow, 'String', 3262);
set(handles.edCol, 'String', 3108);

set(handles.edOffset,   'String', 0);
set(handles.edLimitLb,   'String', -1);
set(handles.edLimitUb,   'String', 2000000);

function setPEInitials(handles)
set(handles.edBeamX,   'String', 2150.0);
set(handles.edBeamxLb, 'String', 2000.0);
set(handles.edBeamxUb, 'String', 2500.0);

set(handles.edBeamY,   'String', 130.0);
set(handles.edBeamyLb, 'String', 90.0);
set(handles.edBeamyUb, 'String', 300.0);

set(handles.edSDD,   'String', 150.0);
set(handles.edSDDLb, 'String', 80.0);
set(handles.edSDDUb, 'String', 400.0);

set(handles.edYaw,   'String', 0.5);
set(handles.edYawLb, 'String', -5.0);
set(handles.edYawUb, 'String', 5.0);

set(handles.edPixelsize, 'String', 0.100);
set(handles.edRow, 'String', 4318);
set(handles.edCol, 'String', 4320);

set(handles.edOffset,   'String', 100);
set(handles.edLimitLb,   'String', -1);
set(handles.edLimitUb,   'String', 100000);

function setMar300Initials(handles)
set(handles.edBeamX,   'String', 821.0);
set(handles.edBeamxLb, 'String', 200.0);
set(handles.edBeamxUb, 'String', 1200.0);

set(handles.edBeamY,   'String', 1682.0);
set(handles.edBeamyLb, 'String', 900.0);
set(handles.edBeamyUb, 'String', 2000.0);

set(handles.edSDD,   'String', 2200.0);
set(handles.edSDDLb, 'String', 800.0);
set(handles.edSDDUb, 'String', 3000.0);

set(handles.edYaw,   'String', 0.5);
set(handles.edYawLb, 'String', -2.0);
set(handles.edYawUb, 'String', 2.0);

set(handles.edPixelsize, 'String', 0.1465);
set(handles.edRow, 'String', 2048);
set(handles.edCol, 'String', 2048);

set(handles.edOffset,   'String', 100);
set(handles.edLimitLb,   'String', -1);
set(handles.edLimitUb,   'String', 100000);

function setPilatus900kInitials(handles)
set(handles.edBeamX,   'String', 821.0);
set(handles.edBeamxLb, 'String', 200.0);
set(handles.edBeamxUb, 'String', 1200.0);

set(handles.edBeamY,   'String', 1682.0);
set(handles.edBeamyLb, 'String', 900.0);
set(handles.edBeamyUb, 'String', 2000.0);

set(handles.edSDD,   'String', 2200.0);
set(handles.edSDDLb, 'String', 800.0);
set(handles.edSDDUb, 'String', 3000.0);

set(handles.edYaw,   'String', 0.5);
set(handles.edYawLb, 'String', -2.0);
set(handles.edYawUb, 'String', 2.0);

set(handles.edPixelsize, 'String', 0.172);
set(handles.edRow, 'String', 1043);
set(handles.edCol, 'String', 981);

set(handles.edOffset,   'String', 0);
set(handles.edLimitLb,   'String', -1);
set(handles.edLimitUb,   'String', 100000);


function setOtherDetectorInitials(handles)
set(handles.edBeamX,   'String', 500.0);
set(handles.edBeamxLb, 'String', 100.0);
set(handles.edBeamxUb, 'String', 800.0);

set(handles.edBeamY,   'String', 500.0);
set(handles.edBeamyLb, 'String', 100.0);
set(handles.edBeamyUb, 'String', 800.0);

set(handles.edSDD,   'String', 2000.0);
set(handles.edSDDLb, 'String', 200.0);
set(handles.edSDDUb, 'String', 3000.0);

set(handles.edYaw,   'String', 0.0);
set(handles.edYawLb, 'String', 0.0);
set(handles.edYawUb, 'String', 0.0);

set(handles.edPixelsize, 'String', 0.100);
set(handles.edRow, 'String', 1024);
set(handles.edCol, 'String', 1024);

set(handles.edOffset,   'String', 0);
set(handles.edLimitLb,   'String', -1);
set(handles.edLimitUb,   'String', 1000000);

function setDetectorInitials(handles, isPilatus2m)
if isPilatus2m == 1
    setPilatus2mInitials(handles);
elseif isPilatus2m == 0
    setPilatus300Initials(handles);
elseif isPilatus2m == 2
    setPEInitials(handles);
elseif isPilatus2m == 3
    %setMar300Initials(handles); 
    setPilatus900kInitials(handles);  
elseif isPilatus2m == 4
    setEiger9mInitials(handles);    
elseif isPilatus2m == -1
    setOtherDetectorInitials(handles);
end

% function setOtherInitials(handles)
%     
% end

function updateUserdata(handles)
%update parameter inputs
qCalData=evalin('base','qCalData');
userData=qCalData.userData;

userData.xeng = str2double(get(handles.edEnergy, 'string'));
userData.wavelength = 12.398418746/userData.xeng;

userData.pSize = str2double(get(handles.edPixelsize, 'string'));

userData.beamx  = str2double(get(handles.edBeamX, 'string'));
userData.beamy  = str2double(get(handles.edBeamY, 'string'));
userData.SDD    = str2double(get(handles.edSDD, 'string'));
userData.yaw = str2double(get(handles.edYaw, 'string'));

userData.SDDLb = str2double(get(handles.edSDDLb, 'string'));
userData.beamxLb = str2double(get(handles.edBeamxLb, 'string'));
userData.beamyLb = str2double(get(handles.edBeamyLb, 'string'));
userData.yawLb = str2double(get(handles.edYawLb, 'string'));

userData.beamxUb = str2double(get(handles.edBeamxUb, 'String'));
userData.beamyUb = str2double(get(handles.edBeamyUb, 'String'));
userData.SDDUb = str2double(get(handles.edSDDUb, 'String'));
userData.yawUb = str2double(get(handles.edYawUb, 'String'));

userData.offset = str2double(get(handles.edOffset,'String'));
limitUb = str2double(get(handles.edLimitUb,'String'));
limitLb = str2double(get(handles.edLimitLb,'String'));
userData.limits=[limitLb limitUb];
userData.absIntCoeff = str2double(get(handles.edAbsCoeff,'String'));

if get(handles.rbPilatus2m,'Value')
    userData.detector = 'Pilatus2M';
    userData.isPilatus2m = 1;
    userData.detRow=1679;
    userData.detCol=1475;
    userData.pSize = 0.172;
elseif get(handles.rbPilatus300,'Value')
    userData.detector = 'Pilatus300K';
    userData.isPilatus2m = 0;
    userData.detRow=619;
    userData.detCol=487;
    userData.pSize = 0.172;
elseif get(handles.rbPE,'Value')
    userData.detector = 'PEXRpad';
    userData.isPilatus2m = 2;
    userData.detRow=4318;
    userData.detCol=4320;    
    userData.pSize = 0.100;
elseif get(handles.rbPilatus900k,'Value')
    %userData.detector = 'Mar300';
    %userData.detRow=2048;
    %userData.detCol=2048;      
    %userData.pSize = 0.1465;    
    userData.detector = 'Pilatus900K';
    userData.isPilatus2m = 3;
    userData.detCol=981; 
    userData.detRow=1043;         
    userData.pSize = 0.172;    
elseif get(handles.rbEiger9m,'Value')
    userData.detector = 'Eiger9M';
    userData.isPilatus2m = 4;
    userData.detRow=3262;
    userData.detCol=3108;
    userData.pSize = 0.075;    
elseif get(handles.rbOther,'Value')
    userData.detector = 'Undefined';
    userData.isPilatus2m = -1;
    userData.detRow=str2num(get(handles.edRow, 'String'));
    userData.detCol=str2num(get(handles.edCol, 'String'));    
end

set(handles.edPixelsize, 'String',userData.pSize );
set(handles.edRow, 'String', userData.detRow);
set(handles.edCol, 'String', userData.detCol)

% userData.isPilatus2m =  get(handles.rbPilatus2m,'Value');
% if userData.isPilatus2m > 0
%     userData.detRow=1679;
%     userData.detCol=1475;
% else
%     userData.detRow=619;
%     userData.detCol=487;
% end

userData.qNum = str2double(get(handles.edPtsNum,'String'));
userData.qMin = str2double(get(handles.edQMin,'String'));
userData.qMax = str2double(get(handles.edQMax,'String'));
userData.qMax = str2double(get(handles.edQMax,'String'));

if get(handles.combLogLinear, 'Value')
    userData.dq1 = str2double(get(handles.cLogQ2,'String'));
    userData.q2 = str2double(get(handles.cLogPts,'String'));
    userData.dq2 = str2double(get(handles.cLogDQ2,'String'));
    userData.q3 = str2double(get(handles.cLogQ3,'String'));
    userData.dq3 = str2double(get(handles.cLogDQ3,'String'));    
else
    userData.dq1 = str2double(get(handles.edDQ1,'String'));
    userData.q2 = str2double(get(handles.edQ2,'String'));
    userData.dq2 = str2double(get(handles.edDQ2,'String'));
    userData.q3 = str2double(get(handles.edQ3,'String'));
    userData.dq3 = str2double(get(handles.edDQ3,'String'));
end
userData.combQs=[userData.dq1 userData.q2 userData.dq2 userData.q3 userData.dq3];  % combined for linear & Log-linear
userData.qChoice = getQChoice(handles);

% index_selected = get(hObject,'Value');
% list = get(hObject,'String');
% item_selected = list{index_selected};
% 
% scatt.absIntStandardIndex =  index_selected;
% scatt.absIntStandard = item_selected;

userData.BeamlineIndex = get(handles.Beamline,'Value');
list = get(handles.Beamline,'String');
userData.Beamline = list{userData.BeamlineIndex};

userData.absIntStandardIndex = get(handles.absIntStandard,'Value');
list = get(handles.absIntStandard,'String');
userData.absIntStandard = list{userData.absIntStandardIndex};

userData.qStandardIndex = get(handles.qStandard,'Value');
list = get(handles.qStandard,'String');
userData.qStandard = list{userData.qStandardIndex};

%average folder
% set(handles.rdAve, 'Value', 1);
% userData.avgFolderChoice = 1;
% userData.avgFolder='Ave';
set(handles.rdAveraged, 'Value', 1);
userData.avgFolderChoice = 2;
userData.avgFolder='Averaged';

userData.sectorA1 = str2num(get(handles.sectorA1,'String'));
userData.sectorA2 = str2num(get(handles.sectorA2,'String'));
userData.sectorMaskOn = get(handles.useSectorMask,'Value');
userData.sectorMaskExclude = get(handles.excludeSectorMask,'Value');

userData.arbitMaskOn = get(handles.arbitMaskOn,'Value') ;
userData.arbitMaskExclude = get(handles.arbitMaskExclude,'Value');


qCalData.userData = userData;
assignin('base','qCalData', qCalData);

function updateScatt(handles)
%qCalData=evalin('base','qCalData');
%userData=qCalData.userData;
%scatt = qCalData.scatt;
[qCalData, userData, scatt, calibrant]= getData;

scatt.limits=userData.limits;
scatt.isPilatus2m = userData.isPilatus2m;
scatt.detector =userData.detector;
scatt.row = userData.detRow;
scatt.col = userData.detCol;
scatt.offset = userData.offset;
scatt.absIntCoeff = userData.absIntCoeff;
scatt.eng = userData.xeng;
scatt.wavelength = userData.wavelength;
scatt.pSize = userData.pSize;
scatt.qNum = userData.qNum;
scatt.qMin = userData.qMin;
scatt.qMax = userData.qMax;
scatt.combQs = userData.combQs;
scatt.qChoice = userData.qChoice;
scatt.SDD = userData.SDD;
scatt.SDDrel = userData.SDD / userData.pSize;
scatt.BeamXY = [userData.beamx userData.beamy];
scatt.yaw = userData.yaw;

scatt.absIntStandardIndex = userData.absIntStandardIndex;
scatt.absIntStandard = userData.absIntStandard;
scatt.qStandardIndex = userData.qStandardIndex;
scatt.qStandard = userData.qStandard;
scatt.BeamlineIndex = userData.BeamlineIndex;
scatt.Beamline = userData.Beamline;

scatt.sectorA1 = userData.sectorA1;
scatt.sectorA2 = userData.sectorA2;
scatt.sectorMaskOn = userData.sectorMaskOn;
scatt.sectorMaskExclude = userData.sectorMaskExclude;
scatt.sectorMask = [];

scatt.arbitMaskOn = userData.arbitMaskOn;
scatt.arbitMask = [];
scatt.arbitMaskExclude = userData.arbitMaskExclude;


setData(qCalData, userData, scatt, calibrant);
%qCalData.scatt = scatt;
%assignin('base','qCalData', qCalData);



function updateAll(handles)
updateUserdata(handles);
updateScatt(handles);

function edOffset_Callback(hObject, eventdata, handles)
% hObject    handle to edOffset (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edOffset as text
%        str2double(get(hObject,'String')) returns contents of edOffset as a double
%handles.offset = str2double(get(handles.edOffset,'String'));
updateAll(handles);



function edLimitUb_Callback(hObject, eventdata, handles)
% hObject    handle to edLimitUb (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edLimitUb as text
%        str2double(get(hObject,'String')) returns contents of edLimitUb as a double
updateAll(handles);


function edLimitLb_Callback(hObject, eventdata, handles)
% hObject    handle to edLimitLb (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edLimitLb as text
%        str2double(get(hObject,'String')) returns contents of edLimitLb as a double
updateAll(handles);


% --- Executes during object creation, after setting all properties.
function edLimitLb_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edLimitLb (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes during object creation, after setting all properties.
function edOffset_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edOffset (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edQMax_Callback(hObject, eventdata, handles)
% hObject    handle to edQMax (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edQMax as text
%        str2double(get(hObject,'String')) returns contents of edQMax as a double
updateAll(handles);


function edQMin_Callback(hObject, eventdata, handles)
% hObject    handle to edQMin (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edQMin as text
%        str2double(get(hObject,'String')) returns contents of edQMin as a double
updateAll(handles);

% --- Executes during object creation, after setting all properties.
function edQMin_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edQMin (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes during object creation, after setting all properties.
function edQMax_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edQMax (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes during object deletion, before destroying properties.
function rbPilatus2m_DeleteFcn(hObject, eventdata, handles)
% hObject    handle to rbPilatus2m (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes during object creation, after setting all properties.
function rbPilatus2m_CreateFcn(hObject, eventdata, handles)
% hObject    handle to rbPilatus2m (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called


% --- Executes during object creation, after setting all properties.
function edPixelsize_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edPixelsize (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes during object creation, after setting all properties.
function edEnergy_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edEnergy (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes during object creation, after setting all properties.
function edLimitUb_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edLimitUb (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes during object creation, after setting all properties.
function edBeamX_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edBeamX (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes during object creation, after setting all properties.
function edBeamY_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edBeamY (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes during object creation, after setting all properties.
function edSDD_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edSDD (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes during object creation, after setting all properties.
function edYaw_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edYaw (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes during object creation, after setting all properties.
function edBeamxLb_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edBeamxLb (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes during object creation, after setting all properties.
function edBeamyLb_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edBeamyLb (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes during object creation, after setting all properties.
function edSDDLb_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edSDDLb (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes during object creation, after setting all properties.
function edYawLb_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edYawLb (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes during object creation, after setting all properties.
function edBeamxUb_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edBeamxUb (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes during object creation, after setting all properties.
function edBeamyUb_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edBeamyUb (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes during object creation, after setting all properties.
function edSDDUb_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edSDDUb (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes during object creation, after setting all properties.
function edYawUb_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edYawUb (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes during object creation, after setting all properties.
function edPtsNum_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edPtsNum (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes during object creation, after setting all properties.
function bgDetector_CreateFcn(hObject, eventdata, handles)
% hObject    handle to bgDetector (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called


% --- Executes during object creation, after setting all properties.
function rbPilatus300_CreateFcn(hObject, eventdata, handles)
% hObject    handle to rbPilatus300 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called


% --- Executes on button press in pbCalMaskQrange.
function pbCalMaskQrange_Callback(hObject, eventdata, handles)
% hObject    handle to pbCalMaskQrange (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
[qCalData, userData, scatt, calibrant]= getData;

setData(qCalData, userData, scatt, calibrant);

% --- Executes on button press in pbLoadMaskFile.
function pbLoadMaskFile_Callback(hObject, eventdata, handles)
% hObject    handle to pbLoadMaskFile (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
[qCalData, userData, scatt, calibrant]= getData;

[file, datadir]=uigetfile({'*.*','All format(*.*)'});
if file==0 return; end
file = [datadir, file];
mask = double(imread(file));
%mask = imgUpsideDn(double(imread(file)));
if scatt.isPilatus2m== -1
    scatt.row = floor(str2num(get(handles.edRow, 'String')));
    scatt.col = floor(str2num(get(handles.edCol, 'String')));
end
if size(mask)==[scatt.row scatt.col]
    scatt.mask = mask;
    assignin('base','qCalData', qCalData);
    disp('Mask is loaded!')
else
    disp('Choose wrong mask file, dimension does not match ....')
end  
setData(qCalData, userData, scatt, calibrant);



% --- Executes when selected object is changed in uipanel13.
function uipanel13_SelectionChangeFcn(hObject, eventdata, handles)
% hObject    handle to the selected object in uipanel13 
% eventdata  structure with the following fields (see UIBUTTONGROUP)
%	EventName: string 'SelectionChanged' (read only)
%	OldValue: handle of the previously selected object or empty if none was selected
%	NewValue: handle of the currently selected object
% handles    structure with handles and user data (see GUIDATA)
%qCalData=evalin('base', 'qCalData');

[qCalData, userData, scatt, calibrant]= getData;
%setData(qCalData, userData, scatt, calibrant);
qChoice = 2;
switch get(eventdata.NewValue, 'Tag')
    
    case 'rbPixelBasedQ'
        row = scatt.row;
        col = scatt.col;
       
        qChoice = 1;
        qNum = pixelBaseQ(row, col, scatt.BeamXY);
        
        scatt.qMin = 0;
        userData.qMin = 0;
        set(handles.edQMin, 'String', 0);
        set(handles.edPtsNum, 'String', qNum);

        
    case 'rbLinearQ'
        qChoice = 2;
        
    case 'rbCombLinearQ'
        qChoice = 3;    
        scatt.qArray = updateComLinearQ(handles);
%         qMin = userData.qMin;
%         qMax = userData.qMax; 
%         
%         if isfield(scatt, 'qArray')
%             scatt.qArray = [];
%         end
%         
%         dq1 = userData.combQs(1);
%         q2 = userData.combQs(2);
%         dq2 = userData.combQs(3);
%         q3 = userData.combQs(4);
%         dq3 = userData.combQs(5);
%         qA1=[];
%         qA2=[];
%         qA3=[];
%         qArray = [];
%         if dq1>0 
%             if q2 < 0
%                 qArray = qMin:dq1:qMax;
%                 qNum = length(qArray);
%                 set(handles.edPtsNum, 'String', qNum);
%             elseif q2 > qMin
%                 qA1=qMin:dq1:q2;
%                 
%                 if dq2 > 0
%                     if q3 < 0
%                        qA2=q2+dq2:dq2:qMax;
%                        qArray = [qA1 qA2];
%                        qNum = length(qArray);
%                        set(handles.edPtsNum, 'String', qNum);
%                     elseif q3 > q2
%                         qA2=q2+dq2:dq2:q3;
%                         
%                         if dq3 > 0
%                             if q3 < qMax
%                                 qA3 = q3+dq3 : dq3:qMax;
%                                 qArray = [qA1 qA2 qA3];
%                                 qNum = length(qArray);
%                                 set(handles.edPtsNum, 'String', qNum);
%                             else
%                                 disp('Setting of Q3 is wrong! Please check!')
%                             end
%                         else
%                             disp('Setting of dQ3 is wrong! Please check!')
%                         end
%                     else
%                         disp('Setting of Q3 is wrong! Please check!')
%                     end
%                 else
%                     disp('Setting of dQ2 is wrong! Please check!')
%                 end
%             else
%                 disp('Setting of Q2 is wrong! Please check!')
%             end
%         else
%             disp('Setting of dQ1 is wrong! Please check!')
%         end
%         
%         scatt.qArray = qArray;
%         
    case 'rbLogQ'
        qChoice = 4;
        
    case 'combLogLinear'
        qChoice = 5;    
        scatt.qArray = updateComLogLinearQ(handles);        
end

userData.qChoice = qChoice;
setData(qCalData, userData, scatt, calibrant);
updateAll(handles);

% % --- dynamically update q array ---
% function updateQArray(handles)
% [qCalData, userData, scatt, calibrant]= getData;
% 
% qChoice = 2;
% if get(handles.rbPixelBasedQ, 'Value')    
%     %case 'rbPixelBasedQ'
%     row = scatt.row;
%     col = scatt.col;
% 
%     qChoice = 1;
%     qNum = pixelBaseQ(row, col, scatt.BeamXY);
% 
%     scatt.qMin = 0;
%     userData.qMin = 0;
%     set(handles.edQMin, 'String', 0);
%     set(handles.edPtsNum, 'String', qNum);
% 
% elseif get(handles.rbLinearQ, 'Value')
%     %case 'rbLinearQ'
%     qChoice = 2;
% 
% case 'rbCombLinearQ'
%     qChoice = 3;    
%     scatt.qArray = updateComLinearQ(handles);
% 
% case 'rbLogQ'
%     qChoice = 4;
% 
% case 'combLogLinear'
%     qChoice = 5;    
%     scatt.qArray = updateComLogLinearQ(handles);        
% end
% 
% 
% scatt.qArray = qArray;
% setData(qCalData, userData, scatt, calibrant);   
% updateAll(handles);
% end

% --- update combined linear q array ---
function qArray = updateComLinearQ(handles)
[qCalData, userData, scatt, calibrant]= getData;

qMin = userData.qMin;
qMax = userData.qMax; 

% if isfield(scatt, 'qArray')
%     scatt.qArray = [];
% end

dq1 = userData.combQs(1);
q2 = userData.combQs(2);
dq2 = userData.combQs(3);
q3 = userData.combQs(4);
dq3 = userData.combQs(5);
qA1=[];
qA2=[];
qA3=[];
qArray = [];
if dq1>0 
    if q2 < 0
        qArray = qMin:dq1:qMax;
        qNum = length(qArray);
        set(handles.edPtsNum, 'String', qNum);
    elseif q2 > qMin
        qA1=qMin:dq1:q2;

        if dq2 > 0
            if q3 < 0
               qA2=q2+dq2:dq2:qMax;
               qArray = [qA1 qA2];
               qNum = length(qArray);
               set(handles.edPtsNum, 'String', qNum);
            elseif q3 > q2
                qA2=q2+dq2:dq2:q3;

                if dq3 > 0
                    if q3 < qMax
                        qA3 = q3+dq3 : dq3:qMax;
                        qArray = [qA1 qA2 qA3];
                        qNum = length(qArray);
                        set(handles.edPtsNum, 'String', qNum);
                    else
                        disp('Setting of Q3 is wrong! Please check!')
                    end
                else
                    disp('Setting of dQ3 is wrong! Please check!')
                end
            else
                disp('Setting of Q3 is wrong! Please check!')
            end
        else
            disp('Setting of dQ2 is wrong! Please check!')
        end
    else
        disp('Setting of Q2 is wrong! Please check!')
    end
else
    disp('Setting of dQ1 is wrong! Please check!')
end


function qArray = updateComLogLinearQ(handles)
[qCalData, userData, scatt, calibrant]= getData;

qMin = userData.qMin;
qMax = userData.qMax; 

% if isfield(scatt, 'qArray')
%     scatt.qArray = [];
% end
% 
% cLogQ2  = userData.cLogLinear(1);
% cLogPts = userData.cLogLinear(2);
% cLogDQ2 = userData.cLogLinear(3);
% cLogQ3  = userData.cLogLinear(4);
% cLogDQ3 = userData.cLogLinear(5);


cLogQ2  = str2num(get(handles.cLogQ2, 'String'));
cLogPts = str2num(get(handles.cLogPts, 'String'));
cLogDQ2 = str2num(get(handles.cLogDQ2, 'String'));
cLogQ3  = str2num(get(handles.cLogQ3, 'String'));
cLogDQ3 = str2num(get(handles.cLogDQ3, 'String'));


% scatt.combQs(1) = cLogQ2;
% scatt.combQs(2) = cLogPts;
% scatt.combQs(3) = cLogDQ2;
% scatt.combQs(4) = cLogQ3;
% scatt.combQs(5) = cLogDQ3;

q2 = cLogQ2;
logPts = cLogPts;

dq2 = cLogDQ2;
q3= cLogQ3;
dq3= cLogDQ3;

qA1=[];
qA2=[];
qA3=[];
qArray = [];
if q2> qMin & q2 < qMax
    if logPts >= 1
        if qMin == 0
            disp('qMin Cannot be Zero for Log scale! Please check!');
            return;
        else
        qA1 = logspace(log10(qMin), log10(q2), round(logPts));
        end
    else
        disp('Setting of Log Pts is wrong! Please check!');
        return;
    end    

else
    disp('Setting of Log Q2 is wrong! Please check!');
    return;
end
    
if dq2 > 0
    if q3 > q2 
        if q3 > qMax
            q3 = qMax;  
            qA2 = q2+dq2:dq2:q3;            
        else
            if dq3 > 0
                qA2 = q2+dq2:dq2:q3;
                qA3 = q3+dq3:dq3:qMax;
                qA2 = [qA2 qA3];
            else
                disp('Setting of Linear dQ3 is wrong! Please check!');
                return;
            end
        end
    else
        disp('Setting of Linear Q3 is wrong! Please check!');
        return;        
    end
else
    disp('Setting of Linear dQ2 is wrong! Please check!');
    return;
end
qArray = [qA1 qA2];
qNum = length(qArray);
if get(handles.combLogLinear, 'Value')
    set(handles.edPtsNum, 'String', qNum);
end
% setData(qCalData, userData, scatt, calibrant);
%updateAll(handles);


% scatt.qArray = qArray;
% setData(qCalData, userData, scatt, calibrant);   
% updateAll(handles);

% --- Executes on button press in pbCalQIndxMap.
function pbCalQIndxMap_Callback(hObject, eventdata, handles)
% hObject    handle to pbCalQIndxMap (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
updateAll(handles);
updateQarray(handles);
[qCalData, userData, scatt, calibrant]= getData;

if isfield(scatt,'qMap') & isfield(scatt, 'qArray') & ~isempty(scatt.qMap) & ~isempty(scatt.qArray)
    disp('Calculating Q index map ...')
    scatt.qRMap = calQRMap(scatt.qMap, scatt.qArray);
    disp('Q index map done!')
else
    disp('Not Ready! Calculate Q and Corr Maps first ....')
end
setData(qCalData, userData, scatt, calibrant);



function qChoice=getQChoice(handles)
qChoice = 0;
if get(handles.rbPixelBasedQ, 'Value')
    qChoice = 1;
elseif get(handles.rbLinearQ, 'Value')
    qChoice = 2;
elseif get(handles.rbCombLinearQ, 'Value')
    %get(handles.rbLogQ, 'Value')
    qChoice = 3;
    
elseif get(handles.rbLogQ, 'Value')
    qChoice = 4;
    
elseif get(handles.combLogLinear, 'Value')
    %get(handles.rbLogQ, 'Value')
    qChoice = 5;   
    
end

function updateQarray(handles)
[qCalData, userData, scatt, calibrant]= getData;

qMin = userData.qMin;
qMax = userData.qMax;
qNum = userData.qNum;
isLog = 0;

row = scatt.row;
col = scatt.col;
SDDrel = scatt.SDDrel;
waveln = scatt.wavelength;

qChoice=userData.qChoice;

if qChoice == 1;
    qNum = pixelBaseQ(row, col, scatt.BeamXY);
    set(handles.edPtsNum, 'String', qNum);
    isLog = 0;
    qMin = 0;
%     scatt.qArray = genQArray(qMin, qMax, qNum, isLog );
    scatt.qArray = linspace(qMin, qMax, qNum);
    scatt.qMode = 'Pixel';
    
elseif qChoice == 2
    isLog = 0;
%     scatt.qArray = genQArray(qMin, qMax, qNum, isLog );
    scatt.qArray = linspace(qMin, qMax, qNum);
    scatt.qMode = 'Linear';
    
elseif qChoice == 3
    isLog = 0;
    %scatt.qArray = genQArray(qMin, qMax, qNum, isLog );
    
    %already bult in 
    % function uipanel13_SelectionChangeFcn(hObject, eventdata, handles)
    scatt.qArray = updateComLinearQ(handles);
    scatt.qNum = length(scatt.qArray);
    set(handles.edPtsNum, 'String', scatt.qNum);
    scatt.qMode = 'CombLinear';
    
elseif qChoice == 4
    isLog = 1;
    if qMin <= 0
        qMin = 0.001; 
        set(handles.edQMin, 'String', qMin);
        disp('Setting of qMin may be wrong! Please check!')
    end
    %scatt.qArray = genQArray(qMin, qMax, qNum, isLog );
      
    scatt.qArray = logspace(log10(qMin), log10(qMax), qNum);
    scatt.qMode = 'Log'; 

elseif qChoice == 5
    isLog = 1;
    %scatt.qArray = genQArray(qMin, qMax, qNum, isLog );
    
    %already bult in 
    % function uipanel13_SelectionChangeFcn(hObject, eventdata, handles)
    scatt.qArray = updateComLogLinearQ(handles);
    scatt.qNum = length(scatt.qArray);
    set(handles.edPtsNum, 'String', scatt.qNum);    
    scatt.qMode = 'CombLogLinear';
        
else
    
    disp('Wrong! should not be here ....')
end
setData(qCalData, userData, scatt, calibrant);



function [qCalData, userData, scatt, calibrant]= getData
qCalData = evalin('base', 'qCalData');
userData = qCalData.userData;
scatt = qCalData.scatt;
calibrant = qCalData.calibrant;

function setData(qCalData, userData, scatt, calibrant)
qCalData.userData = userData;
qCalData.scatt = scatt;
qCalData.calibrant = calibrant;
assignin('base', 'qCalData', qCalData);



%
%[qCalData, userData, scatt, calibrant]= getData;
%setData(qCalData, userData, scatt, calibrant);


% --- Executes on button press in pushbutton16.
function pushbutton16_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton16 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
%disp('here')
[qCalData, userData, scatt, calibrant]= getData;

if ~isempty(scatt.qMap)
    figure('NumberTitle','off' ,'Name','q Map'); imagesc(scatt.qMap);
end

if ~isempty(scatt.qCMap)
    figure('NumberTitle','off' ,'Name','Correction Map'); imagesc(scatt.qCMap);
end

if ~isempty(scatt.qRMap)
    figure('NumberTitle','off' ,'Name','q-index Map'); imagesc(scatt.qRMap);
end

if ~isempty(scatt.mask)
    figure('NumberTitle','off' ,'Name','mask'); imagesc(scatt.mask);
end

% --- Executes on button press in pushbutton17.
function pushbutton17_Callback(hObject, eventdata, handles)
% Process individual image, also It normalization
% hObject    handle to pushbutton17 (see GCBO)
% handles    structure with handles and user data (see GUIDATA)
[qCalData, userData, scatt, calibrant]= getData;
[file, pathname]=uigetfile({'*.tif;*.h5';'*.ccd'});
if file==0 return; end

% outputFileName=[file(1:end-3) 'dat'];
% 
% cfile = fullfile(datadir, file);
% sImg = imgUpsideDn(double(imread(cfile)));

[~, fn9, ext9] = fileparts(file);
%outputFileName=[file(1:end-3) 'dat'];
outputFileName=[fn9 '.dat'];

cfile = fullfile(pathname, file);
disp(cfile);

if strncmpi(ext9, '.tif', 4)      % tif file
    sImg = imgUpsideDn(double(imread(cfile)));
elseif strcmpi(ext9, '.h5')   % hdf5 data
    if contains(fn9, 'master')
        fprintf('%s is a master file, will not be proecssed!\n',cfile);       
    else
        if isfield(scatt, 'Beamline')
            beamline = scatt.Beamline;
        else
            beamline = '12-ID-B';
        end
        h5data = readHdf5(cfile, beamline);
        sImg = imgUpsideDn(double(h5data.data'));
    end
else
    fprintf('%s is not supported and proecssed!\n',cfile);
end    

figure('NumberTitle','off','Name',cfile);imagesc(sImg, [0 500]);

ctmp = strsplit(pathname, filesep);
%tmpParentDir = fullfile(filesep, ctmp{1:end-2});
if isempty(ctmp{1})
    tmpParentDir = fullfile(filesep,ctmp{2:end-2});
else
    tmpParentDir = fullfile(ctmp{1:end-2});
end
if strcmp(pathname(1:2), '\\')
    tmpParentDir = fullfile(filesep,tmpParentDir);
end
logPath = fullfile(tmpParentDir,'Log',filesep);
%logPath = fullfile(pathname(1:end-5),'Log',filesep);
newAvgFolder = fullfile(pathname,'Avg',filesep); 
logFileName = ['L' fn9(2:end) '.meta'];
fullLogFile = fullfile(logPath, logFileName);


[phd, Io, eng, expt, metaInfo] = parseMetafile2(fullLogFile, 0);    
if ~exist(fullLogFile, 'file')
    disp(sprintf('LogFile: %s not found, 1D data will not be normalized by It.', fullLogFile));
end

isReady = checkReadiness;
if isfield(scatt, 'absIntCoeff')
    absIntCoeff = scatt.absIntCoeff;
else
    absIntCoeff = 1.0;
end

% apply sector mask
% if scatt.sectorMaskOn && ~isempty(scatt.sectorMask)
%     mask = scatt.mask & scatt.sectorMask;
%     figure(3003); imagesc(mask); set(gca, 'YDir','normal');
% else
%     mask = scatt.mask;
% end

% apply sector mask
if scatt.sectorMaskOn && ~isempty(scatt.sectorMask)
    mask = scatt.mask & scatt.sectorMask;
    if scatt.arbitMaskOn && ~isempty(scatt.arbitMask)
        mask = mask & scatt.arbitMask;
    end
    figure(3003); imagesc(mask); set(gca, 'YDir','normal');
else
    mask = scatt.mask;
    if scatt.arbitMaskOn && ~isempty(scatt.arbitMask)
        mask = mask & scatt.arbitMask;
    end
    figure(3003); imagesc(mask); set(gca, 'YDir','normal');
end

if isReady
    tic
    %data = circavgnew2(sImg, scatt.mask, scatt.qCMap, scatt.qRMap, scatt.qArray, scatt.offset, scatt.limits);
    data = circavgnew2(sImg, mask, scatt.qCMap, scatt.qRMap, scatt.qArray, scatt.offset, scatt.limits);
    figure('NumberTitle','off','Name',['1d data for: ' cfile]); loglog(data(:,1),data(:,2)*absIntCoeff/phd, '-r');
    
    %[filename, pathname] = uiputfile('test.dat', 'Save file name(*.dat)');
    [filename, pathname] = uiputfile(outputFileName, 'Save file name(*.dat)');
    
    if filename==0
        return;
    else
        outputFile = [pathname filename];    
        errmsg ='';
        [fid5, errmsg]=fopen(outputFile, 'w');
        if fid5 > 0
            fprintf(fid5, '# Image Filename:  %s\n', file);
            fprintf(fid5, '# Original Folder: %s\n', metaInfo.datadir);
            fprintf(fid5, '# Data Collection Time:  %s\n', metaInfo.expttime);
            fprintf(fid5, '# Data Processing Time:  %s\n', datetime());
            fprintf(fid5, '# X-ray Energy (KeV):  %.4f\n', eng);
            fprintf(fid5, '# Beam Center :  %.4f, %.4f\n', scatt.BeamXY(1), scatt.BeamXY(2) );
            fprintf(fid5, '# Sample-to-Detector Distance (SDD) (mm):  %.3f\n', scatt.SDD);
            fprintf(fid5, '# Detector Pixel Size (mm):  %.4f\n', scatt.pSize);
            fprintf(fid5, '# X-ray scattering Intensities are normallized by transmitted flux(It)\n');
            fprintf(fid5, '#    and scaled to absolute values if Absolute Intensity Coefficient provided.\n');
            fprintf(fid5, '# Photodiode Value (It):  %.4f\n', phd);
            fprintf(fid5, '# I0 Value:  %.4f\n', Io);
            fprintf(fid5, '# Exposure Time (s):  %.4f\n', expt);
            fprintf(fid5, '# Scaled by Absolute Intensity Coefficient: %.4f\n', absIntCoeff);
            fprintf(fid5, '# \n');
            fprintf(fid5, '# q_(1/A)    I(q)_(cm^-1)     I(q)_SE\n');
            
            fclose(fid5);
            %disp(sprintf('%s was processed.\n',cfile));            
            if phd<=0
                phd = 1.0;
                %disp(sprintf('Negative phd value!! %s was processed, but not normalized.\n', cfile));
                disp(sprintf('Negative phd value!! %s was processed, but not normalized.', cfile));
            else
                %disp(sprintf('%s was processed.\n', cfile));
                disp(sprintf('%s was processed.', cfile));
            end
            %dlmwrite(outputFile,[data(:,1) data(:,2:3)/phd],'delimiter','\t','precision','%.6f');
            dlmwrite(outputFile,[data(:,1) data(:,2)*absIntCoeff/phd data(:,3)./sqrt(data(:,5))*absIntCoeff/phd data(:,4:6)], '-append', 'delimiter','\t','precision','%.6e');
            %dlmwrite(outputFile,[data(:,1) data(:,2:3)*absIntCoeff/phd data(:,4:6)], '-append', 'delimiter','\t','precision','%.6e');
        else
            disp(sprintf('Error: file: %s cannot open! %s.\n', outputFile, errmsg));
            %disp(errmsg);
        end   
    end  
    
    
%     if filename==0
%         return;
%     else
%         outputFile = [pathname filename];
%         dlmwrite(outputFile,[data(:,1) data(:,2:3)*absIntCoeff/phd],'delimiter','\t','precision','%.6f');
%     end
    toc
end

function isReady = checkReadiness
isReady = 1;
[qCalData, userData, scatt, calibrant]= getData;

if ~isfield(scatt, 'mask') || isempty(scatt.mask)
    disp('Mask is not ready! ...')
    isReady = 0;
    return;
end

if ~isfield(scatt, 'qCMap') || isempty(scatt.qCMap)
    isReady = 0;
    disp('qCMap is not ready! ...')
    return;
end

if ~isfield(scatt, 'qRMap') || isempty(scatt.qRMap)
    isReady = 0;
    disp('qRMap is not ready! ...')
    return;
end

if ~isfield(scatt, 'qArray') || isempty(scatt.qArray)
    isReady = 0;
    disp('qArray is not ready! ...')
    return;
end

if ~isfield(scatt, 'offset') || isempty(scatt.offset)
    isReady = 0
    disp('offset is not available! ...')
    return;
end

if ~isfield(scatt, 'limits') || isempty(scatt.limits)
    isReady = 0;
    disp('limits is not available! ...')
    return;
end


% --- Executes on button press in pushbutton18.
function pushbutton18_Callback(hObject, eventdata, handles)
% Process Multiple Image Function
%
% hObject    handle to pushbutton18 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

[qCalData, userData, scatt, calibrant]= getData;


[filenames, pathname, filterindex] = uigetfile( ...
{  '*.tif;*.h5','SAXS Image File (*.dat)'; ...
   '*.*',  'All Files (*.*)'}, ...
   'Pick a file', ...
   'MultiSelect', 'on');


% filenames=0 if cancelled
if isnumeric(filenames) 
    fprintf('No files selected!\n');
    return;
end
if ischar(filenames)   % only one file was selected.
    filenames={filenames};
end

% % sector mask
% if scatt.sectorMaskOn && ~isempty(scatt.sectorMask)
%     mask = scatt.mask & scatt.sectorMask;
%     figure(3003); imagesc(mask); set(gca, 'YDir','normal');
% else
%     mask = scatt.mask;
% end

% apply sector/arbitary mask
if isfield(scatt, 'sectorMaskOn') && scatt.sectorMaskOn && ~isempty(scatt.sectorMask)
    mask = scatt.mask & scatt.sectorMask;
    if scatt.arbitMaskOn && ~isempty(scatt.arbitMask)
        mask = mask & scatt.arbitMask;
    end
    figure(3003); imagesc(mask); set(gca, 'YDir','normal');
else
    mask = scatt.mask;
    if isfield(scatt, 'arbitMaskOn') && scatt.arbitMaskOn && ~isempty(scatt.arbitMask)
        mask = mask & scatt.arbitMask;
    end
    figure(3003); imagesc(mask); set(gca, 'YDir','normal');
end

for nfile = filenames
    tic
    h5data=[];
    file = nfile{1,1};
    [~, fn9, ext9] = fileparts(file);
    %outputFileName=[file(1:end-3) 'dat'];
    outputFileName=[fn9 '.dat'];

    cfile = fullfile(pathname, file);
    disp(cfile)
    if strncmpi(ext9, '.tif', 4)      % tif file
        sImg = imgUpsideDn(double(imread(cfile)));
    elseif strcmpi(ext9, '.h5')   % hdf5 data
        if contains(fn9, 'master')
            fprintf('%s is a master file, will not be proecssed!\n',cfile);
            continue;            
        else
            %h5data = readHdf5(cfile);
            if isfield(scatt, 'Beamline')
                beamline = scatt.Beamline;
            else
                beamline = '12-ID-B';
            end
            h5data = readHdf5(cfile, beamline);
            sImg = imgUpsideDn(double(h5data.data'));
        end
    else
        fprintf('%s is not supported and proecssed!\n',cfile);
        continue;
    end    
    %figure('NumberTitle','off','Name',cfile);imagesc(sImg);
    ctmp = strsplit(pathname, filesep);
    %tmpParentDir = fullfile(filesep, ctmp{1:end-2});
    if isempty(ctmp{1})
        tmpParentDir = fullfile(filesep,ctmp{2:end-2});
    else
        tmpParentDir = fullfile(ctmp{1:end-2});
    end
    if strcmp(pathname(1:2), '\\')
        tmpParentDir = fullfile(filesep,tmpParentDir);
    end
    logPath = fullfile(tmpParentDir,'Log',filesep);
    %logPath = fullfile(pathname(1:end-5),'Log',filesep);
    if userData.avgFolderChoice < 3
        newAvgFolder = fullfile(pathname,userData.avgFolder,filesep); 
    else
        newAvgFolder = fullfile(userData.avgFolder,filesep);
    end
    
    logFileName = ['L' fn9(2:end) '.meta'];
    fullLogFile = fullfile(logPath, logFileName);
    
    if ~exist(newAvgFolder)
        mkfolder = ['mkdir ' newAvgFolder];
        dos(mkfolder);
    end
    
    %[phd, Io, eng] = parseMetafile(fullLogFile);
    %[phd, Io, eng] = parseMetafile(fullLogFile, 0); % comment out 10/18/2017
    if exist(fullLogFile,'file')==1
        [phd, Io, eng, expt, metaInfo] = parseMetafile2(fullLogFile, -1); 
        expttime = metaInfo.expttime;
    else
        expt = -1;
        expttime = -1;
        eng = scatt.eng;
        phd = -1;
        Io = -1;
        try
            if ~isempty(h5data)
                phd = h5data.It_flux; % read from h5 data
                Io = h5data.Io_flux;
                %phd = h5data.It_phd; % read from h5 data
                %Io = h5data.IC1_phd;
                expttime = h5data.Sample_Time;
                expt=h5data.ExposureTime;
                eng = h5data.monoE;
            else
                disp(sprintf('LogFile: %s and Header Info not found,  1D data will not be normalized by It.', fullLogFile));
                %expt = -1;
                %expttime = -1;
                %eng = scatt.eng;
                %phd = -1;
                %Io = -1;
            end
        catch
            disp('Parameter missed!');
        end
    end
%     if ~exist(fullLogFile, 'file')
%         disp(sprintf('LogFile: %s not found, 1D data will not be normalized by It.', fullLogFile));
%     end
    
    isReady = checkReadiness;
    if isfield(scatt, 'absIntCoeff')
        absIntCoeff = scatt.absIntCoeff;
    else
        absIntCoeff = 1.0;
    end
    toc
    tic
    if isReady
        data = circavgnew2c(sImg, mask, scatt.qCMap, scatt.qRMap, scatt.qArray, scatt.offset, scatt.limits);
        fprintf('circ v2c\n');
            % this is the suggested version
            % data(1): Q; 
            % data(2): I(Q)--averaged from photons, corrected by solid angle, scaled by 1e-6
            % data(3): standard error of the mean of I(Q) among pixels in a given Q-bin
            % data(4): error estimated by: sqrt(I(Q)/N), corrected by solid angle, scaled by 1e-6
            % data(5): total pixel number in a given Q-bin
            % data(6): averaged photon number without solid angle correction , scaled by 1e-6 


        %data = circavgnew2(sImg, mask, scatt.qCMap, scatt.qRMap, scatt.qArray, scatt.offset, scatt.limits);
        %fprintf('circ v2b\n');
            % data(1): Q; 
            % data(2): I(Q)--averaged from photons
            % data(3): error estimated by: sqrt(I(Q)/N) among pixels in a given Q-bin
            % data(4): Average of I(Q)^2 values for a given Q-bin
            % data(5): total pixel number in a given Q-bin
            % data(6): solid angle correction value

        %data = circavgnew3(sImg, mask, scatt.qCMap, scatt.qRMap, scatt.qArray, scatt.offset, scatt.limits);
        %fprintf('circ v3\n');
            % data(1): Q; 
            % data(2): I(Q)--averaged from photons
            % data(3): error estimated by: sqrt(I(Q)/N) among pixels in a given Q-bin
            % data(4): Zero
            % data(5): total pixel number in a given Q-bin
            % data(6): solid angle correction value      
        %data = circavgnew2(sImg, scatt.mask, scatt.qCMap, scatt.qRMap, scatt.qArray, scatt.offset, scatt.limits);
        %figure('NumberTitle','off','Name',['1d data for: ' cfile]); loglog(data(:,1),data(:,2), '-r');
        
        %outputFile = [pathname outputFileName];
        outputFile = [newAvgFolder outputFileName];
        errmsg ='';
        [fid5, errmsg]=fopen(outputFile, 'w');
        if fid5 > 0            
            fprintf(fid5, '# Image Filename:  %s\n', file);
            fprintf(fid5, '# Original Folder: %s\n', pathname);
            %fprintf(fid5, '# Data Collection Time:  %s\n', expttime);
            fprintf(fid5, '# Data Processing Time:  %s\n', datetime());
            fprintf(fid5, '# X-ray Energy (KeV):  %.4f\n', eng);
            fprintf(fid5, '# X-ray Wavelength (A):  %.4f\n', scatt.wavelength);
            fprintf(fid5, '# X-ray Detector: %s\n', scatt.detector);
            fprintf(fid5, '# Detector Pixel Size (mm):  %.4f\n', scatt.pSize); 
            fprintf(fid5, '# Detector Q calibrant: %s\n', scatt.qStandard);
            fprintf(fid5, '# Beam Center :  %.4f, %.4f\n', scatt.BeamXY(1), scatt.BeamXY(2) );
            fprintf(fid5, '# Sample-to-Detector Distance (SDD) (mm):  %.3f\n', scatt.SDD);                                   
            fprintf(fid5, '# Photodiode Value (It):  %.4f\n', phd);
            fprintf(fid5, '# I0 Value:  %.4f\n', Io);
            fprintf(fid5, '# Exposure Time (s):  %.4f\n', expt);
            fprintf(fid5, '# X-ray scattering Intensities are normallized by transmitted flux(i.e. Photodiode Value)\n');
            fprintf(fid5, '# Absolute Intensity was calibrated using: %s.\n', scatt.absIntStandard);
            fprintf(fid5, '# Scaled by Absolute Intensity Coefficient: %.4f\n', absIntCoeff);
            fprintf(fid5, '# \n');
            fprintf(fid5, '# q_(1/A)    I(q)_(cm^-1)     I(q)_SE\n');
            
            fclose(fid5);
            %disp(sprintf('%s was processed.\n',cfile));            
            if phd<=0
                phd = 1.0;
                %disp(sprintf('Negative phd value!! %s was processed, but not normalized.\n', cfile));
                disp(sprintf('Negative phd value!! %s was processed, but not normalized.', cfile));
            else
                %disp(sprintf('%s was processed.\n', cfile));
                disp(sprintf('%s was processed.', cfile));
            end
            %dlmwrite(outputFile,[data(:,1) data(:,2:3)/phd],'delimiter','\t','precision','%.6f');
            % data(1): Q; 
            % data(2): I(Q)--averaged from photons, corrected by solid angle, scaled by 1e-6
            % data(3): standard deviation of I(Q) among pixels in a given Q-bin
            % data(4): total photon number solid angle corrected in a given Q-bin
            % data(5): total pixel number in a given Q-bin
            % data(6): averaged photon number without solid angle correction

            dlmwrite(outputFile,[data(:,1) data(:,2)*absIntCoeff/phd data(:,3)*absIntCoeff/phd data(:,5)], '-append', 'delimiter','\t','precision','%.6e');
            %dlmwrite(outputFile,[data(:,1) data(:,2) data(:,3) data(:,4) data(:,5) data(:,6)], '-append', 'delimiter','\t','precision','%.6e');
            %dlmwrite(outputFile,[data(:,1) data(:,2)*absIntCoeff/phd data(:,3)./sqrt(data(:,5))*absIntCoeff/phd data(:,5)], '-append', 'delimiter','\t','precision','%.6e');            
            
            fprintf('cof=%.5E.\n', absIntCoeff/phd);
            %dlmwrite(outputFile,[data(:,1) data(:,2)*absIntCoeff/phd data(:,3)*absIntCoeff/phd data(:,5)], '-append', 'delimiter','\t','precision','%.6e');
        else
            disp(sprintf('Error: file: %s cannot open! %s.\n', outputFile, errmsg));
            %disp(errmsg);
        end   

    end
    toc
    disp(' ');
end


% --- Executes during object deletion, before destroying properties.
function edBeamyLb_DeleteFcn(hObject, eventdata, handles)
% hObject    handle to edBeamyLb (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)



function edit50_Callback(hObject, eventdata, handles)
% hObject    handle to edit50 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit50 as text
%        str2double(get(hObject,'String')) returns contents of edit50 as a double


% --- Executes during object creation, after setting all properties.
function edit50_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit50 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit51_Callback(hObject, eventdata, handles)
% hObject    handle to edit51 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit51 as text
%        str2double(get(hObject,'String')) returns contents of edit51 as a double


% --- Executes during object creation, after setting all properties.
function edit51_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit51 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit52_Callback(hObject, eventdata, handles)
% hObject    handle to edit52 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit52 as text
%        str2double(get(hObject,'String')) returns contents of edit52 as a double


% --- Executes during object creation, after setting all properties.
function edit52_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit52 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edDQ1_Callback(hObject, eventdata, handles)
% hObject    handle to edDQ1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edDQ1 as text
%        str2double(get(hObject,'String')) returns contents of edDQ1 as a double


% --- Executes during object creation, after setting all properties.
function edDQ1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edDQ1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edQ2_Callback(hObject, eventdata, handles)
% hObject    handle to edQ2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edQ2 as text
%        str2double(get(hObject,'String')) returns contents of edQ2 as a double


% --- Executes during object creation, after setting all properties.
function edQ2_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edQ2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edDQ2_Callback(hObject, eventdata, handles)
% hObject    handle to edDQ2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edDQ2 as text
%        str2double(get(hObject,'String')) returns contents of edDQ2 as a double


% --- Executes during object creation, after setting all properties.
function edDQ2_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edDQ2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edQ3_Callback(hObject, eventdata, handles)
% hObject    handle to edQ3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edQ3 as text
%        str2double(get(hObject,'String')) returns contents of edQ3 as a double


% --- Executes during object creation, after setting all properties.
function edQ3_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edQ3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edDQ3_Callback(hObject, eventdata, handles)
% hObject    handle to edDQ3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edDQ3 as text
%        str2double(get(hObject,'String')) returns contents of edDQ3 as a double


% --- Executes during object creation, after setting all properties.
function edDQ3_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edDQ3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in pbCalExistQIndxMap.
function pbCalExistQIndxMap_Callback(hObject, eventdata, handles)
% hObject    handle to pbCalExistQIndxMap (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% updateQarray(handles)
 [qCalData, userData, scatt, calibrant]= getData;

if isfield(scatt,'qMap') & isfield(scatt, 'qArray') & ~isempty(scatt.qMap) & ~isempty(scatt.qArray)
    scatt.qRMap = calQRMap(scatt.qMap, scatt.qArray);
    disp('Q index map done!')
else
    disp('Not Ready! Calculate Q and Corr Maps first ....')
end
setData(qCalData, userData, scatt, calibrant);


function edRow_Callback(hObject, eventdata, handles)
% hObject    handle to edRow (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edRow as text
%        str2double(get(hObject,'String')) returns contents of edRow as a double
 [qCalData, userData, scatt, calibrant]= getData;
 scatt.row = floor(str2num(get(handles.edRow, 'String')));
 setData(qCalData, userData, scatt, calibrant);
 

% --- Executes during object creation, after setting all properties.
function edRow_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edRow (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edCol_Callback(hObject, eventdata, handles)
% hObject    handle to edCol (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edCol as text
%        str2double(get(hObject,'String')) returns contents of edCol as a double
 [qCalData, userData, scatt, calibrant]= getData;
 scatt.col = floor(str2num(get(handles.edCol, 'String')));
 setData(qCalData, userData, scatt, calibrant);

% --- Executes during object creation, after setting all properties.
function edCol_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edCol (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function cLogQ2_Callback(hObject, eventdata, handles)
% hObject    handle to cLogQ2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of cLogQ2 as text
%        str2double(get(hObject,'String')) returns contents of cLogQ2 as a double


% --- Executes during object creation, after setting all properties.
function cLogQ2_CreateFcn(hObject, eventdata, handles)
% hObject    handle to cLogQ2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function cLogPts_Callback(hObject, eventdata, handles)
% hObject    handle to cLogPts (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of cLogPts as text
%        str2double(get(hObject,'String')) returns contents of cLogPts as a double


% --- Executes during object creation, after setting all properties.
function cLogPts_CreateFcn(hObject, eventdata, handles)
% hObject    handle to cLogPts (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function cLogDQ2_Callback(hObject, eventdata, handles)
% hObject    handle to cLogDQ2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of cLogDQ2 as text
%        str2double(get(hObject,'String')) returns contents of cLogDQ2 as a double


% --- Executes during object creation, after setting all properties.
function cLogDQ2_CreateFcn(hObject, eventdata, handles)
% hObject    handle to cLogDQ2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function cLogQ3_Callback(hObject, eventdata, handles)
% hObject    handle to cLogQ3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of cLogQ3 as text
%        str2double(get(hObject,'String')) returns contents of cLogQ3 as a double


% --- Executes during object creation, after setting all properties.
function cLogQ3_CreateFcn(hObject, eventdata, handles)
% hObject    handle to cLogQ3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function cLogDQ3_Callback(hObject, eventdata, handles)
% hObject    handle to cLogDQ3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of cLogDQ3 as text
%        str2double(get(hObject,'String')) returns contents of cLogDQ3 as a double


% --- Executes during object creation, after setting all properties.
function cLogDQ3_CreateFcn(hObject, eventdata, handles)
% hObject    handle to cLogDQ3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in copyFitting.
function copyFitting_Callback(hObject, eventdata, handles)
% hObject    handle to copyFitting (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
[qCalData, userData, scatt, calibrant]= getData;
dist = 0;

fitRst=getgihandle;
set(handles.edEnergy, 'String', fitRst.xeng);

userData.xeng = str2double(get(handles.edEnergy, 'string'));
userData.wavelength = 12.398418746/userData.xeng;

scatt.wavelength = userData.wavelength;
scatt.xeng = userData.xeng;


qAgbe = str2num(get(handles.qAgbeh, 'String'));
theta = 2.*asin(qAgbe * scatt.wavelength / (4.*pi));

if isempty(fitRst) || ~isfield(fitRst,'center')
    disp('Not Ready! Do the fitting first ....');
    
else    
    scatt.BeamXY = fitRst.center;
    dist = fitRst.px * scatt.pSize;
    scatt.SDD = dist / tan(theta);

    set(handles.edBeamX, 'String', scatt.BeamXY(1));
    set(handles.edBeamY, 'String', scatt.BeamXY(2));
    set(handles.edSDD, 'String', scatt.SDD);

    setData(qCalData, userData, scatt, calibrant);
end



function qAgbeh_Callback(hObject, eventdata, handles)
% hObject    handle to qAgbeh (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of qAgbeh as text
%        str2double(get(hObject,'String')) returns contents of qAgbeh as a double


% --- Executes during object creation, after setting all properties.
function qAgbeh_CreateFcn(hObject, eventdata, handles)
% hObject    handle to qAgbeh (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in combineMask.
function combineMask_Callback(hObject, eventdata, handles)
% hObject    handle to combineMask (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
%
% Note: Mouse: Read mouse clicks: on /off to clear a previous trial
% (1) Mouse: Read Mouse click ON
% (2) press down "Mouse cursor"
% (3) use mouse cursor clicks to make a ROI
% (4) right click mouse then "Toggle off" to finish select
% (5) mouse clicks in "mousepnt.x & .y"

[qCalData, userData, scatt, calibrant]= getData;

% try 
%     mask2 = evalin('base', 'tmpmask');
%     figure; imagesc(mask2); set(gca, 'YData', 'normal');
%     if ~isfield(scatt, 'mask') || isempty(scatt.mask)
%         disp('Initial Mask is not ready! ...')
%         disp('Selected area will be calculated w/o initial mask! ...')
%         scatt.mask = mask2;
%     else
%         scatt.mask = scatt.mask & mask2;
%     end
%     
% catch
%     disp('Select interested area in SAXSimageviewer window ....');
% end
try
    figH = evalin('base', 'SAXSimageviewerhandle');
catch
    fprintf('Start SAXSimageView first...\n');
    return;
end
h2 = guidata(figH);

try
    pnt = evalin('base', 'mousepnt');
    x=[pnt.x pnt.x(1)]; y=[pnt.y pnt.y(1)]; 
    cmd1=sprintf('line(x, y,  ''parent'', h2.ImageAxes, ''Color'', ''m'', ''LineStyle'', ''-'', ''linewidth'', 1);');    
    eval(cmd1); 
catch
    disp('pick the mask / region of interest first');
    return;
end    

imgx = [];
imgy = [];
col = scatt.col; row = scatt.row;
tic
for kk=1:col
    imgx = [imgx ones(1,row)*kk];
    imgy = [imgy 1:row];
end

toc

in = inpolygon(imgx,imgy,x,y);
scatt.arbitMask = reshape(in,[row, col]);
if scatt.arbitMaskExclude
    scatt.arbitMask = 1 - scatt.arbitMask;
end
figure(3004); imagesc(scatt.arbitMask); set(gca, 'YDir', 'normal');

%     if scatt.arbitMaskOn 
%         %mask2 = reshape(in,[row, col]);
%         scatt.arbitMask = reshape(in,[row, col]);
%     else
%         scatt.arbitMask = [];
%     end


setData(qCalData, userData, scatt, calibrant);


% --- Executes on button press in processCurrentImg.
% need work on hdf5 files 1/23/2022
function processCurrentImg_Callback(hObject, eventdata, handles)
% hObject    handle to processCurrentImg (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
[qCalData, userData, scatt, calibrant]= getData;
phd = 0;
cImg = getgihandle;
datadir = cImg.dir;
file = cImg.imagename;

% if file(end-3)~='.'
%     file = [file '.tif'];
% end


%[file, datadir]=uigetfile({'*.tif';'*.ccd'});
if file==0 return; end

% outputFileName=[file(1:end-3) 'dat'];

cfile = fullfile(datadir, file);
[fpath9, fn9, ext9] = fileparts(file);

ctmp = strsplit(datadir, filesep);
tmpParentDir = fullfile(ctmp{1:end-1});
logPath = fullfile(tmpParentDir,'Log',filesep);    
%logPath = fullfile(datadir(1:end-5),'Log',filesep);
logFileName = ['L' fn9(2:end) '.meta'];
fullLogFile = fullfile(logPath, logFileName);

[phd, Io, eng] = parseMetafile(fullLogFile);


% apply sector mask
if scatt.sectorMaskOn && ~isempty(scatt.sectorMask)
    mask = scatt.mask & scatt.sectorMask;
    if scatt.arbitMaskOn && ~isempty(scatt.arbitMask)
        mask = mask & scatt.arbitMask;
    end
    figure(3003); imagesc(mask); set(gca, 'YDir','normal');
else
    mask = scatt.mask;
    if scatt.arbitMaskOn && ~isempty(scatt.arbitMask)
        mask = mask & scatt.arbitMask;
    end
    figure(3003); imagesc(mask); set(gca, 'YDir','normal');
end

%outputFileName=[file(1:end-3) 'dat'];
outputFileName=[fn9 '.dat'];

cfile = fullfile(datadir, file);
disp(cfile);

%if strcmpi(ext9, '.tif')      % tif file
if strncmpi(ext9, '.tif', 4)      % tif file
    sImg = imgUpsideDn(double(imread(cfile)));
    
elseif strcmpi(ext9, '.h5')   % hdf5 data
    if contains(fn9, 'master')
        fprintf('%s is a master file, will not be proecssed!\n',cfile);       
    else
        if isfield(scatt, 'Beamline')
            beamline = scatt.Beamline;
        else
            beamline = '12-ID-B';
        end
        h5data = readHdf5(cfile, beamline);
        %h5data = readHdf5(cfile);
        sImg = imgUpsideDn(double(h5data.data'));
    end
else
    fprintf('%s is not supported and proecssed!\n',cfile);
end    

figure('NumberTitle','off','Name',cfile);imagesc(sImg, [0 500]);

%sImg = imgUpsideDn(double(imread(cfile)));
%figure('NumberTitle','off','Name',cfile);imagesc(sImg);

isReady = checkReadiness;

if isReady
    data = circavgnew2(sImg, mask, scatt.qCMap, scatt.qRMap, scatt.qArray, scatt.offset, scatt.limits);
    %data = circavgnew2(sImg, scatt.mask, scatt.qCMap, scatt.qRMap, scatt.qArray, scatt.offset, scatt.limits);
    figure('NumberTitle','off','Name',['1d data for: ' cfile]); loglog(data(:,1),data(:,2), '-r');
    
    %[filename, pathname] = uiputfile('test.dat', 'Save file name(*.dat)');
    [filename, pathname] = uiputfile(outputFileName, 'Save file name(*.dat)');
    if filename==0
        return;
    else
        outputFile = [pathname filename];
        % [data(:,1) data(:,2)*absIntCoeff/phd data(:,3)./sqrt(data(:,5))*absIntCoeff/phd
        outputData = [data(:,1) data(:,2)/phd data(:,3)./sqrt(data(:,5))/phd];
%         if phd > 0
%             outputData = [data(:,1) data(:,2)/phd data(:,3)/phd];
%         else
%             outputData = data(:,1:3);
%         end
        dlmwrite(outputFile,outputData,'delimiter','\t','precision','%.6f');
    end
    
end



function edit72_Callback(hObject, eventdata, handles)
% hObject    handle to edit72 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit72 as text
%        str2double(get(hObject,'String')) returns contents of edit72 as a double


% --- Executes during object creation, after setting all properties.
function edit72_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit72 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit73_Callback(hObject, eventdata, handles)
% hObject    handle to edit73 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit73 as text
%        str2double(get(hObject,'String')) returns contents of edit73 as a double


% --- Executes during object creation, after setting all properties.
function edit73_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit73 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit74_Callback(hObject, eventdata, handles)
% hObject    handle to edit74 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit74 as text
%        str2double(get(hObject,'String')) returns contents of edit74 as a double


% --- Executes during object creation, after setting all properties.
function edit74_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit74 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit75_Callback(hObject, eventdata, handles)
% hObject    handle to edit75 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit75 as text
%        str2double(get(hObject,'String')) returns contents of edit75 as a double


% --- Executes during object creation, after setting all properties.
function edit75_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit75 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit76_Callback(hObject, eventdata, handles)
% hObject    handle to edit76 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit76 as text
%        str2double(get(hObject,'String')) returns contents of edit76 as a double


% --- Executes during object creation, after setting all properties.
function edit76_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit76 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit77_Callback(hObject, eventdata, handles)
% hObject    handle to edit77 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit77 as text
%        str2double(get(hObject,'String')) returns contents of edit77 as a double


% --- Executes during object creation, after setting all properties.
function edit77_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit77 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit78_Callback(hObject, eventdata, handles)
% hObject    handle to edit78 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit78 as text
%        str2double(get(hObject,'String')) returns contents of edit78 as a double


% --- Executes during object creation, after setting all properties.
function edit78_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit78 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit79_Callback(hObject, eventdata, handles)
% hObject    handle to edit79 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit79 as text
%        str2double(get(hObject,'String')) returns contents of edit79 as a double


% --- Executes during object creation, after setting all properties.
function edit79_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit79 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edCX_Callback(hObject, eventdata, handles)
% hObject    handle to edCX (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edCX as text
%        str2double(get(hObject,'String')) returns contents of edCX as a double


% --- Executes during object creation, after setting all properties.
function edCX_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edCX (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edCY_Callback(hObject, eventdata, handles)
% hObject    handle to edCY (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edCY as text
%        str2double(get(hObject,'String')) returns contents of edCY as a double


% --- Executes during object creation, after setting all properties.
function edCY_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edCY (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edDist_Callback(hObject, eventdata, handles)
% hObject    handle to edDist (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edDist as text
%        str2double(get(hObject,'String')) returns contents of edDist as a double


% --- Executes during object creation, after setting all properties.
function edDist_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edDist (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edQValue_Callback(hObject, eventdata, handles)
% hObject    handle to edQValue (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edQValue as text
%        str2double(get(hObject,'String')) returns contents of edQValue as a double


% --- Executes during object creation, after setting all properties.
function edQValue_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edQValue (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edAstart_Callback(hObject, eventdata, handles)
% hObject    handle to edAstart (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edAstart as text
%        str2double(get(hObject,'String')) returns contents of edAstart as a double


% --- Executes during object creation, after setting all properties.
function edAstart_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edAstart (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edAend_Callback(hObject, eventdata, handles)
% hObject    handle to edAend (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edAend as text
%        str2double(get(hObject,'String')) returns contents of edAend as a double


% --- Executes during object creation, after setting all properties.
function edAend_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edAend (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edDA_Callback(hObject, eventdata, handles)
% hObject    handle to edDA (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edDA as text
%        str2double(get(hObject,'String')) returns contents of edDA as a double


% --- Executes during object creation, after setting all properties.
function edDA_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edDA (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edDPixel_Callback(hObject, eventdata, handles)
% hObject    handle to edDPixel (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edDPixel as text
%        str2double(get(hObject,'String')) returns contents of edDPixel as a double


% --- Executes during object creation, after setting all properties.
function edDPixel_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edDPixel (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in fbPilatus2MCalib.
function fbPilatus2MCalib_Callback(hObject, eventdata, handles)
% hObject    handle to fbPilatus2MCalib (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% This function find multiple pnts on a diffraction ring and stack them in
% the
% updateQarray(handles)
 [qCalData, userData, scatt, calibrant]= getData;
 
cX = str2double(get(handles.edCX, 'string'));
cY = str2double(get(handles.edCY, 'string'));
dist = str2double(get(handles.edDist, 'string'));
qValue = str2double(get(handles.edQValue, 'string'));
Astart = str2double(get(handles.edAsart, 'string'));
Aend = str2double(get(handles.edAend, 'string'));
dA = str2double(get(handles.edDA, 'string'));
dPixel = str2double(get(handles.edDPixel, 'string'));
angles = Astart:dA:Aend;
aLen = length(angles);

try 
    mousePts = evalin('base', 'mousepnt');
    %cImgName = evalin('base', 'cImgName');
    %cImg = imgUpsideDn(double(imread(cImgName)));
    cImgName2 = getgihandle;
    cImg = cImgName2.image;
    [row, col] = size(cImg);
    bnds =[col row];
    halfW = 30;
    pks=[];
    rsqVal = 0.95;
    for kk = 1:length(mousePts.x)

        x=round(mousePts.x(kk));
        y=round(mousePts.y(kk));
        cpt = [x y];
        pkV = findVertPeak(cpt, cImg, halfW, bnds, rsqVal);
        pkH = findHorzPeak(cpt, cImg, halfW, bnds, rsqVal);
        newPKs = [pkV pkH];
        if isempty(newPKs) 
            newPKs = [x;y];
        end
        pks = [pks newPKs];
    end
    
    cdat.xy = sortrows(pks')';
   
    %cdat.xy = [mousePts.x; mousePts.y];
    
    cq = input('Please enter the q-value for the selected ring:');
    cdat.q=cq;

    calibrant.data = [calibrant.data cdat]
    qCalData.calibrant = calibrant;
    assignin('base','qCalData', qCalData);
catch
    disp('Pick calibrant pixels in SAXSimageviewer window ....');
    calibrant.data = [];
    calibrant.wavelength = 0.;
    qCalData.calibrant = qCalData.userData.wavelength;
    assignin('base','qCalData', qCalData);
end



% if isfield(scatt,'qMap') & isfield(scatt, 'qArray') & ~isempty(scatt.qMap) & ~isempty(scatt.qArray)
%     scatt.qRMap = calQRMap(scatt.qMap, scatt.qArray);
%     disp('Q index map done!')
% else
%     disp('Not Ready! Calculate Q and Corr Maps first ....')
% end
setData(qCalData, userData, scatt, calibrant);


% --- Executes on button press in rbPilatus2m.
function rbPilatus2m_Callback(hObject, eventdata, handles)
% hObject    handle to rbPilatus2m (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of rbPilatus2m



function edAbsCoeff_Callback(hObject, eventdata, handles)
% hObject    handle to edAbsCoeff (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edAbsCoeff as text
%        str2double(get(hObject,'String')) returns contents of edAbsCoeff as a double
updateAll(handles);


% --- Executes during object creation, after setting all properties.
function edAbsCoeff_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edAbsCoeff (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in pbPilatus300Calib.
function pbPilatus300Calib_Callback(hObject, eventdata, handles)
% hObject    handle to pbPilatus300Calib (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
agbh_findpeaks_12IDB


% --- Executes on button press in pbDetQCalib.
function pbDetQCalib_Callback(hObject, eventdata, handles)
% hObject    handle to pbDetQCalib (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
detQCalib;
% global detQCalibApp;
% detQCalibApp;
% if isempty(detQCalibApp)
%     detQCalibApp = detQCalib;
% end

% --- Executes on button press in pbMar300Calib.
function pbMar300Calib_Callback(hObject, eventdata, handles)
% hObject    handle to pbMar300Calib (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in pbEiger9mCalib.
function pbEiger9mCalib_Callback(hObject, eventdata, handles)
% hObject    handle to pbEiger9mCalib (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in rbPilatus300.
function rbPilatus300_Callback(hObject, eventdata, handles)
% hObject    handle to rbPilatus300 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of rbPilatus300


% --- Executes during object deletion, before destroying properties.
function rbPilatus300_DeleteFcn(hObject, eventdata, handles)
% hObject    handle to rbPilatus300 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% function writeSetup2Text(fileName, scatt)
% fid = fopen(fileName, 'w');
% fields = fieldnames(scatt);
% for kk =1:length(fields)
%     fieldVal = getfield(scatt, fields{kk});
%     if (strcmp(fields{kk}, 'limits') | strcmp(fields{kk}, 'combQs') |  strcmp(fields{kk}, 'BeamXY'))
%         fprintf(fid, '%+15s: [%s]\n', fields{kk}, num2str(fieldVal));
%     elseif (strcmp(fields{kk}, 'qMap') | strcmp(fields{kk}, 'qCMap') |  strcmp(fields{kk}, 'qArray') |(strcmp(fields{kk}, 'qRMap') | strcmp(fields{kk}, 'mask') )
%         
%     end
%     fprintf(fid, '%+15s: %s\n', fields{kk}, num2str(fieldVal));
% end
% fclose(fid);


% --- Executes on selection change in absIntStandard.
function absIntStandard_Callback(hObject, eventdata, handles)
% hObject    handle to absIntStandard (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns absIntStandard contents as cell array
%        contents{get(hObject,'Value')} returns selected item from absIntStandard
[qCalData, userData, scatt, calibrant]= getData;

index_selected = get(hObject,'Value');
list = get(hObject,'String');
item_selected = list{index_selected};

scatt.absIntStandardIndex =  index_selected;
scatt.absIntStandard = item_selected;

setData(qCalData, userData, scatt, calibrant);

% --- Executes during object creation, after setting all properties.
function absIntStandard_CreateFcn(hObject, eventdata, handles)
% hObject    handle to absIntStandard (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: listbox controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in qStandard.
function qStandard_Callback(hObject, eventdata, handles)
% hObject    handle to qStandard (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns qStandard contents as cell array
%        contents{get(hObject,'Value')} returns selected item from qStandard
[qCalData, userData, scatt, calibrant]= getData;

index_selected = get(hObject,'Value');
list = get(hObject,'String');
item_selected = list{index_selected};

scatt.qStandardIndex =  index_selected;
scatt.qStandard = item_selected;

setData(qCalData, userData, scatt, calibrant);

% --- Executes during object creation, after setting all properties.
function qStandard_CreateFcn(hObject, eventdata, handles)
% hObject    handle to qStandard (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: listbox controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in Beamline.
function Beamline_Callback(hObject, eventdata, handles)
% hObject    handle to Beamline (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns Beamline contents as cell array
%        contents{get(hObject,'Value')} returns selected item from Beamline
[qCalData, userData, scatt, calibrant]= getData;

index_selected = get(hObject,'Value');
list = get(hObject,'String');
item_selected = list{index_selected};

scatt.BeamlineIndex =  index_selected;
scatt.Beamline = item_selected;

setData(qCalData, userData, scatt, calibrant);

% --- Executes during object creation, after setting all properties.
function Beamline_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Beamline (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: listbox controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in selectAvgFolder.
function selectAvgFolder_Callback(hObject, eventdata, handles)
% hObject    handle to selectAvgFolder (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in rdAve.
function rdAve_Callback(hObject, eventdata, handles)
% hObject    handle to rdAve (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of rdAve
[qCalData, userData, scatt, calibrant]= getData;

set(handles.rdAve, 'Value', 1);
userData.avgFolderChoice = 1;
userData.avgFolder='Ave';
set(handles.folderStr, 'String', userData.avgFolder);

setData(qCalData, userData, scatt, calibrant);


% --- Executes on button press in rdAveraged.
function rdAveraged_Callback(hObject, eventdata, handles)
% hObject    handle to rdAveraged (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of rdAveraged
[qCalData, userData, scatt, calibrant]= getData;

set(handles.rdAveraged, 'Value', 1);
userData.avgFolderChoice = 2;
userData.avgFolder='Averaged';
set(handles.folderStr, 'String', userData.avgFolder);

setData(qCalData, userData, scatt, calibrant);


% --- Executes on button press in rdSelectAvgFolder.
function rdSelectAvgFolder_Callback(hObject, eventdata, handles)
% hObject    handle to rdSelectAvgFolder (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of rdSelectAvgFolder
[qCalData, userData, scatt, calibrant]= getData;
avgFolder = uigetdir;
if isnumeric(avgFolder)    
    fprintf('No folder selected!\n');
    if userData.avgFolderChoice==1
        set(handles.rdAve, 'Value', 1);
    else
        set(handles.rdAveraged, 'Value', 1);
    end
    return;
else
    userData.avgFolder= avgFolder;
    set(handles.rdSelectAvgFolder, 'Value', 1);
    userData.avgFolderChoice = 3;
    set(handles.folderStr, 'String', userData.avgFolder);
end

setData(qCalData, userData, scatt, calibrant);


% --- Executes on selection change in AbsIntlistbox.
function AbsIntlistbox_Callback(hObject, eventdata, handles)
% hObject    handle to AbsIntlistbox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns AbsIntlistbox contents as cell array
%        contents{get(hObject,'Value')} returns selected item from AbsIntlistbox
[qCalData, userData, scatt, calibrant]= getData;

index_selected = get(hObject,'Value');
list = get(hObject,'String');
item_selected = list{index_selected};

userData.absIntCalibrantIndex =  index_selected;
userData.absIntCalibrant = item_selected;

switch lower(userData.absIntCalibrant) 
    case 'water'
        set(handles.ed_absIntQmin, 'String', 0.01);
        set(handles.ed_absIntQmax, 'String', 0.8);
    case 'glassy carbon'
        set(handles.ed_absIntQmin, 'String', 0.02);
        set(handles.ed_absIntQmax, 'String', 0.25);        
end

setData(qCalData, userData, scatt, calibrant);

% --- Executes during object creation, after setting all properties.
function AbsIntlistbox_CreateFcn(hObject, eventdata, handles)
% hObject    handle to AbsIntlistbox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: listbox controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function ed_absIntQmax_Callback(hObject, eventdata, handles)
% hObject    handle to ed_absIntQmax (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of ed_absIntQmax as text
%        str2double(get(hObject,'String')) returns contents of ed_absIntQmax as a double


% --- Executes during object creation, after setting all properties.
function ed_absIntQmax_CreateFcn(hObject, eventdata, handles)
% hObject    handle to ed_absIntQmax (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function ed_absIntQmin_Callback(hObject, eventdata, handles)
% hObject    handle to ed_absIntQmin (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of ed_absIntQmin as text
%        str2double(get(hObject,'String')) returns contents of ed_absIntQmin as a double


% --- Executes during object creation, after setting all properties.
function ed_absIntQmin_CreateFcn(hObject, eventdata, handles)
% hObject    handle to ed_absIntQmin (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in absIntSelectData.
function absIntSelectData_Callback(hObject, eventdata, handles)
% hObject    handle to absIntSelectData (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
[qCalData, userData, scatt, calibrant]= getData;
[filename, pathname] = uigetfile( ...
       {'*.dat;*.avg;*.txt', 'All Data Files (*.dat)'}, ...
        'Open as');
calibIntDataFile = [pathname, filename]; 
if calibIntDataFile==0 return; end
Io_water= 0.01632;

qMin = str2num(get(handles.ed_absIntQmin, 'String'));
qMax = str2num(get(handles.ed_absIntQmax, 'String'));

calibData = readData(calibIntDataFile);
if ~isfield(userData,'absIntCalibrantIndex')
    fprintf('Select Calibrant First!\n');
    return;
end

if userData.absIntCalibrantIndex == 1 % None
    userData.absIntCoeff = 1.0;
    set(handles.absIntCoefftext, 'String', sprintf('Coeff.: %.5f', userData.absIntCoeff)); 
elseif userData.absIntCalibrantIndex == 2 % water
    %qMin = get(handles.ed_absIntQmin, 'String')
    [qIndexMin, qIndexMax] = findQrange(calibData(:,1), qMin, qMax);
    x=calibData(qIndexMin:qIndexMax,1);
    y=calibData(qIndexMin:qIndexMax,2);
    p=polyfit(x,y,3);
    disp(sprintf('P3-P0: [%s]; P0=%f.', num2str(p), p(end)));
    yfit=polyval(p,x);
    figure(801); semilogx(calibData(:,1),calibData(:,2),'-r', x,yfit,'-b'); 
    title("Water as Absolute Intensity Calibrant"); xlabel('q, 1/A'); ylabel('scattering intensity');
    userData.absIntCoeff = Io_water/p(end);    
elseif userData.absIntCalibrantIndex == 3 % Glassy Carbon
    gcData= readData(fullfile('GlassyCarbon','GlassyCarbon_Laverage.dat'));
    [qIndexMin, qIndexMax] = findQrange(gcData(:,1), qMin, qMax);
    qgcData=gcData(qIndexMin:qIndexMax,1);
    caliDataInterp = interp1(calibData(:,1), calibData(:,2), qgcData);
    userData.absIntCoeff = mean(gcData(qIndexMin:qIndexMax,2)./caliDataInterp);
    figure(802); loglog(calibData(:,1),calibData(:,2), '-m', calibData(:,1),calibData(:,2)*userData.absIntCoeff,'-r', gcData(:,1),gcData(:,2),'-b'); 
    title("Glassy Carbon as Absolute Intensity Calibrant"); xlabel('q, 1/A'); ylabel('scattering intensity');
end
set(handles.absIntStandard, 'Value', userData.absIntCalibrantIndex);
set(handles.edAbsCoeff, 'String', userData.absIntCoeff);  
scatt.absIntCoeff = userData.absIntCoeff;
list = get(handles.absIntStandard,'String');
%item_selected = list{index_selected};
scatt.absIntStandard = list{userData.absIntCalibrantIndex};
setData(qCalData, userData, scatt, calibrant);


% --- Executes on button press in makeMask.
function makeMask_Callback(hObject, eventdata, handles)
% hObject    handle to makeMask (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
maskmaker


% --- Executes on button press in CopyFromDetQCalib.
function CopyFromDetQCalib_Callback(hObject, eventdata, handles)
% hObject    handle to CopyFromDetQCalib (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
[qCalData, userData, scatt, calibrant]= getData;

set(handles.edEnergy, 'String', scatt.eng);
setDetectorSelection(handles, scatt.isPilatus2m);
setDetectorInitials(handles, scatt.isPilatus2m);

set(handles.edBeamX, 'String', scatt.BeamXY(1));
set(handles.edBeamY, 'String', scatt.BeamXY(2));
set(handles.edSDD, 'String', scatt.SDD);
set(handles.edYaw, 'String', scatt.yaw);

updateAll(handles)
%setData(qCalData, userData, scatt, calibrant);

function setDetectorSelection(handles, isPilatus2m)
if isPilatus2m == 1
    set(handles.rbPilatus2m, 'Value', 1);
elseif isPilatus2m == 0
    set(handles.rbPilatus300, 'Value', 1);
elseif isPilatus2m == 2
    set(handles.rbPE, 'Value', 1);
elseif isPilatus2m == 3
    set(handles.rbPilatus900k, 'Value', 1);
elseif isPilatus2m == 4
    set(handles.rbEiger9m, 'Value', 1);
elseif isPilatus2m == -1
    set(handles.rbOther, 'Value', 1);
end


% --- Executes on button press in ConvH5toTiff.
function ConvH5toTiff_Callback(hObject, eventdata, handles)
% hObject    handle to ConvH5toTiff (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

[filenames, pathname, filterindex] = uigetfile( ...
{  '*.h5','SAXS Image File (*.dat)'; ...
   '*.*',  'All Files (*.*)'}, ...
   'Pick a file', ...
   'MultiSelect', 'on');
% filenames=0 if cancelled
if isnumeric(filenames) 
    fprintf('No files selected!\n');
    return;
end
if ischar(filenames)   % only one file was selected.
    filenames={filenames};
end


for nfile = filenames
    tic
    %h5data=[];
    file = nfile{1,1};
%     [~, fn9, ext9] = fileparts(file);
%     %outputFileName=[file(1:end-3) 'dat'];
%     outputFileName=[fn9 '.dat'];

    cfile = fullfile(pathname, file);
    disp(cfile);
    h5conv2tiff2(cfile)
    toc
end



function sectorA1_Callback(hObject, eventdata, handles)
% hObject    handle to sectorA1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of sectorA1 as text
%        str2double(get(hObject,'String')) returns contents of sectorA1 as a double
[qCalData, userData, scatt, calibrant]= getData;
userData.sectorA1 = str2num(get(handles.sectorA1,'String'));
scatt.sectorA1 = userData.sectorA1;
setData(qCalData, userData, scatt, calibrant);

% --- Executes during object creation, after setting all properties.
function sectorA1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to sectorA1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function sectorA2_Callback(hObject, eventdata, handles)
% hObject    handle to sectorA2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of sectorA2 as text
%        str2double(get(hObject,'String')) returns contents of sectorA2 as a double
[qCalData, userData, scatt, calibrant]= getData;
userData.sectorA2 = str2num(get(handles.sectorA2,'String'));
scatt.sectorA2 = userData.sectorA2;
setData(qCalData, userData, scatt, calibrant);

% --- Executes during object creation, after setting all properties.
function sectorA2_CreateFcn(hObject, eventdata, handles)
% hObject    handle to sectorA2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in plotSector.
function plotSector_Callback(hObject, eventdata, handles)
% hObject    handle to plotSector (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
[qCalData, userData, scatt, calibrant]= getData;
try
    figH = evalin('base', 'SAXSimageviewerhandle');
catch
    fprintf('Start SAXSimageView first...\n');
    return;
end
h2 = guidata(figH);
% ImageView handle" h2.ImageAxes
imgH = getgihandle;
imgHrow = imgH.imgsize(1); imgHcol = imgH.imgsize(2);
row = scatt.row; col = scatt.col;
if (row ~= imgHrow || col ~= imgHcol)
    disp('image size not match!');
end

beamXY = scatt.BeamXY;
nPts = 100;

[xy1, vx1]=genSectorLine(beamXY, imgH.imgsize, scatt.sectorA1, nPts);
[xy2, vx2]=genSectorLine(beamXY, imgH.imgsize, scatt.sectorA2, nPts);

cmd1=sprintf('line(xy1(:,1), xy1(:,2),  ''parent'', h2.ImageAxes, ''Color'', ''m'', ''LineStyle'', ''-'', ''linewidth'', 1);');
cmd2=sprintf('line(xy2(:,1), xy2(:,2),  ''parent'', h2.ImageAxes, ''Color'', ''m'', ''LineStyle'', ''-'', ''linewidth'', 1);');
eval(cmd1); eval(cmd2);
setData(qCalData, userData, scatt, calibrant);

function [xy, vx] = genSectorLine(beamxy, imgsize, angle, nPts)
    %vx = []; vx polygon voxels in anti-clock directions, 
    col=imgsize(2); row = imgsize(1);
    theta1=getSectorAngle(beamxy,[col row]);
    theta2=getSectorAngle(beamxy,[1 row]);
    theta3=getSectorAngle(beamxy,[1 1]);
    theta4=getSectorAngle(beamxy,[col 1]);
    %vx = [theta1 col row; theta2 1 row; theta3 1 1; theta4 col 1 ];
    
    if angle <= theta1
        x=col; y = ceil(beamxy(2) + (x-beamxy(1))*tand(angle));
        vx = [-1 beamxy(1) beamxy(2); angle x y;  theta1 col row; theta2 1 row; theta3 1 1; theta4 col 1 ];
    elseif angle <= theta2
        if angle == 90
            y = row; x=beamxy(1);
        else
            y = row; x = round(beamxy(1) + (y-beamxy(2))/tand(angle));
        end
        vx = [-1 beamxy(1) beamxy(2); angle x y;   theta2 1 row; theta3 1 1; theta4 col 1; theta1 col row];
    elseif angle <= theta3
        x = 1; y = ceil(beamxy(2) + (x-beamxy(1))*tand(angle)); 
        vx = [-1 beamxy(1) beamxy(2); angle x y;   theta3 1 1; theta4 col 1; theta1 col row; theta2 1 row];
    elseif angle <= theta4
        if angle == 270
            y = 1; x=beamxy(1);
        else
            y = 1; x = round(beamxy(1) + (y-beamxy(2))/tand(angle));
        end      
        vx = [-1 beamxy(1) beamxy(2); angle x y;  theta4 col 1; theta1 col row; theta2 1 row;  theta3 1 1];
    else
        x=col; y = ceil(beamxy(2) + (x-beamxy(1))*tand(angle));
        vx = [-1 beamxy(1) beamxy(2); angle x y;  theta1 col row; theta2 1 row;  theta3 1 1; theta4 col 1];
    end
    
    xx=linspace(beamxy(1), x, nPts);
    yy=linspace(beamxy(2), y, nPts);
    xy = [xx' yy'];
    
% calculate sector angle for a certain point(x1,y1)
function angdeg = getSectorAngle(beamxy, pnt)
    angdeg = 0;
    x0 = beamxy(1); y0 = beamxy(2);
    x1 = pnt(1);    y1 = pnt(2);
    if x1 == x0 
        if y1 > y0 
            angdeg = 90;
        else
            angdeg = 270;
        end
    else
        if (x1-x0)>0 && (y1-y0)>=0
            angdeg = atand((y1-y0)/(x1-x0));
        elseif (x1-x0)<0 && (y1-y0)>=0
            angdeg = atand((y1-y0)/(x1-x0)) + 180;
        elseif (x1-x0)<0 && (y1-y0)<=0
            angdeg = atand((y1-y0)/(x1-x0)) + 180;
        elseif (x1-x0)>0 && (y1-y0)<=0
            angdeg = atand((y1-y0)/(x1-x0)) + 360;
        end
    end



% --- Executes on button press in useSectorMask.
function useSectorMask_Callback(hObject, eventdata, handles)
% hObject    handle to useSectorMask (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of useSectorMask
[qCalData, userData, scatt, calibrant]= getData;
userData.sectorMaskOn = get(handles.useSectorMask,'Value');
scatt.sectorMaskOn = userData.sectorMaskOn;
%imgH = getgihandle;
setData(qCalData, userData, scatt, calibrant);

genSectorMask(handles);

% --- Executes on button press in excludeSectorMask.
function excludeSectorMask_Callback(hObject, eventdata, handles)
% hObject    handle to excludeSectorMask (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of excludeSectorMask
[qCalData, userData, scatt, calibrant]= getData;
userData.sectorMaskExclude = get(handles.excludeSectorMask,'Value');
scatt.sectorMaskExclude = userData.sectorMaskExclude;
setData(qCalData, userData, scatt, calibrant);

genSectorMask(handles);


% function genSectorMask3(voxels)
% [qCalData, userData, scatt, calibrant]= getData;
% imgx = [];
% imgy = [];
% xv = voxels(:,1); yv =voxels(:,2);
% col = scatt.col; row = scatt.row;
% for kk=1:col
%     imgx = [imgx ones(1,row)*kk];
%     imgy = [imgy 1:row];
% end
% 
% in = inpolygon(imgx,imgy,xv,yv);
% mask2 = reshape(in,[row, col]);
% 
% if scatt.sectorMaskOn
%     if scatt.sectorMaskExclude
%         scatt.sectorMask = 1-mask2;
%     else
%         scatt.sectorMask = mask2;
%     end
% else
%     scatt.sectorMask = [];
% end
% setData(qCalData, userData, scatt, calibrant);

% generate Sector mask
function genSectorMask(handles)
[qCalData, userData, scatt, calibrant]= getData;
tic
sectorMaskOn = scatt.sectorMaskOn;
if ~sectorMaskOn
    scatt.sectorMask = [];
else
    try
        figH = evalin('base', 'SAXSimageviewerhandle');
    catch
        fprintf('Start SAXSimageView first...\n');
        return;
    end
    %h2 = guidata(figH);
    % ImageView handle" h2.ImageAxes
    imgH = getgihandle;
    imgHrow = imgH.imgsize(1); imgHcol = imgH.imgsize(2);
    row = scatt.row; col = scatt.col;
    if (row ~= imgHrow || col ~= imgHcol)
        disp('image size not match! use qCalibration2 parameters!');
    end

    beamXY = scatt.BeamXY;
    nPts = 100;

    [~, vx1]=genSectorLine(beamXY, imgH.imgsize, scatt.sectorA1, nPts);
    [~, vx2]=genSectorLine(beamXY, imgH.imgsize, scatt.sectorA2, nPts);

    voxels = genPolygonVoxels(scatt.sectorA1, scatt.sectorA2, vx1, vx2 );
    tic
    imgx = [];
    imgy = [];
    xv = voxels(:,1); yv =voxels(:,2);
    col = scatt.col; row = scatt.row;
    for kk=1:col
        imgx = [imgx ones(1,row)*kk];
        imgy = [imgy 1:row];
    end
    toc
    in = inpolygon(imgx,imgy,xv,yv);
    mask2 = reshape(in,[row, col]);

% if scatt.sectorMaskOn
    if scatt.sectorMaskExclude
        scatt.sectorMask = 1-mask2;
        figure(3001); imagesc(scatt.sectorMask); set(gca, 'YDir','normal');
    else
        scatt.sectorMask = mask2;
        figure(3002); imagesc(scatt.sectorMask); set(gca, 'YDir','normal');
    end
% else
%     scatt.sectorMask = [];
end
toc
setData(qCalData, userData, scatt, calibrant);


function voxels=genPolygonVoxels(A1, A2, vx1, vx2 )
ang1 = vx1(:,1);
ang2 = vx2(:,1);
AA1 = vx1(3,1);
AA2 = vx2(3,1);
if A1 >= A2
    kk= find(ang1==AA2);
    voxels = [vx1(1:kk-1,2:3); vx2(2,2:3) ];
else
    kk= find(ang2==AA1);
    voxels = [vx2(1:kk-1,2:3); vx1(2,2:3) ];
end


% --- Executes on button press in arbitMaskOn.
function arbitMaskOn_Callback(hObject, eventdata, handles)
% hObject    handle to arbitMaskOn (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of arbitMaskOn

[qCalData, userData, scatt, calibrant]= getData;
userData.arbitMaskOn = get(handles.arbitMaskOn,'Value');
scatt.arbitMaskOn = userData.arbitMaskOn;
setData(qCalData, userData, scatt, calibrant);


% --- Executes on button press in arbitMaskExclude.
function arbitMaskExclude_Callback(hObject, eventdata, handles)
% hObject    handle to arbitMaskExclude (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of arbitMaskExclude

[qCalData, userData, scatt, calibrant]= getData;
userData.arbitMaskExclude = get(handles.arbitMaskExclude,'Value');
scatt.arbitMaskExclude = userData.arbitMaskExclude;
setData(qCalData, userData, scatt, calibrant);


% --- Executes on button press in rbPilatus900k.
function rbPilatus900k_Callback(hObject, eventdata, handles)
% hObject    handle to rbPilatus900k (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of rbPilatus900k
