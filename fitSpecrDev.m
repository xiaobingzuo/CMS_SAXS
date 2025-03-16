function varargout = fitSpecrDev(varargin)
% FITSPECRDEV MATLAB code for fitSpecrDev.fig
%      FITSPECRDEV, by itself, creates a new FITSPECRDEV or raises the existing
%      singleton*.
%
%      H = FITSPECRDEV returns the handle to a new FITSPECRDEV or the handle to
%      the existing singleton*.
%
%      FITSPECRDEV('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in FITSPECRDEV.M with the given input arguments.
%
%      FITSPECRDEV('Property','Value',...) creates a new FITSPECRDEV or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before fitSpecrDev_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to fitSpecrDev_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help fitSpecrDev

% Last Modified by GUIDE v2.5 30-Nov-2024 16:05:09

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @fitSpecrDev_OpeningFcn, ...
                   'gui_OutputFcn',  @fitSpecrDev_OutputFcn, ...
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


% --- Executes just before fitSpecrDev is made visible.
function fitSpecrDev_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to fitSpecrDev (see VARARGIN)

% Choose default command line output for fitSpecrDev
handles.output = hObject;

handles.dat={};
handles.dat.XData = [];
handles.dat.YData = [];
handles.dat.Xrange=[];
handles.dat.XrangeIndex=[];
handles.dat.useDataRange = 0;

handles.dat.xp = [];
handles.dat.xp_lb = [];
handles.dat.xp_ub = [];
handles.dat.yfit=[];
handles.dat.fit_r= 0;
handles.dat.figTag ='specr_Fig'; % default for specr fig
handles.dat.figType = 1; % 1, specr; 2, SAXS_Lee; 3, generic.
handles.dat.curveTag='';
handles.dat.hFig = [];
handles.dat.hFig = findall(0, 'Tag', handles.dat.figTag);

handles.dat.fitFig_tag = 'fit_curve_fig';


handles=update_handles(hObject, eventdata, handles);

% Update handles structure
guidata(hObject, handles);
% UIWAIT makes fitSpecrDev wait for user response (see UIRESUME)
% uiwait(handles.figure1);


function handles=update_handles(hObject, eventdata, handles, curveTag)
% hObject    handle to bkg (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
if ~isempty(handles.dat.hFig)
    if nargin<4
        curveTag='';
    end
    childIndex=0;
    %hFigSpecr = findall(0,'Tag','specr_Fig');

    %hFigSpecr = findall(0,'Tag', handles.dat.figTag);
    %if ~isempty(hFigSpecr)
    %hFigSpecr =
    h = handles.dat.hFig.CurrentAxes.Children;    
    childIndex=0;
    if isempty(curveTag)
        childIndex=1;
        curveTag = h.Tag;
    else
        for j=1:length(h)
            if strcmp(h(j).Tag, curveTag)
                childIndex = j;
            end
        end
    end

    if childIndex>0 
        x=h(childIndex).XData;
        y=h(childIndex).YData;
        curveTag = h(childIndex).Tag;
        handles.dat.XData=x;
        handles.dat.YData=y;
        handles.dat.curveTag = curveTag;

        %C = str2double(regexprep(C,'[a-zA-Z^()]','')).' 
        %D = str2double( C(ismember(C, '0':'9')) );
        if (handles.dat.figType==1)
            curveType = regexprep(curveTag, '\d*', '');
            serN = round(str2double(regexprep(curveTag,'[a-zA-Z^()]','')));
            handles.dat.curveType = curveType;
            handles.dat.serN = serN;
            if strcmp(curveType(6:end), 'Dev')
                set(handles.rdDev, 'Value', 1);
                set(handles.rdScan, 'Value', 0);
            else
                set(handles.rdScan, 'Value', 1);
                set(handles.rdDev, 'Value', 0);
            end
            
            set(handles.scanNum, 'String', serN);
        end

        f=figure(666); 
        plot(x,y,'-ob', 'LineWidth',1.5);
        title(sprintf('Gaussian Fit for %s', curveTag(1:(min(15,length(curveTag))))));
        f.Tag=handles.dat.fitFig_tag;

        
        % initial parameters
        handles.dat.Xrange=[min(x) max(x)];
        handles.dat.XrangeIndex=[1 length(x)];

        xspan = abs(max(x) - min(x));
        yspan = abs(max(y) -min(y));

        x0 = mean(x);
        bkg = (y(1) + y(end))/2.;
        if min(y)+max(y) < 2*mean(y) % mean(y) - min(y) > max(y)- mean(y)
            a = -1.*yspan;
        else
            a = yspan;
        end
        sigma = xspan / 4.;

        xp = [a x0 sigma bkg];

        x0_lb = x0 - xspan /3.;
        x0_ub = x0 + xspan /3.;
        
        a_lb = a - yspan /2.;
        a_ub = a + yspan /2.;

        sigma_lb = 0.;
        sigma_ub = xspan / 2.;
        
        bkg_lb = min(y);
        bkg_ub = max(y);

        set(handles.a, 'String', sprintf('%.2f',a));
        set(handles.a_lb, 'String', a_lb);
        set(handles.a_ub, 'String', a_ub);

        set(handles.x0, 'String', x0);
        set(handles.x0_lb, 'String', x0_lb);
        set(handles.x0_ub, 'String', x0_ub);

        set(handles.sigma, 'String', sigma);
        set(handles.sigma_lb, 'String', sigma_lb);
        set(handles.sigma_ub, 'String', sigma_ub);

        set(handles.bkg, 'String', bkg);
        set(handles.bkg_lb, 'String', bkg_lb);
        set(handles.bkg_ub, 'String', bkg_ub);
        
        xp_lb = [a_lb x0_lb sigma_lb bkg_lb];
        xp_ub = [a_ub x0_ub sigma_ub bkg_ub];

        handles.dat.xp = xp;
        handles.dat.xp_lb = xp_lb;
        handles.dat.xp_ub = xp_ub;

        figure(f);
%         else
%     
%         end
    end
end
guidata(hObject, handles);


% --- Outputs from this function are returned to the command line.
function varargout = fitSpecrDev_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on button press in PlotSelectedCurve.
function PlotSelectedCurve_Callback(hObject, eventdata, handles)
% hObject    handle to PlotSelectedCurve (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
if (handles.dat.figType==1)
    serN = round(str2num(get(handles.scanNum, 'String')));
    if get(handles.rdDev, 'Value')
        prefix ='specrDev';
    else
        prefix ='specrScan';
    end
    curveTag = sprintf('%s%d', prefix, serN);
    handles = update_handles(hObject, eventdata, handles, curveTag);
else
    serN = round(str2double(get(handles.curveN, 'String')));
    h = handles.dat.hFig.CurrentAxes.Children; 
    try
        curveTag = h(serN).Tag;
        handles = update_handles(hObject, eventdata, handles, curveTag);
    catch
        fprintf('something wrong! 2222');
    end
end
guidata(hObject, handles);


% hFigSpecr = findall(0,'Tag','specr_Fig');
% if ~isempty(hFigSpecr)
%     h = hFigSpecr.CurrentAxes.Children;    
%     index=0;
%     for j=1:length(h)
% 	    if strcmp(h(j).Tag, curveStr)
%             index = j;
%         end
%     end
%     if index>0 
%         x=h(index).XData;
%         y=h(index).YData;
%         f=figure(666); 
%         plot(x,y,'-ob', 'LineWidth',1.5);
%         title(sprintf('Gaussian Fit for %s', curveStr));
%         f.Tag='specr_fig_fit';
%     end
% 
% end


function scanNum_Callback(hObject, eventdata, handles)
% hObject    handle to scanNum (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of scanNum as text
%        str2double(get(hObject,'String')) returns contents of scanNum as a double


% --- Executes during object creation, after setting all properties.
function scanNum_CreateFcn(hObject, eventdata, handles)
% hObject    handle to scanNum (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function bkg_Callback(hObject, eventdata, handles)
% hObject    handle to bkg (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of bkg as text
%        str2double(get(hObject,'String')) returns contents of bkg as a double


% --- Executes during object creation, after setting all properties.
function bkg_CreateFcn(hObject, eventdata, handles)
% hObject    handle to bkg (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function x0_Callback(hObject, eventdata, handles)
% hObject    handle to x0 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of x0 as text
%        str2double(get(hObject,'String')) returns contents of x0 as a double


% --- Executes during object creation, after setting all properties.
function x0_CreateFcn(hObject, eventdata, handles)
% hObject    handle to x0 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function sigma_Callback(hObject, eventdata, handles)
% hObject    handle to sigma (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of sigma as text
%        str2double(get(hObject,'String')) returns contents of sigma as a double


% --- Executes during object creation, after setting all properties.
function sigma_CreateFcn(hObject, eventdata, handles)
% hObject    handle to sigma (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function a_Callback(hObject, eventdata, handles)
% hObject    handle to a (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of a as text
%        str2double(get(hObject,'String')) returns contents of a as a double


% --- Executes during object creation, after setting all properties.
function a_CreateFcn(hObject, eventdata, handles)
% hObject    handle to a (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in FitWithGaussian.
function FitWithGaussian_Callback(hObject, eventdata, handles)
% hObject    handle to FitWithGaussian (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

bkg = str2double(get(handles.bkg, 'String'));
x0 = str2double(get(handles.x0, 'String'));
sigma = str2double(get(handles.sigma, 'String'));
a = str2double(get(handles.a, 'String'));
c1 = str2double(get(handles.c1, 'String'));

bkg_lb = str2double(get(handles.bkg_lb, 'String'));
x0_lb = str2double(get(handles.x0_lb, 'String'));
sigma_lb = str2double(get(handles.sigma_lb, 'String'));
a_lb = str2double(get(handles.a_lb, 'String'));
c1_lb = str2double(get(handles.c1_lb, 'String'));

bkg_ub = str2double(get(handles.bkg_ub, 'String'));
x0_ub = str2double(get(handles.x0_ub, 'String'));
sigma_ub = str2double(get(handles.sigma_ub, 'String'));
a_ub = str2double(get(handles.a_ub, 'String'));
c1_ub = str2double(get(handles.c1_ub, 'String'));

% h = findall(0,'Tag', handles.dat.fitFig_tag); %'specr_fig_fit');
% figure(h);
xp0 =[a x0 sigma bkg c1];
xp_lb =[a_lb x0_lb sigma_lb bkg_lb c1_lb];
xp_ub =[a_ub x0_ub sigma_ub bkg_ub c1_ub];

if ~isempty(handles.dat.XData)
    %dat=h.CurrentAxes.Children;
    %x= dat(length(dat)).XData;
    %y= dat(length(dat)).YData;
    
    if handles.dat.useDataRange
        d1=handles.dat.XrangeIndex(1);
        d2=handles.dat.XrangeIndex(2);
        x=handles.dat.XData(d1:d2);
        y=handles.dat.YData(d1:d2);
    else
        x=handles.dat.XData; 
        y=handles.dat.YData;         
    end

    % y(x) = a*exp(-(x-x0)^2/(2*s^2)) + bkg + c1*x
%     rst= sum((y-(xp0(1)*exp(-(x-xp0(2)).^2/(2*xp0(3)^2)) + xp0(4))).^2)
%     % 
%     fun = @(xp)sum((y-(xp(1)*exp(-(x-xp(2)).^2/(2*xp(3)^2)) + xp(4))).^2);
%     [xp, fval] = fminsearchcon(fun, xp0, xp_lb, xp_ub);
%     
%     
%     set(handles.bkg, 'String', sprintf('%.2f', xp(4)));
%     set(handles.x0, 'String', sprintf('%.3f', xp(2)));
%     set(handles.a, 'String', sprintf('%.2f', xp(1)));
%     set(handles.sigma, 'String', sprintf('%.4f', xp(3)));

    rst= sum((y-(xp0(1)*exp(-(x-xp0(2)).^2/(2*xp0(3)^2)) + xp0(4)+ xp0(5)*x)).^2)
    % 
    fun = @(xp)sum((y-(xp(1)*exp(-(x-xp(2)).^2/(2*xp(3)^2)) + xp(4) + xp(5)*x)).^2);
    [xp, fval] = fminsearchcon(fun, xp0, xp_lb, xp_ub);
    
    
    set(handles.bkg, 'String', sprintf('%.4f', xp(4)));
    set(handles.x0, 'String', sprintf('%.4f', xp(2)));
    set(handles.a, 'String', sprintf('%.4f', xp(1)));
    set(handles.sigma, 'String', sprintf('%.4f', xp(3)));
    set(handles.c1, 'String', sprintf('%.4f', xp(5)));


    h = findall(0,'Tag', handles.dat.fitFig_tag); %'specr_fig_fit');
    figure(h);
    
    hold on;
    x_sim=linspace(min(x), max(x), 200);
    y_sim=xp(1)*exp(-(x_sim-xp(2)).^2/(2*xp(3)^2)) + xp(4) + xp(5)*x_sim;
    y_fit = xp(1)*exp(-(x-xp(2)).^2/(2*xp(3)^2)) + xp(4) + xp(5)*x;
    ssr = sum(abs(y_fit-mean(y)).^2);
    sse = sum(abs(y_fit-y).^2);
    r=sqrt(ssr/(sse+ssr));

    plot(x_sim,y_sim, '-m', 'LineWidth',1);
    hold on;
    plot(x,y_fit, '-r', 'LineWidth',2);
    hold off;




    % legend({});
    f.name = 'a*exp(-(x-x0)^2/(2*sigma^2)) + bkg + c1*x'
    % f.yvar='y';
    % f.xvar='x';
    % f.fitmode='fminsearch';
    f.eq = 'y=a*exp(-(x-x0)^2/(2*sigma^2)) + bkg + c1*x';
    f.param = {'a';  'x0';  'sigma';   'bkg'; 'c1'};
    f.m = xp;
    % f.m0=xp0;
    % f.x =x;
    % f.y = y1;
    f.r = r;
    rst_a_str= sprintf('  %s=%.4f',f.param{1}, f.m(1))
    rst_x0_str= sprintf('  %s=%.4f',f.param{2}, f.m(2));
    rst_sigma_str= sprintf('  %s=%.5f',f.param{3}, f.m(3));
    rst_bkg_str= sprintf('  %s=%.4f',f.param{4}, f.m(4));
    rst_c1_str= sprintf('  %s=%.4f',f.param{5}, f.m(5));
    rst_r_str = sprintf('  R=%.4f', f.r);
    rst_fwhm_str= sprintf('  FWHM=%.5f', f.m(3)*2.355);

    fitRst={'Result of Gaussian Fit:'; f.eq; rst_a_str; rst_x0_str; rst_sigma_str; rst_bkg_str; rst_c1_str; rst_r_str; rst_fwhm_str};

    annotation('textbox',...
    [0.55 0.45 0.1 0.1],...
    'String',fitRst,...
    'FontSize',10,...
    'FontName','Arial',...
    'LineStyle','-',...
    'LineWidth',0.5,...
    'Color','Red');
    
    %showfit(f,'fitcolor', 'red', 'fitlinewidth', 1);

    % 
    % fitStr = sprintf('y(x) = a*exp(-(x-x0)^2/(2*s^2)) + bkg; x0=%f; s=%f; bkg=%f; a=%f', x0, sigma, bkg, a);
    % sprintf(fitStr);
    % f = ezfit(x,y, fitStr);
    % 
    % showfit(f, 'fitcolor', 'red', 'fitlinewidth', 1);
    % 
    % htb=findall(gcf,'Type','TextBox')
    % %htb(length(htb)).Position(1)=0.55
    % %htb(length(htb)).Position(2)=0.2
    % 
    % htb(1).Position(1)=0.55
    % htb(1).Position(2)=0.2
    
    handles.dat.xp=xp;
    handles.dat.xp_lb=xp_lb;
    handles.dat.xp_ub=xp_ub;
    handles.dat.yfit=y_fit;
    handles.dat.fit_r=r;
end

guidata(hObject, handles);



% --- Executes on button press in RemoveTextBox.
function RemoveTextBox_Callback(hObject, eventdata, handles)
% hObject    handle to RemoveTextBox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
h = findall(0,'Tag', handles.dat.fitFig_tag); %'specr_fig_fit');

figure(h);
htb=findall(gcf,'Type','TextBox');
delete(htb);


function bkg_lb_Callback(hObject, eventdata, handles)
% hObject    handle to bkg_lb (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of bkg_lb as text
%        str2double(get(hObject,'String')) returns contents of bkg_lb as a double


% --- Executes during object creation, after setting all properties.
function bkg_lb_CreateFcn(hObject, eventdata, handles)
% hObject    handle to bkg_lb (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function x0_lb_Callback(hObject, eventdata, handles)
% hObject    handle to x0_lb (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of x0_lb as text
%        str2double(get(hObject,'String')) returns contents of x0_lb as a double


% --- Executes during object creation, after setting all properties.
function x0_lb_CreateFcn(hObject, eventdata, handles)
% hObject    handle to x0_lb (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function sigma_lb_Callback(hObject, eventdata, handles)
% hObject    handle to sigma_lb (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of sigma_lb as text
%        str2double(get(hObject,'String')) returns contents of sigma_lb as a double


% --- Executes during object creation, after setting all properties.
function sigma_lb_CreateFcn(hObject, eventdata, handles)
% hObject    handle to sigma_lb (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function a_lb_Callback(hObject, eventdata, handles)
% hObject    handle to a_lb (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of a_lb as text
%        str2double(get(hObject,'String')) returns contents of a_lb as a double


% --- Executes during object creation, after setting all properties.
function a_lb_CreateFcn(hObject, eventdata, handles)
% hObject    handle to a_lb (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function bkg_ub_Callback(hObject, eventdata, handles)
% hObject    handle to bkg_ub (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of bkg_ub as text
%        str2double(get(hObject,'String')) returns contents of bkg_ub as a double


% --- Executes during object creation, after setting all properties.
function bkg_ub_CreateFcn(hObject, eventdata, handles)
% hObject    handle to bkg_ub (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function x0_ub_Callback(hObject, eventdata, handles)
% hObject    handle to x0_ub (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of x0_ub as text
%        str2double(get(hObject,'String')) returns contents of x0_ub as a double


% --- Executes during object creation, after setting all properties.
function x0_ub_CreateFcn(hObject, eventdata, handles)
% hObject    handle to x0_ub (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function sigma_ub_Callback(hObject, eventdata, handles)
% hObject    handle to sigma_ub (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of sigma_ub as text
%        str2double(get(hObject,'String')) returns contents of sigma_ub as a double


% --- Executes during object creation, after setting all properties.
function sigma_ub_CreateFcn(hObject, eventdata, handles)
% hObject    handle to sigma_ub (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function a_ub_Callback(hObject, eventdata, handles)
% hObject    handle to a_ub (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of a_ub as text
%        str2double(get(hObject,'String')) returns contents of a_ub as a double


% --- Executes during object creation, after setting all properties.
function a_ub_CreateFcn(hObject, eventdata, handles)
% hObject    handle to a_ub (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in xplbub.
function xplbub_Callback(hObject, eventdata, handles)
% hObject    handle to xplbub (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of xplbub
xp_name={'a'; 'x0'; 'sigma'; 'bkg'};
v0=-9999.0;
p_range= 0.3;
datRange = 0;
h = findall(0,'Tag', handles.dat.fitFig_tag); %'specr_fig_fit');
figure(h);
if ~isempty(h)
    dat=h.CurrentAxes.Children;
    x= dat(length(dat)).XData;
    y= dat(length(dat)).YData;
    xrange=max(x)-min(x);
    yrange=max(y)-min(y);

    if get(handles.xplbub, 'Value')
        for j=1:4
            eval(sprintf("v0 = str2double(get(handles.%s, 'String'));", xp_name{j}));
            %v0 = str2double(get(handles.bkg, 'String')); 
            if (strcmp(xp_name{j}, 'a') || strcmp(xp_name{j}, 'bkg'))
                datRange = yrange;
            elseif (strcmp(xp_name{j}, 'sigma') || strcmp(xp_name{j}, 'x0'))
                datRange = xrange;
            end
    
            lb = min(v0 + datRange * p_range * [-1. 1.]);
            ub = max(v0 + datRange * p_range * [-1. 1.]);
            
            %set(handles.bkg, 'String', sprintf('%.2f', xp(4)));
            eval(sprintf("set(handles.%s_lb, 'String', %.4f);", xp_name{j}, lb));
            eval(sprintf("set(handles.%s_ub, 'String', %.4f);", xp_name{j}, ub));
        end
    end
end


% --- If Enable == 'on', executes on mouse press in 5 pixel border.
% --- Otherwise, executes on mouse press in 5 pixel border or over xplbub.
function xplbub_ButtonDownFcn(hObject, eventdata, handles)
% hObject    handle to xplbub (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)



% --- Executes on button press in pushbutton4.
function pushbutton4_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton4 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
h = findall(0,'Tag', handles.dat.fitFig_tag); %'specr_fig_fit');
figure(h);
if ~isempty(h)
    if length(h.CurrentAxes.Children) > 1
        delete(h.CurrentAxes.Children(1))
    end
end


% --- Executes on button press in useDataRange.
function useDataRange_Callback(hObject, eventdata, handles)
% hObject    handle to useDataRange (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of useDataRange

if get(hObject, 'Value')
    dataRange = str2num(get(handles.dataRange, 'String'));
    handles.dat.useDataRange = 1;
    if ~isempty(handles.dat.XData)
        [d1, d2] = findQrange(handles.dat.XData, min(dataRange), max(dataRange));
        handles.dat.XrangeIndex=[d1 d2];
        handles = update_xp_useDataRange(hObject, eventdata, handles);
    end
else
    handles.dat.useDataRange = 0;
end
guidata(hObject, handles);


function handles=update_xp_useDataRange(hObject, eventdata, handles)
% hObject    handle to useDataRange (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
if ~isempty(handles.dat.XData) && ~isempty(handles.dat.YData) && handles.dat.useDataRange
        

    curveTag = handles.dat.curveTag;
    
    % initial parameters
    x_all = handles.dat.XData;
    y_all = handles.dat.YData;
    d1=handles.dat.XrangeIndex(1);
    d2=handles.dat.XrangeIndex(2);

    x = x_all(d1:d2);
    y = y_all(d1:d2);


    xspan = abs(max(x) - min(x));
    yspan = abs(max(y) -min(y));

    x0 = mean(x);
    bkg = (y(1) + y(end))/2.;
    if min(y)+max(y) < 2*mean(y) % mean(y) - min(y) > max(y)- mean(y)
        a = -1.*yspan;
    else
        a = yspan;
    end
    sigma = xspan / 4.;

    xp = [a x0 sigma bkg];

    x0_lb = x0 - xspan /3.;
    x0_ub = x0 + xspan /3.;
    
    a_lb = a - yspan /2.;
    a_ub = a + yspan /2.;

    sigma_lb = 0.;
    sigma_ub = xspan / 2.;
    
    bkg_lb = min(y);
    bkg_ub = max(y);

    set(handles.a, 'String', sprintf('%.2f',a));
    set(handles.a_lb, 'String', a_lb);
    set(handles.a_ub, 'String', a_ub);

    set(handles.x0, 'String', x0);
    set(handles.x0_lb, 'String', x0_lb);
    set(handles.x0_ub, 'String', x0_ub);

    set(handles.sigma, 'String', sigma);
    set(handles.sigma_lb, 'String', sigma_lb);
    set(handles.sigma_ub, 'String', sigma_ub);

    set(handles.bkg, 'String', bkg);
    set(handles.bkg_lb, 'String', bkg_lb);
    set(handles.bkg_ub, 'String', bkg_ub);
    
    xp_lb = [a_lb x0_lb sigma_lb bkg_lb];
    xp_ub = [a_ub x0_ub sigma_ub bkg_ub];

    handles.dat.xp = xp;
    handles.dat.xp_lb = xp_lb;
    handles.dat.xp_ub = xp_ub;

    f=figure(666); 
    plot(x_all,y_all,'-ob', 'LineWidth',1.5); 
    xlim([min(x)-0.05*xspan max(x)+0.05*xspan]);
    ylim([min(y)-0.10*yspan max(y)+0.10*yspan]);
    title(sprintf('Gaussian Fit for %s', curveTag(1:(min(15,length(curveTag))))));
    f.Tag=handles.dat.fitFig_tag;

    figure(f);
end

guidata(hObject, handles);

function dataRange_Callback(hObject, eventdata, handles)
% hObject    handle to dataRange (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of dataRange as text
%        str2double(get(hObject,'String')) returns contents of dataRange as a double


% --- Executes during object creation, after setting all properties.
function dataRange_CreateFcn(hObject, eventdata, handles)
% hObject    handle to dataRange (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in selectFirstOrOnlyCurve.
function selectFirstOrOnlyCurve_Callback(hObject, eventdata, handles)
% hObject    handle to selectFirstOrOnlyCurve (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
%handles.output = hObject;

%handles.dat={};
handles=update_handles(hObject, eventdata, handles);

% Update handles structure
guidata(hObject, handles);


% --- Executes when selected object is changed in uibuttongroup2.
function uibuttongroup2_SelectionChangedFcn(hObject, eventdata, handles)
% hObject    handle to the selected object in uibuttongroup2 
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
if (get(handles.rbSpecr, 'Value'))
    handles.dat.figTag ='specr_Fig'; % default for specr fig
    handles.dat.figType = 1; % 1, specr; 2, SAXS_Lee; 3, generic.
elseif (get(handles.rbSAXSLee, 'Value'))
    handles.dat.figTag ='SAXSLee_Fig'; % default for specr fig
    handles.dat.figType = 2; % 1, specr; 2, SAXS_Lee; 3, generic.    
elseif (get(handles.rbOtherFig, 'Value'))
    %handles.dat.figTag ='other_Fig'; % default for specr fig
    handles.dat.figType = 3; % 1, specr; 2, SAXS_Lee; 3, generic.    
    handles.dat.figTag = get(handles.other_fig_tag, 'String');    
else
    disp('Should not be here! 333');
end
hFig = findall(0,'Tag', handles.dat.figTag);
if isempty(hFig)
    fprintf('Something wrong with the figure %s, please check!', handles.dat.figTag);
end
handles.dat.hFig = hFig;
guidata(hObject, handles);

function c1_Callback(hObject, eventdata, handles)
% hObject    handle to c1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of c1 as text
%        str2double(get(hObject,'String')) returns contents of c1 as a double


% --- Executes during object creation, after setting all properties.
function c1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to c1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function c1_lb_Callback(hObject, eventdata, handles)
% hObject    handle to c1_lb (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of c1_lb as text
%        str2double(get(hObject,'String')) returns contents of c1_lb as a double


% --- Executes during object creation, after setting all properties.
function c1_lb_CreateFcn(hObject, eventdata, handles)
% hObject    handle to c1_lb (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function c1_ub_Callback(hObject, eventdata, handles)
% hObject    handle to c1_ub (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of c1_ub as text
%        str2double(get(hObject,'String')) returns contents of c1_ub as a double


% --- Executes during object creation, after setting all properties.
function c1_ub_CreateFcn(hObject, eventdata, handles)
% hObject    handle to c1_ub (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function curveN_Callback(hObject, eventdata, handles)
% hObject    handle to curveN (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of curveN as text
%        str2double(get(hObject,'String')) returns contents of curveN as a double


% --- Executes during object creation, after setting all properties.
function curveN_CreateFcn(hObject, eventdata, handles)
% hObject    handle to curveN (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function other_fig_tag_Callback(hObject, eventdata, handles)
% hObject    handle to other_fig_tag (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of other_fig_tag as text
%        str2double(get(hObject,'String')) returns contents of other_fig_tag as a double


% --- Executes during object creation, after setting all properties.
function other_fig_tag_CreateFcn(hObject, eventdata, handles)
% hObject    handle to other_fig_tag (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- If Enable == 'on', executes on mouse press in 5 pixel border.
% --- Otherwise, executes on mouse press in 5 pixel border or over useDataRange.
function useDataRange_ButtonDownFcn(hObject, eventdata, handles)
% hObject    handle to useDataRange (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
dataRange = str2num(get(handles.dataRange, 'String'));
if (get(handles.useDataRange, 'Value'))
    handles.dat.useDataRange = 1;
    if ~isempty(handles.dat.XData)
        [d1, d2] = findQrange(handles.dat.XData, min(dataRange), max(dataRange));
        handles.dat.XrangeIndex=[d1 d2];
        update_xp_useDataRange(hObject, eventdata, handles);
    end
else
    handles.dat.useDataRange = 0;
end
guidata(hObject, handles);


% --- Executes on key press with focus on useDataRange and none of its controls.
function useDataRange_KeyPressFcn(hObject, eventdata, handles)
% hObject    handle to useDataRange (see GCBO)
% eventdata  structure with the following fields (see MATLAB.UI.CONTROL.UICONTROL)
%	Key: name of the key that was pressed, in lower case
%	Character: character interpretation of the key(s) that was pressed
%	Modifier: name(s) of the modifier key(s) (i.e., control, shift) pressed
% handles    structure with handles and user data (see GUIDATA)
dataRange = str2num(get(handles.dataRange, 'String'));
if (get(handles.useDataRange, 'Value'))
    handles.dat.useDataRange = 1;
    if ~isempty(handles.dat.XData)
        [d1, d2] = findQrange(handles.dat.XData, min(dataRange), max(dataRange));
        handles.dat.XrangeIndex=[d1 d2];
        update_xp_useDataRange(hObject, eventdata, handles)
    end
else
    handles.dat.useDataRange = 0;
end
guidata(hObject, handles);
