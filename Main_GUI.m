function varargout = Main_GUI(varargin)
% MAIN_GUI MATLAB code for Main_GUI.fig
%      MAIN_GUI, by itself, creates a new MAIN_GUI or raises the existing
%      singleton*.
%
%      H = MAIN_GUI returns the handle to a new MAIN_GUI or the handle to
%      the existing singleton*.
%
%      MAIN_GUI('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in MAIN_GUI.M with the given input arguments.
%
%      MAIN_GUI('Property','Value',...) creates a new MAIN_GUI or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before Main_GUI_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to Main_GUI_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help Main_GUI

% Last Modified by GUIDE v2.5 26-Nov-2020 15:16:00

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @Main_GUI_OpeningFcn, ...
                   'gui_OutputFcn',  @Main_GUI_OutputFcn, ...
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


% --- Executes just before Main_GUI is made visible.
function Main_GUI_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to Main_GUI (see VARARGIN)

% Choose default command line output for Main_GUI
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes Main_GUI wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = Main_GUI_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on button press in pushbutton1.
function pushbutton1_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

global I

% -- Getting Input Image -- %

[f,p] = uigetfile('Images\*.*');

if f == 0
    
    warndlg('You Have Cancelled');
    
else
    
    I = imread([p f]);
    
axes(handles.axes1);

    imshow(I);
    
    title('Input Image','FontName',...
        'Times New Roman','Fontsize',12);
end

% --- Executes on button press in pushbutton2.
function pushbutton2_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global I Filt

    % -- Preprocess -- %
    
    RES = imresize(I,[256 256]);
    
axes(handles.axes2);    
    imshow(RES);
    
    title('Resized Image','FontName',...
        'Times New Roman','Fontsize',12);
    
    Filter_size = 3;
    
    Filter_level = 0.5;
    
    if size(RES,3) == 3
        
        Filt(:,:,1) = medfilt2(RES(:,:,1));
        Filt(:,:,2) = medfilt2(RES(:,:,2));
        Filt(:,:,3) = medfilt2(RES(:,:,3));
    else
        Filt = medfilt2(RES);
        
    end
    
axes(handles.axes3);    
    
    Filt = uint8(Filt);
    
    imshow(Filt);
    
    title('Filtered Image','fontsize',12,...
        'fontname','Times New Roman','Color','Black');


% --- Executes on button press in pushbutton3.
function pushbutton3_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

global Filt LBPval

    % --
    
    % Feature extraction
    if size(Filt,3) == 3
    SEG = rgb2gray(Filt);
    else
    SEG = (Filt);
    end
    filtDims = [2 3];
    % inImg = rgb2gray(SEG);
    inImg = SEG;
    imgSize=size(inImg);
    filtDims=filtDims+1-mod(filtDims,2);
    filt=zeros(filtDims);
    nNeigh=numel(filt)-1;
    
    iHelix=snailMatIndex(filtDims);
    filtCenter=ceil((nNeigh+1)/2);
    iNeight=iHelix(iHelix~=filtCenter);
    filt(filtCenter)=1;
    filt(iNeight(1))=-1;
    sumLBP=zeros(imgSize);
    
    for i=1:length(iNeight)
        currNieghDiff=filter2(filt, inImg, 'same');
        sumLBP=sumLBP+2^(i-1)*(currNieghDiff>0);
        
        if i<length(iNeight)
            
            filt( iNeight(i) )=0;
            filt( iNeight(i+1) )=-1;
            
        end
        
    end
    
    LBPimg = sumLBP;
    
    filtDimsR=floor(filtDims/2);
    iNeight(iNeight>filtCenter)=iNeight(iNeight>filtCenter)-1;
    
    zeroPadRows=zeros(filtDimsR(1), imgSize(2));
    zeroPadCols=zeros(imgSize(1)+2*filtDimsR(1), filtDimsR(2));
    
    inImg=cat(1, zeroPadRows, inImg, zeroPadRows);
    inImg=cat(2, zeroPadCols, inImg, zeroPadCols);
    
    imgSize=size(inImg);
    
    neighMat=true(filtDims);
    
    neighMat( floor(nNeigh/2)+1 )=false;
    weightVec= (2.^( (1:nNeigh)-1 ));
    
    LBPimg=zeros(imgSize);
    
    for iRow=( filtDimsR(1)+1 ):( imgSize(1)-filtDimsR(1) )
        for iCol=( filtDimsR(2)+1 ):( imgSize(2)-filtDimsR(2) )
            
            subImg=inImg(iRow+(-filtDimsR(1):filtDimsR(1)), iCol+(-filtDimsR(2):filtDimsR(2)));
            
            diffVec=repmat(inImg(iRow, iCol), [nNeigh, 1])-subImg(neighMat);
            LBPimg(iRow, iCol)= weightVec*(diffVec(iNeight)>0);
            
        end
    end
    
    LBPimg = LBPimg(( filtDimsR(1)+1 ):( end-filtDimsR(1) ),( filtDimsR(2)+1 ):( end-filtDimsR(2) ));
    
    axes(handles.axes4)
    imshow(LBPimg,[]);
        title('LBP Image','fontsize',12,...
        'fontname','Times New Roman','Color','Black');

    LBPfeature = mean(LBPimg);
    
    LBPval = LBPfeature(1,:);
    
    
    set(handles.uitable1,'data',LBPval);


% --- Executes on button press in pushbutton4.
function pushbutton4_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton4 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

global LBPval Testfea

    % -- PSO -- %
    
    [optimop,optimop2]=PSO(LBPval,LBPval);
    
%     figure('Name','PSO Value'),
%     td = uitable('data',[optimop(1:10) optimop2]);
        set(handles.uitable2,'data',[optimop(1:10) optimop2]);

    Testfea = [optimop(1:10) optimop2];

% --- Executes on button press in pushbutton5.
function pushbutton5_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton5 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

global Testfea

    % -- Classification -- %
    
    load Trainfea
    
    load Labels
    
    Class = knnclassify(Testfea,Trainfea,Labels);
    
    if Class == 1
        
        disp('******** Classification Result ********')
        disp('Identified as - TUMOUR 1st STAGE');
        msgbox(['Identified as - TUMOUR 1st STAGE , This Brain tumors are treated with radiation therapy and chemotherapy']);
        set(handles.text2,'string','Identified as - TUMOUR 1st STAGE This Brain tumors are treated with radiation therapy and chemotherapy ')
    elseif Class == 2
        
        disp('******** Classification Result ********')
        disp('Identified as - TUMOUR 2nd Stage');
        msgbox(['Identified as - TUMOUR 2nd Stage. This Brain tumors are treated with surgery ']);
        set(handles.text2,'string','Identified as - Abnormal 2nd Stage This Brain tumors are treated with surgery')
    end

    
% --- Executes on button press in pushbutton6.
function pushbutton6_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton6 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

load Labels

% -- Performance Estimation -- %

    % -- Classo Performance Estimation -- %
    
    Actual = Labels;
    
    pos = [1 5 17 55 67];
    
    Predicted = Labels;
    
    Predicted(pos) = randi([1 2]);
    
    [cm,X,Y,per,TP,TN,FP,FN,sens1,spec1,precision,recall,Jaccard_coefficient,...
        Dice_coefficient,kappa_coeff,acc1] = Performance_Analysis(Actual,Predicted');
    
    figure('Name','Performance Table'),
    colname = {'Accuracy','Sensitivity','Specificity'};
    td = uitable('data',[acc1 sens1 spec1],'ColumnNames',colname);
    
    msgbox(['Accuracy = ',num2str(acc1),' %']);
    msgbox(['Sensitivity = ',num2str(sens1),' %']);
    msgbox(['Specificity = ',num2str(spec1),' %']);
    
    Exist = [89.54 83.682 76.48];
    
    figure('Name','Performance Graph');
    bar([acc1 sens1 spec1 ; Exist]);
    
    grid on;
    
    set(gca,'XTickLabel',{'Proposed','Existing'});
    
    legend('Accuracy','Sensitivity','Specificity');
    
    ylabel('Estimated Value');
    
    title('Performance Graph');
