function varargout = DBJH_CAL(varargin)
% CYLINDER_PORE MATLAB code for DBJH_CAL.fig
%      CYLINDER_PORE, by itself, creates a new CYLINDER_PORE or raises the existing
%      singleton*.
%
%      H = DBJH-CAL returns the handle to a new DBJH-CAL or the handle to
%      the existing singleton*.
%
%      CYLINDER_PORE('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in CYLINDER_PORE.M with the given input arguments.
%
%      CYLINDER_PORE('Property','Value',...) creates a new CYLINDER_PORE or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before cylinder_pore_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to cylinder_pore_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help cylinder_pore

% Last Modified by GUIDE v2.5 06-Jul-2020 00:42:06

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @DBJH_CAL_OpeningFcn, ...
                   'gui_OutputFcn',  @DBJH_CAL_OutputFcn, ...
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


% --- Executes just before cylinder_pore is made visible.
function DBJH_CAL_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to cylinder_pore (see VARARGIN)

% Choose default command line output for cylinder_pore
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes cylinder_pore wait for user response (see UIRESUME)
% uiwait(handles.figure_dbjhpsdcal);
setappdata(handles.figure_dbjhpsdcal,'pore_size',8);
setappdata(handles.figure_dbjhpsdcal,'isotherm',[]);
setappdata(handles.figure_dbjhpsdcal,'isofilename',0);
setappdata(handles.figure_dbjhpsdcal,'min_psd',0);
setappdata(handles.figure_dbjhpsdcal,'max_psd',50);
setappdata(handles.figure_dbjhpsdcal,'integ_psd',1);
%origin set
set(handles.pan_psdcal,'visible','on');
set(handles.cal_cylind,'Value',1);
set(handles.cal_slit,'Value',1);
setappdata(handles.figure_dbjhpsdcal,'tpcyld',1);
setappdata(handles.figure_dbjhpsdcal,'tpslit',1);
setappdata(handles.figure_dbjhpsdcal,'isotherm_fa',[]);
setappdata(handles.figure_dbjhpsdcal,'isotherm_fd',[]);
setappdata(handles.figure_dbjhpsdcal,'Iso_fa',[]);
setappdata(handles.figure_dbjhpsdcal,'Iso_fd',[]);
setappdata(handles.figure_dbjhpsdcal,'Iso_ca',[]);
setappdata(handles.figure_dbjhpsdcal,'Iso_cd',[]);
setappdata(handles.figure_dbjhpsdcal,'Iso_sa',[]);
setappdata(handles.figure_dbjhpsdcal,'Iso_sd',[]);
setappdata(handles.figure_dbjhpsdcal,'psd_s',[]);
setappdata(handles.figure_dbjhpsdcal,'psd_c',[]);

% --- Outputs from this function are returned to the command line.
function varargout = DBJH_CAL_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;

% --- Executes on button press in calculation.
function calculation_Callback(hObject, eventdata, handles)
% hObject    handle to calculation (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
isotherm=getappdata(handles.figure_dbjhpsdcal,'isotherm');
min_psd=getappdata(handles.figure_dbjhpsdcal,'min_psd');
max_psd=getappdata(handles.figure_dbjhpsdcal,'max_psd');
integ_psd=getappdata(handles.figure_dbjhpsdcal,'integ_psd');
% pore_type=getappdata(handles.figure_dbjhpsdcal,'cylindpore_type')+2*getappdata(handles.figure_dbjhpsdcal,'inkbotpore_type')+5*getappdata(handles.figure_dbjhpsdcal,'slitpore_type');
n_point=round((max_psd-min_psd)/integ_psd);

tleng=length(isotherm);
[~,maxidn]=max(isotherm(:,1));
isotherm_fa=zeros(n_point-1,2);
isotherm_fd=isotherm_fa;

isotherm_a=isotherm(1:maxidn,:);
isotherm_d=isotherm(maxidn:tleng,:);
r_ca_p=isotherm_fa;
r_sa_p=isotherm_fa;
r_sd_p=isotherm_fa;
r_cd_p=isotherm_fa;

% 构建圆柱孔和狭缝孔  注意吸附脱附不一致
for nn=2:n_point
    ppa=(nn-0.5)/n_point;
    ppd=(ppa+0.2)/1.2; %让脱附支的最小压力从0.16开始，符合实际
    tma=-0.37*nthroot((5/log(ppa)),3);
    tmd=-0.37*nthroot((5/log(ppd)),3);
    double_rka=-0.9533/log(ppa);% return natural logarithm ln(x) 
    single_rka=double_rka/2;
    double_rkd=-0.9533/log(ppd);% return natural logarithm ln(x) 
    single_rkd=double_rkd/2;
    % 等温线插值
    isotherm_fa(nn-1,2)=mediff(isotherm_a,ppa);
    isotherm_fd(nn-1,2)=mediff(isotherm_d,ppd);
    isotherm_fa(nn-1,1)=ppa;
    isotherm_fd(nn-1,1)=ppd;
    
    r_ca_p(nn-1,2)=single_rka+tma;
    r_sa_p(nn-1,2)=tma;
    r_cd_p(nn-1,2)=double_rkd+tmd;
    r_sd_p(nn-1,2)=single_rkd+tmd;
    r_ca_p(nn-1,1)=ppa;
    r_sa_p(nn-1,1)=ppa;
    r_cd_p(nn-1,1)=ppd;
    r_sd_p(nn-1,1)=ppd;
    %r_ca_p/r_sa_p/r_cd_p/r_sd_p就是构建好的临界孔径/压力对应曲线
end

setappdata(handles.figure_dbjhpsdcal,'isotherm_fa',isotherm_fa);
setappdata(handles.figure_dbjhpsdcal,'isotherm_fd',isotherm_fd);
% figure,plot(isotherm_fa(:,1),isotherm_fa(:,2))
% figure,plot(isotherm_fd(:,1),isotherm_fd(:,2))

dr_dp_ca=gradient(r_ca_p(:,2),r_ca_p(:,1));
dr_dp_cd=gradient(r_cd_p(:,2),r_cd_p(:,1));
dr_dp_sa=gradient(r_sa_p(:,2),r_sa_p(:,1));
dr_dp_sd=gradient(r_sd_p(:,2),r_sd_p(:,1));
%dr/dp数据也构建完成 单列向量

disoa=gradient(isotherm_fa(:,2),isotherm_fa(:,1));
disod=gradient(isotherm_fd(:,2),isotherm_fd(:,1));
%diso/dp也构建完成 单列向量

poretp=getappdata(handles.figure_dbjhpsdcal,'tpcyld')+2*getappdata(handles.figure_dbjhpsdcal,'tpslit');% 1--cylind;2--slit;3--both
dvsa_dr(:,1)=r_sa_p(:,2);
dvsa_dr(:,2)=disoa./dr_dp_sa;
if poretp==1
    dvca_dr(:,1)=r_ca_p(:,2);
    dvca_dr(:,2)=disoa./dr_dp_ca;
    dvcd_dr(:,1)=r_cd_p(:,2);
    dvcd_dr(:,2)=disod./dr_dp_cd;
    dvsd_dr(:,1)=r_sd_p(:,2);
    dvsd_dr(:,2)=0;
    dvsa_dr(:,1)=r_sa_p(:,2);
    dvsa_dr(:,2)=0;
elseif poretp==2
    dvsd_dr(:,1)=r_sd_p(:,2);
    dvsd_dr(:,2)=disod./dr_dp_sd;
    dvca_dr(:,1)=r_ca_p(:,2);
    dvca_dr(:,2)=0;
    dvcd_dr(:,1)=r_cd_p(:,2);
    dvcd_dr(:,2)=0;
elseif poretp==3
    %构建0初值情况加入迭代计算
    dvsa_dr(:,2)=disoa./dr_dp_ca;
    dvca_dr(:,1)=r_ca_p(:,2);
    dvca_dr(:,2)=0;
    dvcd_dr2(:,1)=r_ca_p(:,2);
    dvcd_dr2(:,2)=0;
    sctp=1;%sctp--1,sa起算；sctp--2,ca起算
    erros=500;
    
    while (erros>0.1)
        [dvsd_dr,dvcd_dr,dvsa_dr,dvca_dr]=itera_elem(disoa,disod,dvsa_dr,dvca_dr,dr_dp_sa,dr_dp_ca,dr_dp_sd,dr_dp_cd,r_ca_p,r_cd_p,r_sa_p,r_sd_p,sctp);
        dvsa_dr2=multidiff(dvsd_dr,r_sa_p(:,2));
        erros=sum((dvcd_dr2(:,2)-dvcd_dr(:,2)).*(dvcd_dr2(:,2)-dvcd_dr(:,2)));
        dvsa_dr=dvsa_dr2;
%         dvca_dr2=multidiff(dvcd_dr,r_ca_p(:,2));
        dvcd_dr2=dvcd_dr;
    end
end
dvs_dr=[dvsa_dr;dvsd_dr];
dvs_drtt=sortrows(dvs_dr);
dvn=length(dvs_drtt);
idx=(1:2:dvn);
dvc_dr=[dvca_dr;dvcd_dr];
dvc_drtt=sortrows(dvc_dr);

dvs_drt=dvs_drtt(idx,:);
dvc_drt=dvc_drtt(idx,:);

%反算等温线
% dppa=1/n_point;
% dppd=1/(1.2*n_point); 
dIsoa(:,1)=isotherm_fa(:,1);
dIsoa(:,2)=dvsa_dr(:,2).*dr_dp_sa+dvca_dr(:,2).*dr_dp_ca;
dIsod(:,1)=isotherm_fd(:,1);
dIsod(:,2)=dvsd_dr(:,2).*dr_dp_sd+dvcd_dr(:,2).*dr_dp_cd;

dIsoca(:,1)=isotherm_fa(:,1);
dIsoca(:,2)=dvca_dr(:,2).*dr_dp_ca;
dIsocd(:,1)=isotherm_fd(:,1);
dIsocd(:,2)=dvcd_dr(:,2).*dr_dp_cd;
dIsosa(:,1)=isotherm_fa(:,1);
dIsosa(:,2)=dvsa_dr(:,2).*dr_dp_sa;
dIsosd(:,1)=isotherm_fd(:,1);
dIsosd(:,2)=dvsd_dr(:,2).*dr_dp_sd;

Isoa=integ(dIsoa);
Isoa(:,2)=Isoa(:,2)+isotherm_fa(1,2);
Isod=integ(dIsod);
Isod(:,2)=Isod(:,2)+isotherm_fd(1,2);
Isoca=integ(dIsoca);
Isoca(:,2)=Isoca(:,2)+0.6*isotherm_fa(1,2);
Isosa=integ(dIsosa);
Isosa(:,2)=Isosa(:,2)+0.4*isotherm_fa(1,2);
Isocd=integ(dIsocd);
Isocd(:,2)=Isocd(:,2)+0.6*isotherm_fd(1,2);
Isosd=integ(dIsosd);
Isosd(:,2)=Isosd(:,2)+0.4*isotherm_fd(1,2);

setappdata(handles.figure_dbjhpsdcal,'psd_c',dvc_drt);
setappdata(handles.figure_dbjhpsdcal,'psd_s',dvs_drt);
setappdata(handles.figure_dbjhpsdcal,'Iso_fa',Isoa);
setappdata(handles.figure_dbjhpsdcal,'Iso_fd',Isod);
setappdata(handles.figure_dbjhpsdcal,'Iso_ca',Isoca);
setappdata(handles.figure_dbjhpsdcal,'Iso_cd',Isocd);
setappdata(handles.figure_dbjhpsdcal,'Iso_sa',Isosa);
setappdata(handles.figure_dbjhpsdcal,'Iso_sd',Isosd);

axes(handles.axes_psd_slit)
plot(dvs_drt(:,1)*2,dvs_drt(:,2)/2);
xlim([min_psd max_psd]);
axes(handles.axes_psd_cyld)
plot(dvc_drt(:,1)*2,dvc_drt(:,2)/2);
xlim([min_psd max_psd]);

axes(handles.axes_recal_iso)
plot(Isoa(:,1),Isoa(:,2));
hold on
plot(Isod(:,1),Isod(:,2));
hold on

axes(handles.axes_multiPSD)
plot(Isoca(:,1),Isoca(:,2)),title('圆柱孔/狭缝孔各自等温线');
hold on
plot(Isocd(:,1),Isocd(:,2));
hold on
plot(Isosa(:,1),Isosa(:,2));
hold on
plot(Isosd(:,1),Isosd(:,2));
legend('圆柱孔吸附支','圆柱孔脱附支','狭缝孔吸附支','狭缝孔脱附支');

% --- Executes on button press in data_out.
function data_out_Callback(hObject, eventdata, handles)
% hObject    handle to data_out (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
isotherm_fa=getappdata(handles.figure_dbjhpsdcal,'isotherm_fa');
isotherm_fd=getappdata(handles.figure_dbjhpsdcal,'isotherm_fd');
psd_c=getappdata(handles.figure_dbjhpsdcal,'psd_c');
psd_s=getappdata(handles.figure_dbjhpsdcal,'psd_s');
Isoa=getappdata(handles.figure_dbjhpsdcal,'Iso_fa');
Isod=getappdata(handles.figure_dbjhpsdcal,'Iso_fd');
Isoca=getappdata(handles.figure_dbjhpsdcal,'Iso_ca');
Isocd=getappdata(handles.figure_dbjhpsdcal,'Iso_cd');
Isosa=getappdata(handles.figure_dbjhpsdcal,'Iso_sa');
Isosd=getappdata(handles.figure_dbjhpsdcal,'Iso_sd');
ins=length(isotherm_fa);

filespec_user = [pwd '\DBJH模型计算结果.xlsx'];
try
    Excel = actxGetRunningServer('Excel.Application');
catch
    Excel = actxserver('Excel.Application'); 
end;
Excel.Visible = 0;    % set(Excel, 'Visible', 1); 
if ~exist(filespec_user,'file'); 
    Workbook = Excel.Workbooks.Add;
    Workbook.SaveAs(filespec_user);
else
    Workbook = Excel.Workbooks.Open(filespec_user);
end
% Sheets = Excel.ActiveWorkbook.Sheets;    % Sheets = Workbook.Sheets;
% Sheet1 = Sheets.Item(1);
% Sheet1.Activate;
Workbook.Save   % 保存文档
Workbook.Close

A = xlsread('DBJH模型计算结果.xlsx');
[row,~]=size(A);
if row==0
    stringa=strcat('A',num2str(1),':T',num2str(1));
%     xlswrite('吸附等温线分类.xlsx','等温线类型','Sheet1','A1');
    result={'实测isop','实测isov','实测dsop','实测dsov','预测isop','预测isov','预测dsop','预测dsov','psd_cr','psd_cv','psd_sr','psd_sv','iso_cap','iso_cav','iso_cdp','iso_cdv','iso_sap','iso_sav','iso_sdp','iso_sdv'};
    sheetname='双孔模型测试';
    xlswrite('DBJH模型计算结果.xlsx',result,sheetname,stringa);
end

stringa=strcat('A',num2str(2),':T',num2str(ins+1));

result=[isotherm_fa,isotherm_fd,Isoa,Isod,psd_c,psd_s,Isoca,Isocd,Isosa,Isosd];
sheetname='双孔模型测试';
xlswrite('DBJH模型计算结果.xlsx',result,sheetname,stringa);

%给定dvsa_dr，求迭代后的dvc_dr和dvs_dr，dvsa_dr与dvsd_dr本质上是一样的，不同根本在于坐标点
function [dvsd_dr,dvcd_dr,dvsa_dr,dvca_dr]=itera_elem(disoa,disod,dvsa_dr,dvca_dr,dr_sa_dp,dr_ca_dp,dr_sd_dp,dr_cd_dp,r_ca_p,r_cd_p,r_sa_p,r_sd_p,sctp)%迭代基元
if sctp==1 %sa初始化
    dvca_dr(:,1)=r_ca_p(:,2);
    tpca=(disoa-dvsa_dr(:,2).*dr_sa_dp)./dr_ca_dp;
    tpca(tpca<0)=0;
    dvca_dr(:,2)=tpca;
    %求 dvcd_dr
    dvcd_dr=multidiff(dvca_dr,r_cd_p(:,2));
    %求 dvsd_dr
    dvsd_dr(:,1)=r_sd_p(:,2);
    tpsd=(disod-dvcd_dr(:,2).*dr_cd_dp)./dr_sd_dp;
    tpsd(tpsd<0)=0;
    dvsd_dr(:,2)=tpsd;
else %ca初始化
    dvsa_dr(:,1)=r_sa_p(:,2);
    tpsa=(disoa-dvca_dr(:,2).*dr_ca_dp)./dr_sa_dp;
    tpsa(tpsa<0)=0;
    dvsa_dr(:,2)=tpsa;
    %求 dvsd_dr
    dvsd_dr=multidiff(dvsa_dr,r_sd_p(:,2));
    %求 dvcd_dr
    dvcd_dr(:,1)=r_cd_p(:,2);
    tpcd=(disod-dvsd_dr(:,2).*dr_sd_dp)./dr_cd_dp;
    tpcd(tpcd<0)=0;
    dvcd_dr(:,2)=tpcd;
end

function integva=integ(vector)%给定一条向量，第一列为x，第二列为y，按x进行y的积分。
integva=vector;
nl=length(vector);
intx=[vector(1,1);vector(:,1);vector(nl,1)];
for ii=1:length(vector)
    dy=0;
    for ij=1:ii
        dx=(intx(ij+2)-intx(ij))/2;
        dy=dy+dx*vector(ij,2);
    end
    integva(ii,2)=dy;
end

function filansw=filtp(vector) %取正插值
    len=length(vector);
    pidx=find(vector>0);%大于0的索引
    plen=length(pidx);%索引向量的长度
    ve=vector;
    for lg=1:len
        if(vector(lg)<=0)
            upv=find(lg<pidx);
            dwv=find(lg>pidx);%索引值的位置索引
            if isempty(upv)
                upv=plen;
            end
            if isempty(dwv)
                dwv=1;
            end
            dwidx=length(dwv);
            upvalue=pidx(upv(1));
            dwvalue=pidx(dwv(dwidx));
            ve(lg)=vector(upvalue)+(vector(upvalue)-vector(dwvalue))*(upvalue-lg)/(upvalue-dwvalue);
        end
    end
    filansw=ve; 

function diffansw=mediff1(vector,xvalue) %单插值函数，vector 第一列为x,第二列为y
    upvectxo=find(vector(:,1)>xvalue);
    dwvectxo=find(vector(:,1)<xvalue);
    [sizedw,~]=size(dwvectxo);
    [sizeup,~]=size(upvectxo);

    % 考虑不在中间的情况，进行补充
    if ~isempty(dwvectxo)
        dwn=dwvectxo(sizedw);
        rr=0;
    else
        dwn=upvectxo(1);
        rr=1;
    end
    if ~isempty(upvectxo)
        upn=upvectxo(1+rr);
    else
        upn=dwvectxo(sizedw);
        dwn=dwvectxo(sizedw-1);
    end

    if(upn<dwn) %反序
        upidx=sizeup;
        fupidx=sizeup-1;
        dwidx=1;
        fdwidx=2;
    else %正序
        upidx=1;
        fupidx=2;
        dwidx=sizedw;
        fdwidx=sizedw-1;
    end

    if ~isempty(upvectxo)
        upvectx=upvectxo(upidx);
    else
        upvectx=dwvectxo(fdwidx);
    end
    if ~isempty(dwvectxo)
        dwvectx=dwvectxo(dwidx);
    else
        dwvectx=upvectxo(fupidx);
    end
    diffansw=((vector(upvectx,2)-vector(dwvectx,2))*(xvalue-vector(dwvectx,1)))/(vector(upvectx,1)-vector(dwvectx,1))+vector(dwvectx,2);

function multasw=multidiff(vector,xvect)
xvlen=length(xvect);
multasw=zeros(xvlen,2);
multasw(:,1)=xvect;
for mi=1:xvlen
    multasw(mi,2)=mediff(vector,xvect(mi));
end
        
function medasw=mediff(vector,xvalue)% 造点法插值,可能更可好
    nl=length(vector);
    seriord=vector(nl)-vector(1);%>0，正序，<0,倒序，==0，单数
    xmax=max(vector(:,1));
    xmin=min(vector(:,1));
    xtop=max(xvalue,xmax)+abs(seriord/nl)*0.2;
    xbot=min(xvalue,xmin)*0.8;
    if seriord>0
        topx=vector(nl,1)-vector(nl-1,1);
        topy=vector(nl,2)-vector(nl-1,2);
        topy(topx==0)=vector(nl,2)-vector(nl-2,2);
        topx(topx==0)=vector(nl,1)-vector(nl-2,1);
        
        botx=vector(2,1)-vector(1,1);
        boty=vector(2,2)-vector(1,2);
        boty(botx==0)=vector(3,2)-vector(1,2);
        botx(botx==0)=vector(3,1)-vector(1,1);
        
        vtop=(xtop-vector(nl,1))*topy/topx+vector(nl,2);
        vbot=(xbot-vector(1,1))*boty/botx+vector(1,2);
        vectorf=[[xbot,vbot];vector;[xtop,vtop]];
    elseif seriord<0
        botx=vector(nl,1)-vector(nl-1,1);
        boty=vector(nl,2)-vector(nl-1,2);
        boty(botx==0)=vector(nl,2)-vector(nl-2,2);
        botx(botx==0)=vector(nl,1)-vector(nl-2,1);
        
        topx=vector(2,1)-vector(1,1);
        topy=vector(2,2)-vector(1,2);
        topy(topx==0)=vector(3,2)-vector(1,2);
        topx(topx==0)=vector(3,1)-vector(1,1);
        
        vbot=(xbot-vector(nl,1))*boty/botx+vector(nl,2);
        vtop=(xtop-vector(1,1))*topy/topx+vector(1,2);
        vectorf=[[xtop,vtop];vector;[xbot,vbot]];
    end
    upidx=find(vectorf(:,1)>xvalue);
    dwidx=find(vectorf(:,1)<xvalue);
    uplg=length(upidx);
    dwlg=length(dwidx);
    if seriord>0
        xup=vectorf(upidx(1),1);
        yup=vectorf(upidx(1),2);
        xdw=vectorf(dwidx(dwlg),1);
        ydw=vectorf(dwidx(dwlg),2);
    else
        xup=vectorf(upidx(uplg),1);
        yup=vectorf(upidx(uplg),2);
        xdw=vectorf(dwidx(1),1);
        ydw=vectorf(dwidx(1),2);
    end
    medasw=(xvalue-xup)*(yup-ydw)/(xup-xdw)+yup;
    
% --- Executes on button press in cal_siglepore.
function cal_siglepore_Callback(hObject, eventdata, handles)
% hObject    handle to cal_siglepore (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
set(handles.pan_nitro_para,'visible','off');
set(handles.pan_pore_size,'visible','off');

% --- Executes on button press in pore_size_dis.
function pore_size_dis_Callback(hObject, eventdata, handles)
% hObject    handle to pore_size_dis (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
[filelist,filepath]=uigetfile('*.txt;*.mat;*.xlsx','Pick an isotherm');
if isequal(filelist,0)||isequal(filepath,0)   % filelist have  '\'
    return;
end

fpath=[filepath filelist];
filename=filelist(1,:);
setappdata(handles.figure_dbjhpsdcal,'filename',filename);

aa=strfind(fpath,'xlsx');
exceltag(~isempty(aa))=1;
setappdata(handles.figure_dbjhpsdcal,'exceltag',exceltag);

if getappdata(handles.figure_dbjhpsdcal,'exceltag')==1
    psd_data=xlsread(fpath,'吸附分支');
    [row,~]=size(psd_data);
    psd_data2(:,1)=psd_data(2:row,3);
    psd_data2(:,2)=psd_data(2:row,6);
else
    psd_data2=load(fpath); 
end

[ndxr,~]=size(psd_data2);  %孔径分布的孔径数目
ndxr(ndxr==0)=1;
if ndxr>1
    axes(handles.axes_iso);
    plot(psd_data2(:,1),psd_data2(:,2));
end
setappdata(handles.figure_dbjhpsdcal,'pore_size_dis',psd_data2);


function cylind_diam_Callback(hObject, eventdata, handles)
% hObject    handle to cylind_diam (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% Hints: get(hObject,'String') returns contents of cylind_diam as text
%        str2double(get(hObject,'String')) returns contents of cylind_diam as a double
val=str2double(get(hObject,'String'));
setappdata(handles.figure_dbjhpsdcal,'cylind_diam',val);

% --- Executes during object creation, after setting all properties.
function cylind_diam_CreateFcn(hObject, eventdata, handles)
% hObject    handle to cylind_diam (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called
% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



% --- Executes on button press in cylind_pore.
function cylind_pore_Callback(hObject, eventdata, handles)
% hObject    handle to cylind_pore (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
val=get(handles.slit_pore,'Value');
setappdata(handles.figure_dbjhpsdcal,'cylind_pore',val);

% Hint: get(hObject,'Value') returns toggle state of cylind_pore
set(handles.text_cylind_diam,'enable','on');
set(handles.cylind_diam,'enable','on');
set(handles.text_cylind_unit,'enable','on');
set(handles.text_cylind_length,'enable','on');
set(handles.ratio_cylind_length,'enable','on');
set(handles.para_confirm,'enable','on');

% --- Executes on button press in slit_pore.
function slit_pore_Callback(hObject, eventdata, handles)
% hObject    handle to slit_pore (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
val=get(handles.slit_pore,'Value');
setappdata(handles.figure_dbjhpsdcal,'slit_pore',val);
% Hint: get(hObject,'Value') returns toggle state of slit_pore
set(handles.text_width_slit,'enable','on');
set(handles.slit_width,'enable','on');
set(handles.slit_width_unit,'enable','on');
set(handles.text_slit_length,'enable','on');
set(handles.ratio_slit_length,'enable','on');
set(handles.text_slit_height,'enable','on');
set(handles.ratio_slit_height,'enable','on');
set(handles.para_confirm,'enable','on');


% --- Executes on button press in radio_cylind.
function radio_cylind_Callback(hObject, eventdata, handles)
% hObject    handle to radio_cylind (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of radio_cylind
set(handles.radio_inkbot,'Value',0);
set(handles.radio_slit,'Value',0);
val=1*get(hObject,'Value');
setappdata(handles.figure_dbjhpsdcal,'cylindpore_type',val);

% --- Executes on button press in radio_slit.
function radio_slit_Callback(hObject, eventdata, handles)
% hObject    handle to radio_slit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of radio_slit
set(handles.radio_inkbot,'Value',0);
set(handles.radio_cylind,'Value',0);
val=get(hObject,'Value');
setappdata(handles.figure_dbjhpsdcal,'slitpore_type',val);

function pore_size_Callback(hObject, eventdata, handles)
% hObject    handle to pore_size (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of pore_size as text
%        str2double(get(hObject,'String')) returns contents of pore_size as a double
val=str2double(get(hObject,'String'));
setappdata(handles.figure_dbjhpsdcal,'pore_size',val);

% --- Executes during object creation, after setting all properties.
function pore_size_CreateFcn(hObject, eventdata, handles)
% hObject    handle to pore_size (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% --- Executes on button press in iso_input.
function iso_input_Callback(hObject, eventdata, handles)
% hObject    handle to iso_input (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global filepath filelist;
[filelist,filepath]=uigetfile('*.txt;*.mat;*.xlsx','Pick an isotherm');
if isequal(filelist,0)||isequal(filepath,0)   % filelist have  '\'
    return;
end
fpath=[filepath filelist];
filename=filelist(1,:);
setappdata(handles.figure_dbjhpsdcal,'isofilename',filename);
set(handles.iso_address,'String',fpath);
aa=strfind(fpath,'xlsx');
exceltag(~isempty(aa))=1;
setappdata(handles.figure_dbjhpsdcal,'exceltag',exceltag);

if getappdata(handles.figure_dbjhpsdcal,'exceltag')==1
    iso_data=xlsread(fpath,'吸附等温线');
    [row,~]=size(iso_data);
    iso_data2(:,1)=iso_data(2:row,1);
    iso_data2(:,2)=iso_data(2:row,2);
else
    iso_data2=load(fpath); 
end

[ndxr,~]=size(iso_data2);  %孔径分布的孔径数目
ndxr(ndxr==0)=1;
if ndxr>1
    axes(handles.axes_iso);
    plot(iso_data2(:,1),iso_data2(:,2));
end
setappdata(handles.figure_dbjhpsdcal,'isotherm',iso_data2);

% --- Executes during object creation, after setting all properties.
function iso_address_CreateFcn(hObject, eventdata, handles)
% hObject    handle to iso_address (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% --- Executes on button press in cal_psd.
function cal_psd_Callback(hObject, eventdata, handles)
% hObject    handle to cal_psd (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
tpslit=getappdata(handles.figure_dbjhpsdcal,'tpslit');
tpcyld=getappdata(handles.figure_dbjhpsdcal,'tpcyld');
if (tpslit+tpcyld==0)
    uiwait(msgbox('您未选择孔型，请重新选择','error','modal'));
    return
else
    calculation_Callback(hObject, eventdata, handles);
end


% --- Executes on slider movement.
function slider_min_pore_Callback(hObject, eventdata, handles)
% hObject    handle to slider_min_pore (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider
val=get(hObject,'Value');
max_psd=getappdata(handles.figure_dbjhpsdcal,'max_psd');
val(val<max_psd)=max_psd;
setappdata(handles.figure_dbjhpsdcal,'min_psd',val);
set(handles.dw_psd,'String',num2str(val));

% --- Executes during object creation, after setting all properties.
function slider_min_pore_CreateFcn(hObject, eventdata, handles)
% hObject    handle to slider_min_pore (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end

% --- Executes on slider movement.
function slider_max_pore_Callback(hObject, eventdata, handles)
% hObject    handle to slider_max_pore (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider
val=get(hObject,'Value');
min_psd=getappdata(handles.figure_dbjhpsdcal,'min_psd');
val(val<min_psd)=min_psd;
set(handles.slider_max_pore,'Value',val);
setappdata(handles.figure_dbjhpsdcal,'max_psd',val);
set(handles.up_psd,'String',num2str(val));

% --- Executes during object creation, after setting all properties.
function slider_max_pore_CreateFcn(hObject, eventdata, handles)
% hObject    handle to slider_max_pore (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end

% --- Executes on slider movement.
function slider_integ_Callback(hObject, eventdata, handles)
% hObject    handle to slider_integ (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider
val=get(hObject,'Value');
setappdata(handles.figure_dbjhpsdcal,'integ_psd',val);
set(handles.intg_psd,'String',num2str(val));

% --- Executes during object creation, after setting all properties.
function slider_integ_CreateFcn(hObject, eventdata, handles)
% hObject    handle to slider_integ (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.

if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end

function dw_psd_Callback(hObject, eventdata, handles)
% hObject    handle to dw_psd (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of dw_psd as text
%        str2double(get(hObject,'String')) returns contents of dw_psd as a double
val=str2double(get(hObject,'String'));
setappdata(handles.figure_dbjhpsdcal,'min_psd',val);
set(handles.slider_min_pore,'Value',val);

% --- Executes during object creation, after setting all properties.
function dw_psd_CreateFcn(hObject, eventdata, handles)
% hObject    handle to dw_psd (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function up_psd_Callback(hObject, eventdata, handles)
% hObject    handle to up_psd (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of up_psd as text
%        str2double(get(hObject,'String')) returns contents of up_psd as a double
val=str2double(get(hObject,'String'));
setappdata(handles.figure_dbjhpsdcal,'max_psd',val);
set(handles.slider_max_pore,'Value',val);

% --- Executes during object creation, after setting all properties.
function up_psd_CreateFcn(hObject, eventdata, handles)
% hObject    handle to up_psd (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function intg_psd_Callback(hObject, eventdata, handles)
% hObject    handle to intg_psd (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of intg_psd as text
%        str2double(get(hObject,'String')) returns contents of intg_psd as a double
val=str2double(get(hObject,'String'));
setappdata(handles.figure_dbjhpsdcal,'integ_psd',val);
set(handles.slider_integ,'Value',val);

% --- Executes during object creation, after setting all properties.
function intg_psd_CreateFcn(hObject, eventdata, handles)
% hObject    handle to intg_psd (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% --- Executes on button press in cal_cylind.
function cal_cylind_Callback(hObject, eventdata, handles)
% hObject    handle to cal_cylind (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of cal_cylind
val=get(hObject,'Value');
setappdata(handles.figure_dbjhpsdcal,'tpcyld',val);


% --- Executes on button press in cal_slit.
function cal_slit_Callback(hObject, eventdata, handles)
% hObject    handle to cal_slit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of cal_slit
val=get(hObject,'Value');
setappdata(handles.figure_dbjhpsdcal,'tpslit',val);
