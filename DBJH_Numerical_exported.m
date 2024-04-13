classdef DBJH_Numerical_exported < matlab.apps.AppBase

    % Properties that correspond to app components
    properties (Access = public)
        DBJH_FIG                      matlab.ui.Figure
        OperatePanel                  matlab.ui.container.Panel
        OpenFileButton                matlab.ui.control.Button
        CalculateButton               matlab.ui.control.Button
        ExportDataButton              matlab.ui.control.Button
        FileAddress_TXT               matlab.ui.control.Label
        FileAddress                   matlab.ui.control.EditField
        DisplayPanel                  matlab.ui.container.Panel
        UISO                          matlab.ui.control.UIAxes
        UIBJHPSD                      matlab.ui.control.UIAxes
        UISlope                       matlab.ui.control.UIAxes
        UIPSD                         matlab.ui.control.UIAxes
        SwitchYaxisformButtonGroup    matlab.ui.container.ButtonGroup
        NormalY                       matlab.ui.control.StateButton
        LogY                          matlab.ui.control.StateButton
        Material_Property_Panal       matlab.ui.container.Panel
        MolecWtEditFieldLabel         matlab.ui.control.Label
        MolecWt_val                   matlab.ui.control.EditField
        LiquidDensityEditFieldLabel   matlab.ui.control.Label
        LiquidDensity_val             matlab.ui.control.EditField
        CrossSectionEditFieldLabel    matlab.ui.control.Label
        CrossSection_val              matlab.ui.control.EditField
        KerogenDensityEditFieldLabel  matlab.ui.control.Label
        KerogenDensity_val            matlab.ui.control.EditField
        KerogenDensity_Unit           matlab.ui.control.Label
        CrossSection_Unit             matlab.ui.control.Label
        LiquidDensity_Unit            matlab.ui.control.Label
        MolecWt_Unit                  matlab.ui.control.Label
        ParametersforPorousStructurePanel  matlab.ui.container.Panel
        SurfaceAreaEditFieldLabel     matlab.ui.control.Label
        SurfaceArea_BJH               matlab.ui.control.EditField
        PoreVolumeEditFieldLabel      matlab.ui.control.Label
        PoreVolume_BJH                matlab.ui.control.EditField
        MeanRadiusEditFieldLabel      matlab.ui.control.Label
        MeanRadius_BJH                matlab.ui.control.EditField
        PorosityEditFieldLabel        matlab.ui.control.Label
        Porosity_BJH                  matlab.ui.control.EditField
        PeakRadiusEditFieldLabel      matlab.ui.control.Label
        PeakRadius_BJH                matlab.ui.control.EditField
        SurfaceArea_DBJH              matlab.ui.control.EditField
        PoreVolume_DBJH               matlab.ui.control.EditField
        MeanRadius_DBJH               matlab.ui.control.EditField
        Porosity_DBJH                 matlab.ui.control.EditField
        PeakRadius_DBJH               matlab.ui.control.EditField
        Lable_BJH                     matlab.ui.control.Label
        Lable_DBJH                    matlab.ui.control.Label
        MeanRadius_Unit               matlab.ui.control.Label
        Porosity_Unit                 matlab.ui.control.Label
        PoreVolume_Unit               matlab.ui.control.Label
        SurfaceArea_Unit              matlab.ui.control.Label
        PeakRadius_Unit               matlab.ui.control.Label
        CopyRightMingYANGinSUSTechEMailyoung_94126comLabel  matlab.ui.control.Label
    end

    
    properties (Access = private)          % Global Variables
        isofile   % Name of Isoterm File
        isopath   % Path of Isotherm File
        LEN       % Length of Total Adsorption Data
        AdsN      % Number of Adsorption Branch Points
        DesN      % Number of Desorption Branch Points
        
        partAdsorV    % Adsorption Branch Volume Original
        partDesorV    % Desorption Branch Volume Original
        partAdsorp    % Adsorption Branch Relative Pressure Original
        partDesorp    % Desorption Branch Relative Pressure Original
        
        AdsorV    % Adsorption Branch Volume
        DesorV    % Desorption Branch Volume
        Adsorp    % Adsorption Branch Relative Pressure
        Desorp    % Desorption Branch Relative Pressure
        
        AdsorR    % Adsorption Branch Radii
        DesorR    % Desorption Branch Radii
        Adsort    % Adsorption Branch Multilayer Thickness
        Desort    % Desorption Branch Multilayer Thickness
        AdsorRk    % Adsorption Branch Kelvin Radii
        DesorRk    % Desorption Branch Kelvin Radii
        
        AdsorS    % Adsorption Branch Slope
        DesorS    % Desorption Branch Slope
        
        AdsorLr
        DesorLr
        AdsorPSD  % Adsorption Branch Pore Size Distribution
        DesorPSD  % Desorption Branch Pore Size Distribution
        TotalVAD   % Total Volume Accumulative Distribution
        CapDesPSD    % PSD only consider Capillary
        
        BJH_PSD_V
        BJH_PSD_R
        
        Kerogen_Density
        Molec_Wt
        Liquid_Density        
        Molar_Volume
    end
    
    methods (Access = private)
        
        function desor_t = fun_desor_halsey(~,desor_p)
            desor_t = 0.354*(5./(-log(desor_p)).^(1/3));
        end
        function desor_rk = fun_desor_kelvin(~,desor_p)
            desor_rk = -0.9533./(log(desor_p));
        end
        
        function adsor_t = fun_adsor_halsey(~,adsor_p)
            adsor_t = 0.354*(5./(-log(adsor_p)).^(1/3));
        end
        function adsor_rk = fun_adsor_kelvin(~,adsor_p)
            adsor_rk = -0.9533./(log(adsor_p));
        end
        
        function grad_result = Fun_Grad_Iso(app)
            grad_result=0;
            
            % Thickness & Kelvin Radii of Desorption Branch
            app.Desort = fun_desor_halsey(app,app.Desorp);
            app.DesorRk = fun_desor_kelvin(app,app.Desorp);
            app.DesorR = app.Desort + app.DesorRk;
            
            % Thickness & Kelvin Radii of Adsorption Branch
            app.Adsort = 0.354*(5./(-log(app.Adsorp)).^(1/3));
            app.AdsorRk = -0.9533./(log(app.Adsorp));
            app.AdsorR = app.Adsort + app.AdsorRk;
            
            % Slope
            app.DesorS = gradient(app.DesorV,app.Desorp);
            app.AdsorS = gradient(app.AdsorV,app.Adsorp);
            
            % log plot of Slope: Both Adsorption and Desorption
            semilogy(app.UISlope,app.AdsorR,app.AdsorS,...
                app.DesorR,app.DesorS);
            legend(app.UISlope,{'Adsorption Branch','Desorption Branch'},'Location','southeast')
            xlim(app.UISlope,[0 50]);
        end
        
        function psd_result = Fun_CAL_Desor_PSD(app)  % Function for Calculating PSD From Desorption Branch
            rho_liquid_N2 = app.Liquid_Density; % g/L
            rho_gas_N2 = 0.001*app.Molec_Wt/app.Molar_Volume; % g/L
            volume_trans_const = rho_gas_N2/rho_liquid_N2;   % volume transform ratio
            psd_result = 0;
            Mat_a = zeros(0,0);      % Creat a 2D Empty Matrix
            Mat_b = Mat_a;
            a_v1  = zeros(app.DesN,1);
            b_v1 = zeros(app.DesN,1);
            for di = 2:app.DesN
                ab_0 = zeros(di - 2, 1);
                
                % Deduction of Matrix [a]
                a_idx1 = di:app.DesN;
                a_idx1 = [a_idx1, app.DesN]';
                a_idx2 = (di-1):(app.DesN - 1);
                a_idx2 = [di-1, a_idx2]';
                ai = [ab_0; (app.DesorR(a_idx1)-app.DesorR(a_idx2))];
                Mat_a = [Mat_a, ai];
                
                % Deduction of Matrix [b]
                b_idx1 = di:app.DesN;
                b_idx2 = b_idx1 - 1;
                rb_1 = app.DesorR(b_idx1');
                rb_2 = app.DesorR(b_idx2');
                b1 = [rb_1.*rb_1 + rb_1.*rb_2 - 2*(rb_2.*rb_2);0];
                b2 = [0; 2*(rb_1.*rb_1) - rb_1.*rb_2 - rb_2.*rb_2];
                bi = [ab_0; b1 + b2];
                Mat_b = [Mat_b, bi];
            end
            % Because the adsorption volume contribution in the largest
            % pressure is very small, set it to zero.
            a_v1 = Mat_a(:,1);
            b_v1 = Mat_b(:,1);
            %             a_1 = app.DesorR(1)/2;
            %             b_1 = app.DesorR(1)*app.DesorR(1)/3;
            %             a_v1(1) = a_v1(1) + a_1;
            %             b_v1(1) = b_v1(1) + b_1;
            Mat_a = [a_v1,Mat_a];
            Mat_b = [b_v1,Mat_b];
            %             a_N = zeros(app.DesN,1);
            %             b_N = zeros(app.DesN,1);
            %             Mat_a = [Mat_a,a_N];
            %             Mat_b = [Mat_b,b_N];
            
            % generate dt/dp and Matrix [c] & [ct] & [C]
            vect_dtp = gradient(app.Desort, app.Desorp);
            vect_dtp_t = vect_dtp.*app.Desort;
            vect_drp = gradient(app.DesorR, app.Desorp);
            DesorRk2 = app.DesorRk;  %[0.4;app.DesorRk(1:app.DesN-1)];
            
            vect_m = (app.DesorRk.*DesorRk2).*vect_drp;
            Mat_dtp = pi * diag(vect_dtp);
            Mat_dtp_t = pi * diag(vect_dtp_t);
            Mat_m = pi * diag(vect_m);
            Mat_TC = (1/3) * (Mat_dtp*(Mat_b')) - 2*Mat_dtp_t*(Mat_a');
            %Mat_TC(:,app.DesN-2) = 1E-3*Mat_TC(:,app.DesN-2);
            %Mat_TC(:,app.DesN-1) = 1E-3*Mat_TC(:,app.DesN-1);
            %Mat_TC(:,app.DesN) = 1E-3*Mat_TC(:,app.DesN);
            Mat_TA = Mat_m + Mat_TC;
            
            %app.DesorLr = Mat_TA\app.DesorS;
            opts.UT = true;
            app.DesorLr = linsolve(Mat_TA, volume_trans_const*app.DesorS, opts);    % Fast Solution Method
            app.DesorPSD = pi * ((app.DesorR.*app.DesorR).*app.DesorLr);
            LR = Mat_m\(volume_trans_const*app.DesorS);
            app.CapDesPSD = pi * ((app.DesorR.*app.DesorR).*LR);
            delta_Rp = app.DesorR(2:app.DesN) - app.DesorR(1:app.DesN-1);
            delta_Rp= [app.DesorR(1);delta_Rp];
            app.TotalVAD = tril(ones(app.DesN))*(delta_Rp.*app.DesorPSD);
            
            DBJH_SAD = 2*pi * (app.DesorR.*app.DesorLr);
            DBJH_TAD = tril(ones(app.DesN))*(delta_Rp.*DBJH_SAD);
            DBJH_SSA = num2str(max(100*DBJH_TAD));
            DBJH_PV = num2str(max(app.TotalVAD)-app.TotalVAD(1));
            DBJH_Porosity = num2str(100*(max(app.TotalVAD)-app.TotalVAD(1))*app.Kerogen_Density);
            DBJH_mean_R = num2str(sqrt(abs(sum(delta_Rp.*app.DesorPSD)/sum(delta_Rp.*app.DesorLr)/pi)));
            [~,idx] = max(app.DesorPSD);
            DBJH_peak_R = num2str(app.DesorR(idx));
            
            app.SurfaceArea_DBJH.Value = DBJH_SSA;
            app.PoreVolume_DBJH.Value = DBJH_PV;
            app.Porosity_DBJH.Value = DBJH_Porosity;
            app.MeanRadius_DBJH.Value = DBJH_mean_R;
            app.PeakRadius_DBJH.Value = DBJH_peak_R;
            % log plot of PSD: Desorption
            %semilogy(app.UIPSD,app.DesorR,app.DesorPSD);
            %             yyaxis(app.UIPSD,"left");
            semilogy(app.UIPSD, app.DesorR, app.DesorPSD,...
                app.DesorR, app.CapDesPSD);
            %             yyaxis(app.UIPSD,"right");
            %             plot(app.UIPSD, app.DesorR, app.TotalVAD);
            legend(app.UIPSD,{'D-BJH PSD','Capillary PSD'},'Location','northeast')
            xlim(app.UIPSD,[0 50]);
        end
        
        function [BrDesorp, BrDesorV] = lagrange_intsect(~, Desorp, DesorV, LocalN)     % also Lagrange Intersection with vector
            amplifyP = 1E+5;
            amplifyV = 1E+2;
            Desorp = amplifyP*Desorp;
            DesorV = amplifyV*DesorV;
            OriDesN = length(Desorp);
            DesMinP = Desorp(1);
            DesMaxP = Desorp(OriDesN);
            IncP = (DesMaxP - DesMinP)/(LocalN - 1);
            BrDesorp = (DesMinP:IncP:DesMaxP)';
            BrDesorV = zeros(LocalN,1);
            for Li = 1:LocalN
                localP = BrDesorp(Li);
                LL = 1:OriDesN;
                for ii = LL
                    Lidx = find(LL~=ii);
                    BrDesorV(Li) = BrDesorV(Li) + DesorV(ii)*(prod(localP - Desorp(Lidx))/prod(Desorp(ii) - Desorp(Lidx)));
                end
            end
            BrDesorp = BrDesorp/amplifyP;
            BrDesorV = BrDesorV/amplifyV;
        end
        
        function [BrDesorp, BrDesorV]=fun_intsect_iso(~, x0,y0,localDesN)  % Lagrange Vector Intersection
            OriDesN = length(x0);
            DesMinP = x0(1);
            DesMaxP = x0(OriDesN);
            IncP = (DesMaxP - DesMinP)/(localDesN-1);
            BrDesorp = (DesMinP:IncP:DesMaxP)';
            BrDesorV = zeros(localDesN,1);
            ii=2:(localDesN-1);
            BrDesorV(1) = y0(1);
            BrDesorV(localDesN) = y0(OriDesN);
            for Li=ii
                localP = BrDesorp(Li);
                BrDesorV(Li) = interp1(x0,y0,localP,'pchip');     % pchip is the best method for this problem
            end
        end
        
        function [BJH_PSD_r,BJH_PSD_v] = fun_bjh_psd_cal(app, DeSorpt_p, DeSorpt_v)
            rho_liquid_N2 = app.Liquid_Density; % g/L
            rho_gas_N2 = app.Molec_Wt/app.Molar_Volume; % g/L
            BJH_V_trans = 0.001*rho_gas_N2/rho_liquid_N2;   % volume transform ratio
            
            [DeSorpt_p,idx] = sort(DeSorpt_p,'descend');
            DeSorpt_v = DeSorpt_v(idx);           
            len_des = length(DeSorpt_p);
            bjh_rk = fun_desor_kelvin(app,DeSorpt_p);
            mean_bjh_rk = (bjh_rk(1:len_des-1)+bjh_rk(2:len_des))/2;
            bjh_t = fun_desor_halsey(app,DeSorpt_p);
            bjh_delta_t = bjh_t(1:len_des-1) - bjh_t(2:len_des);
            bjh_rp = bjh_rk + bjh_t;
            mean_bjh_rp = (bjh_rp(1:len_des-1)+bjh_rp(2:len_des))/2;
            bjh_delta_rp = bjh_rp(1:len_des-1) - bjh_rp(2:len_des);
            %             bjh_dp = 2.0*bjh_rp;
            %             mean_bjh_dp = 2.0*mean_bjh_rp;
            bjh_Q = (mean_bjh_rp./(mean_bjh_rk + bjh_delta_t)).^2;
            bjh_delta_v_liquid = BJH_V_trans * (DeSorpt_v(1:len_des-1) - DeSorpt_v(2:len_des));
            bjh_delta_vt = zeros(len_des-1,1);
            bjh_delta_vk = bjh_delta_vt;
            bjh_delta_vp = bjh_delta_vt;
            bjh_delta_ap = bjh_delta_vt;   % Pore Surface Area Distribution
            bjh_sum_vp = bjh_delta_vt;
            bjh_sum_ap = bjh_delta_vt;
            
            bjh_delta_vt(1) = 0;
            for ii = 1:len_des-1
                delta_vk = bjh_delta_v_liquid(ii) - bjh_delta_vt(ii);
                delta_vk(delta_vk < 0) = 0;
                bjh_delta_vk(ii) = delta_vk;
                bjh_delta_vp(ii) = bjh_delta_vk(ii)*bjh_Q(ii);
                bjh_delta_ap(ii) = 2.0*1000.0*bjh_delta_vp(ii)/bjh_rp(ii+1);  % notice cm^3/g  /  nm  == 10^-6/10^-9 ; so ap is in m^2/g
                if(ii==1)
                    bjh_sum_vp(ii) = bjh_delta_vp(ii);
                    bjh_sum_ap(ii) = bjh_delta_ap(ii);  % so sum_ap is also in m^2/g
                else
                    bjh_sum_vp(ii) = bjh_sum_vp(ii-1) + bjh_delta_vp(ii);
                    bjh_sum_ap(ii) = bjh_sum_ap(ii-1) + bjh_delta_ap(ii);
                end
                
                if(ii<len_des-1)
                    bjh_delta_vt(ii+1) = 0.85*0.001 * bjh_delta_t(ii)*bjh_sum_ap(ii);
                end
            end
            BJH_PSD_v = bjh_delta_vp./bjh_delta_rp;
            BJH_PSD_r = mean_bjh_rp;
            BJH_SAD = bjh_delta_ap.*bjh_delta_rp;
            BJH_SSA = num2str(max(bjh_sum_ap));
            BJH_PV = num2str(sum(bjh_delta_vp));
            
            BJH_Porosity = num2str(100*sum(bjh_delta_vp)*app.Kerogen_Density);
            BJH_LR = BJH_PSD_v./(BJH_PSD_r.*BJH_PSD_r)/pi;
            BJH_mean_R = num2str(sqrt(sum(bjh_delta_rp.*BJH_PSD_v)/sum(bjh_delta_rp.*BJH_LR)/pi));
            %             BJH_mean_R = num2str((bjh_delta_vp'*BJH_PSD_r)/sum(bjh_delta_vp));
            [~,idx] = max(BJH_PSD_v);
            BJH_peak_R = num2str(BJH_PSD_r(idx));
            
            app.SurfaceArea_BJH.Value = BJH_SSA;
            app.PoreVolume_BJH.Value = BJH_PV;
            app.Porosity_BJH.Value = BJH_Porosity;
            app.MeanRadius_BJH.Value = BJH_mean_R;
            app.PeakRadius_BJH.Value = BJH_peak_R;
           
            BJH_VAD = triu(ones(len_des-1))*bjh_delta_vp;
            %             yyaxis(app.UIBJHPSD,"left");
            semilogy(app.UIBJHPSD,BJH_PSD_r,BJH_PSD_v,...
                app.DesorR, app.DesorPSD);
            %             yyaxis(app.UIBJHPSD,"right");
            %             plot(app.UIBJHPSD, BJH_PSD_r, BJH_VAD);
            legend(app.UIBJHPSD,{'BJH PSD','D-BJH PSD'},'Location','northeast')
            xlim(app.UIBJHPSD,[0 50]);
        end
    end

    % Callbacks that handle component events
    methods (Access = private)

        % Button pushed function: OpenFileButton
        function OpenFileButtonPushed(app, event)
            [app.isofile,app.isopath]=uigetfile('*.txt;*.mat;*.xlsx;*.xls','Pick an isotherm');
            if isequal(app.isofile,0)||isequal(app.isopath,0)
                disp('User selected Cancel');
            else
                disp(['User selected ', fullfile(app.isopath,app.isofile)]);
            end
            fpath = fullfile(app.isopath,app.isofile);
            app.FileAddress.Value = fpath;
            iso_data = zeros(500,2);
            if contains(fpath,'xls')
                ori_data=xlsread(fpath,'吸附等温线');
                if isempty(ori_data)
                    ori_data = xlsread(fpath,'sheet1');
                end
                [row,~]=size(ori_data);
                iso_data(1:row-1,1)=ori_data(2:row,1);
                iso_data(1:row-1,2)=ori_data(2:row,2);
            else
                iso_data = load(fpath);
            end
            I=find(iso_data(:,1) ~= 0);
            localLEN = length(I);
            [~,localAdsN] = max(iso_data(:,1));
            
            [app.partAdsorp, idx] = unique(iso_data(1:localAdsN,1));
            app.partAdsorV = iso_data(1:localAdsN,2);
            app.partAdsorV = app.partAdsorV(idx);
            
            app.partDesorp = sort(unique(iso_data(localAdsN:localLEN,1)));
            app.partDesorV = sort(unique(iso_data(localAdsN:localLEN,2)));
            
            plot(app.UISO,app.partAdsorp,app.partAdsorV,...
                app.partDesorp,app.partDesorV);
            legend(app.UISO,{'Adsorption Branch','Desorption Branch'},'Location','northwest')
            xlim(app.UISO,[0 1]);
            
            %             Expansion
            app.AdsN = 200;
            app.DesN = 100;
            app.LEN = app.AdsN + app.DesN;
            [app.Adsorp, app.AdsorV] = fun_intsect_iso(app, app.partAdsorp,app.partAdsorV,app.AdsN);
            [app.Desorp, app.DesorV] = fun_intsect_iso(app, app.partDesorp,app.partDesorV,app.DesN);
            
            % Without Expansion  % this way is not good enough!
            %             app.AdsN = localAdsN;
            %             app.LEN = localLEN;
            %             app.DesN = localLEN - localAdsN;
            %             app.Adsorp = app.partAdsorp;
            %             app.AdsorV = app.partAdsorV;
            %             app.Desorp = app.partDesorp;
            %             app.DesorV = app.partDesorp;
        end

        % Button pushed function: CalculateButton
        function CalculateButtonPushed(app, event)
            app.Kerogen_Density = str2double(app.KerogenDensity_val.Value);
            app.Molec_Wt = str2double(app.MolecWt_val.Value);
            app.Liquid_Density = str2double(app.LiquidDensity_val.Value);
            app.Molar_Volume = 22.4; %mol/L
            re = Fun_Grad_Iso(app);
            re2 = Fun_CAL_Desor_PSD(app);
            [app.BJH_PSD_R, app.BJH_PSD_V] = fun_bjh_psd_cal(app, app.Desorp,app.DesorV);
        end

        % Button pushed function: ExportDataButton
        function ExportDataButtonPushed(app, event)
            filename = 'DBJH_Model_for_Desorption_Branch.xlsx';
            default_file = fullfile(app.isopath,filename);
            [saveisofile,saveisopath]=uiputfile('*.xls;*.xlsx;*.mat','Choose a Direct and Set the File Name',default_file);
            save_file = fullfile(saveisopath,saveisofile);
            data_des_save = [app.Desorp,app.Desort,app.DesorRk,app.DesorR,app.DesorS,app.DesorPSD,app.CapDesPSD,app.BJH_PSD_R,app.BJH_PSD_V];
            [m,n] = size(data_des_save);
            data_des_save_cell = mat2cell(data_des_save,ones(m,1),ones(n,1)); % cut the data mat into cell with the size of m*n
            data_des_title = {'Reletive_Pressure','Thickness_Multi_Adsor','Kelvin_Radii','Pore_Radii','Splope_Isotherm','PSD','Cap_PSD','BJH_R','BJH_PSD'};
            data_des = [data_des_title;data_des_save_cell];
            xlswrite(save_file,data_des,'desorption');
        end

        % Value changed function: NormalY
        function NormalYValueChanged(app, event)
            value = app.NormalY.Value;
            cla(app.UISO,'reset');
            cla(app.UIBJHPSD,'reset');
            cla(app.UIPSD,'reset');
            cla(app.UISlope,'reset');
            % Isotherm
            plot(app.UISO,app.partAdsorp,app.partAdsorV,...
                app.partDesorp,app.partDesorV);
            legend(app.UISO,{'Adsorption Branch','Desorption Branch'},'Location','northwest')
            xlim(app.UISO,[0 1]);
            % BJH PSD
            plot(app.UIBJHPSD,app.BJH_PSD_R,app.BJH_PSD_V,...
                app.DesorR, app.DesorPSD);
            legend(app.UIBJHPSD,{'BJH PSD','D-BJH PSD'},'Location','northeast')
            xlim(app.UIBJHPSD,[0 50]);
            % Pore Size Distribution
            plot(app.UIPSD,app.DesorR,app.DesorPSD,...
                app.DesorR,app.CapDesPSD);
            legend(app.UIPSD,{'Desorption Branch','Capillary PSD'},'Location','northeast')
            xlim(app.UIPSD,[0 50]);
            % Slope
            plot(app.UISlope,app.AdsorR,app.AdsorS,...
                app.DesorR,app.DesorS);
            legend(app.UISlope,{'Adsorption Branch','Desorption Branch'},'Location','southeast')
            xlim(app.UISlope,[0 50]);
            app.NormalY.BackgroundColor = 'cyan';
            app.LogY.BackgroundColor = 'w';
        end

        % Value changed function: LogY
        function LogYValueChanged(app, event)
            value = app.LogY.Value;
            % Isotherm
            semilogy(app.UISO,app.partAdsorp,app.partAdsorV,...
                app.partDesorp,app.partDesorV);
            legend(app.UISO,{'Adsorption Branch','Desorption Branch'},'Location','northwest')
            xlim(app.UISO,[0 1]);
            % Re Isotherm
            semilogy(app.UIBJHPSD,app.BJH_PSD_R,app.BJH_PSD_V,...
                app.DesorR, app.DesorPSD);
            legend(app.UIBJHPSD,{'BJH PSD','D-BJH PSD'},'Location','northeast')
            xlim(app.UIBJHPSD,[0 50]);
            % Pore Size Distribution
            semilogy(app.UIPSD,app.DesorR,app.DesorPSD,...
                app.DesorR,app.CapDesPSD);
            legend(app.UIPSD,{'Desorption Branch','Capillary PSD'},'Location','northeast')
            xlim(app.UIPSD,[0 50]);
            % Slope
            semilogy(app.UISlope,app.AdsorR,app.AdsorS,...
                app.DesorR,app.DesorS);
            legend(app.UISlope,{'Adsorption Branch','Desorption Branch'},'Location','southeast')
            xlim(app.UISlope,[0 50]);
            app.NormalY.BackgroundColor = 'w';
            app.LogY.BackgroundColor = 'cyan';
        end
    end

    % Component initialization
    methods (Access = private)

        % Create UIFigure and components
        function createComponents(app)

            % Create DBJH_FIG and hide until all components are created
            app.DBJH_FIG = uifigure('Visible', 'off');
            app.DBJH_FIG.Color = [1 1 1];
            app.DBJH_FIG.Position = [100 100 966 749];
            app.DBJH_FIG.Name = 'MATLAB App';

            % Create OperatePanel
            app.OperatePanel = uipanel(app.DBJH_FIG);
            app.OperatePanel.ForegroundColor = [0.502 0.502 0.502];
            app.OperatePanel.Title = 'Operate Panel';
            app.OperatePanel.BackgroundColor = [0.902 0.902 0.902];
            app.OperatePanel.Position = [12 671 942 67];

            % Create OpenFileButton
            app.OpenFileButton = uibutton(app.OperatePanel, 'push');
            app.OpenFileButton.ButtonPushedFcn = createCallbackFcn(app, @OpenFileButtonPushed, true);
            app.OpenFileButton.Position = [689 11 72 28];
            app.OpenFileButton.Text = 'Open File';

            % Create CalculateButton
            app.CalculateButton = uibutton(app.OperatePanel, 'push');
            app.CalculateButton.ButtonPushedFcn = createCallbackFcn(app, @CalculateButtonPushed, true);
            app.CalculateButton.BackgroundColor = [0.851 0.3255 0.098];
            app.CalculateButton.Position = [770 11 72 28];
            app.CalculateButton.Text = 'Calculate';

            % Create ExportDataButton
            app.ExportDataButton = uibutton(app.OperatePanel, 'push');
            app.ExportDataButton.ButtonPushedFcn = createCallbackFcn(app, @ExportDataButtonPushed, true);
            app.ExportDataButton.BackgroundColor = [0.0745 0.6235 1];
            app.ExportDataButton.FontColor = [0.149 0.149 0.149];
            app.ExportDataButton.Position = [853 11 79 28];
            app.ExportDataButton.Text = 'Export Data';

            % Create FileAddress_TXT
            app.FileAddress_TXT = uilabel(app.OperatePanel);
            app.FileAddress_TXT.HorizontalAlignment = 'right';
            app.FileAddress_TXT.Position = [11 14 72 22];
            app.FileAddress_TXT.Text = 'File Address';

            % Create FileAddress
            app.FileAddress = uieditfield(app.OperatePanel, 'text');
            app.FileAddress.Position = [99 14 583 22];
            app.FileAddress.Value = 'choose a file ... ...';

            % Create DisplayPanel
            app.DisplayPanel = uipanel(app.DBJH_FIG);
            app.DisplayPanel.Title = 'Display Panel';
            app.DisplayPanel.BackgroundColor = [0.9412 0.9412 0.9412];
            app.DisplayPanel.Position = [12 14 722 650];

            % Create UISO
            app.UISO = uiaxes(app.DisplayPanel);
            title(app.UISO, 'Measured Isotherm Data')
            xlabel(app.UISO, 'relative pressure p [-]')
            ylabel(app.UISO, 'adsorption volume V [ml/g]')
            zlabel(app.UISO, 'Z')
            app.UISO.FontName = 'Times New Roman';
            app.UISO.XColor = [0 0 0];
            app.UISO.YColor = [0 0 0];
            app.UISO.ZColor = [0 0 0];
            app.UISO.FontSize = 12;
            app.UISO.Position = [14 331 341 278];

            % Create UIBJHPSD
            app.UIBJHPSD = uiaxes(app.DisplayPanel);
            title(app.UIBJHPSD, 'BJH Pore Size Distribution')
            xlabel(app.UIBJHPSD, 'pore size R [nm]')
            ylabel(app.UIBJHPSD, 'PSD dV/dr  [ml/(g·nm)]')
            zlabel(app.UIBJHPSD, 'Z')
            app.UIBJHPSD.FontName = 'Times New Roman';
            app.UIBJHPSD.XLim = [0 50];
            app.UIBJHPSD.YScale = 'log';
            app.UIBJHPSD.YMinorTick = 'on';
            app.UIBJHPSD.Position = [16 30 339 278];

            % Create UISlope
            app.UISlope = uiaxes(app.DisplayPanel);
            title(app.UISlope, 'Slope of Adsorption Isotherm')
            xlabel(app.UISlope, 'pore size R  [nm]')
            ylabel(app.UISlope, 'Slope  dV/dp  [ml/g]')
            zlabel(app.UISlope, 'Z')
            app.UISlope.FontName = 'Times New Roman';
            app.UISlope.XLim = [0 50];
            app.UISlope.YScale = 'log';
            app.UISlope.YMinorTick = 'on';
            app.UISlope.Position = [369 331 336 278];

            % Create UIPSD
            app.UIPSD = uiaxes(app.DisplayPanel);
            title(app.UIPSD, 'Pore Size Distribution')
            xlabel(app.UIPSD, 'pore size R  [nm]')
            ylabel(app.UIPSD, 'PSD dV/dr  [ml/(g·nm)]')
            zlabel(app.UIPSD, 'Z')
            app.UIPSD.FontName = 'Times New Roman';
            app.UIPSD.XLim = [0 50];
            app.UIPSD.YScale = 'log';
            app.UIPSD.YMinorTick = 'on';
            app.UIPSD.Position = [369 30 336 278];

            % Create SwitchYaxisformButtonGroup
            app.SwitchYaxisformButtonGroup = uibuttongroup(app.DBJH_FIG);
            app.SwitchYaxisformButtonGroup.Title = 'Switch Y-axis form';
            app.SwitchYaxisformButtonGroup.FontName = 'Times New Roman';
            app.SwitchYaxisformButtonGroup.Position = [743 562 211 102];

            % Create NormalY
            app.NormalY = uibutton(app.SwitchYaxisformButtonGroup, 'state');
            app.NormalY.ValueChangedFcn = createCallbackFcn(app, @NormalYValueChanged, true);
            app.NormalY.Text = 'Normal Y Axis';
            app.NormalY.BackgroundColor = [1 1 1];
            app.NormalY.FontName = 'Times New Roman';
            app.NormalY.Position = [8 45 195 26];

            % Create LogY
            app.LogY = uibutton(app.SwitchYaxisformButtonGroup, 'state');
            app.LogY.ValueChangedFcn = createCallbackFcn(app, @LogYValueChanged, true);
            app.LogY.Text = 'Log Y Axis';
            app.LogY.BackgroundColor = [1 1 1];
            app.LogY.FontName = 'Times New Roman';
            app.LogY.Position = [8 11 195 26];

            % Create Material_Property_Panal
            app.Material_Property_Panal = uipanel(app.DBJH_FIG);
            app.Material_Property_Panal.Title = 'Materials Property : N_2 & Kerogen';
            app.Material_Property_Panal.Position = [743 384 211 171];

            % Create MolecWtEditFieldLabel
            app.MolecWtEditFieldLabel = uilabel(app.Material_Property_Panal);
            app.MolecWtEditFieldLabel.FontName = 'Times New Roman';
            app.MolecWtEditFieldLabel.Position = [11 118 54 22];
            app.MolecWtEditFieldLabel.Text = 'Molec.Wt';

            % Create MolecWt_val
            app.MolecWt_val = uieditfield(app.Material_Property_Panal, 'text');
            app.MolecWt_val.FontName = 'Times New Roman';
            app.MolecWt_val.Position = [100 118 60 22];
            app.MolecWt_val.Value = '28.013';

            % Create LiquidDensityEditFieldLabel
            app.LiquidDensityEditFieldLabel = uilabel(app.Material_Property_Panal);
            app.LiquidDensityEditFieldLabel.FontName = 'Times New Roman';
            app.LiquidDensityEditFieldLabel.Position = [11 82 78 22];
            app.LiquidDensityEditFieldLabel.Text = 'Liquid Density';

            % Create LiquidDensity_val
            app.LiquidDensity_val = uieditfield(app.Material_Property_Panal, 'text');
            app.LiquidDensity_val.FontName = 'Times New Roman';
            app.LiquidDensity_val.Position = [100 82 60 22];
            app.LiquidDensity_val.Value = '0.8083';

            % Create CrossSectionEditFieldLabel
            app.CrossSectionEditFieldLabel = uilabel(app.Material_Property_Panal);
            app.CrossSectionEditFieldLabel.FontName = 'Times New Roman';
            app.CrossSectionEditFieldLabel.Position = [11 46 72 22];
            app.CrossSectionEditFieldLabel.Text = 'Cross Section';

            % Create CrossSection_val
            app.CrossSection_val = uieditfield(app.Material_Property_Panal, 'text');
            app.CrossSection_val.FontName = 'Times New Roman';
            app.CrossSection_val.Position = [100 46 60 22];
            app.CrossSection_val.Value = '16.200';

            % Create KerogenDensityEditFieldLabel
            app.KerogenDensityEditFieldLabel = uilabel(app.Material_Property_Panal);
            app.KerogenDensityEditFieldLabel.FontName = 'Times New Roman';
            app.KerogenDensityEditFieldLabel.Position = [11 11 87 22];
            app.KerogenDensityEditFieldLabel.Text = 'Kerogen Density';

            % Create KerogenDensity_val
            app.KerogenDensity_val = uieditfield(app.Material_Property_Panal, 'text');
            app.KerogenDensity_val.FontName = 'Times New Roman';
            app.KerogenDensity_val.Position = [100 11 60 22];
            app.KerogenDensity_val.Value = '1.65';

            % Create KerogenDensity_Unit
            app.KerogenDensity_Unit = uilabel(app.Material_Property_Panal);
            app.KerogenDensity_Unit.BackgroundColor = [0.9412 0.9412 0.9412];
            app.KerogenDensity_Unit.FontName = 'Times New Roman';
            app.KerogenDensity_Unit.Position = [162 11 41 22];
            app.KerogenDensity_Unit.Text = 'g/cm^3';

            % Create CrossSection_Unit
            app.CrossSection_Unit = uilabel(app.Material_Property_Panal);
            app.CrossSection_Unit.BackgroundColor = [0.9412 0.9412 0.9412];
            app.CrossSection_Unit.FontName = 'Times New Roman';
            app.CrossSection_Unit.Position = [162 46 41 22];
            app.CrossSection_Unit.Text = 'Å^2';

            % Create LiquidDensity_Unit
            app.LiquidDensity_Unit = uilabel(app.Material_Property_Panal);
            app.LiquidDensity_Unit.BackgroundColor = [0.9412 0.9412 0.9412];
            app.LiquidDensity_Unit.FontName = 'Times New Roman';
            app.LiquidDensity_Unit.Position = [161 82 42 22];
            app.LiquidDensity_Unit.Text = 'g/ml';

            % Create MolecWt_Unit
            app.MolecWt_Unit = uilabel(app.Material_Property_Panal);
            app.MolecWt_Unit.BackgroundColor = [0.9412 0.9412 0.9412];
            app.MolecWt_Unit.FontName = 'Times New Roman';
            app.MolecWt_Unit.Position = [162 118 41 22];
            app.MolecWt_Unit.Text = 'g/mol';

            % Create ParametersforPorousStructurePanel
            app.ParametersforPorousStructurePanel = uipanel(app.DBJH_FIG);
            app.ParametersforPorousStructurePanel.Title = 'Parameters for Porous Structure';
            app.ParametersforPorousStructurePanel.Position = [744 143 210 234];

            % Create SurfaceAreaEditFieldLabel
            app.SurfaceAreaEditFieldLabel = uilabel(app.ParametersforPorousStructurePanel);
            app.SurfaceAreaEditFieldLabel.FontName = 'Times New Roman';
            app.SurfaceAreaEditFieldLabel.Position = [7 154 68 22];
            app.SurfaceAreaEditFieldLabel.Text = 'Surface Area';

            % Create SurfaceArea_BJH
            app.SurfaceArea_BJH = uieditfield(app.ParametersforPorousStructurePanel, 'text');
            app.SurfaceArea_BJH.FontName = 'Times New Roman';
            app.SurfaceArea_BJH.Position = [80 154 42 22];

            % Create PoreVolumeEditFieldLabel
            app.PoreVolumeEditFieldLabel = uilabel(app.ParametersforPorousStructurePanel);
            app.PoreVolumeEditFieldLabel.FontName = 'Times New Roman';
            app.PoreVolumeEditFieldLabel.Position = [7 118 68 22];
            app.PoreVolumeEditFieldLabel.Text = 'Pore Volume';

            % Create PoreVolume_BJH
            app.PoreVolume_BJH = uieditfield(app.ParametersforPorousStructurePanel, 'text');
            app.PoreVolume_BJH.FontName = 'Times New Roman';
            app.PoreVolume_BJH.Position = [80 118 42 22];

            % Create MeanRadiusEditFieldLabel
            app.MeanRadiusEditFieldLabel = uilabel(app.ParametersforPorousStructurePanel);
            app.MeanRadiusEditFieldLabel.FontName = 'Times New Roman';
            app.MeanRadiusEditFieldLabel.Position = [7 47 69 22];
            app.MeanRadiusEditFieldLabel.Text = 'Mean Radius';

            % Create MeanRadius_BJH
            app.MeanRadius_BJH = uieditfield(app.ParametersforPorousStructurePanel, 'text');
            app.MeanRadius_BJH.FontName = 'Times New Roman';
            app.MeanRadius_BJH.Position = [80 47 42 22];

            % Create PorosityEditFieldLabel
            app.PorosityEditFieldLabel = uilabel(app.ParametersforPorousStructurePanel);
            app.PorosityEditFieldLabel.FontName = 'Times New Roman';
            app.PorosityEditFieldLabel.Position = [7 83 46 22];
            app.PorosityEditFieldLabel.Text = 'Porosity';

            % Create Porosity_BJH
            app.Porosity_BJH = uieditfield(app.ParametersforPorousStructurePanel, 'text');
            app.Porosity_BJH.FontName = 'Times New Roman';
            app.Porosity_BJH.Position = [80 83 42 22];

            % Create PeakRadiusEditFieldLabel
            app.PeakRadiusEditFieldLabel = uilabel(app.ParametersforPorousStructurePanel);
            app.PeakRadiusEditFieldLabel.FontName = 'Times New Roman';
            app.PeakRadiusEditFieldLabel.Position = [7 10 65 22];
            app.PeakRadiusEditFieldLabel.Text = 'Peak Radius';

            % Create PeakRadius_BJH
            app.PeakRadius_BJH = uieditfield(app.ParametersforPorousStructurePanel, 'text');
            app.PeakRadius_BJH.FontName = 'Times New Roman';
            app.PeakRadius_BJH.Position = [80 10 42 22];

            % Create SurfaceArea_DBJH
            app.SurfaceArea_DBJH = uieditfield(app.ParametersforPorousStructurePanel, 'text');
            app.SurfaceArea_DBJH.FontName = 'Times New Roman';
            app.SurfaceArea_DBJH.Position = [125 154 42 22];

            % Create PoreVolume_DBJH
            app.PoreVolume_DBJH = uieditfield(app.ParametersforPorousStructurePanel, 'text');
            app.PoreVolume_DBJH.FontName = 'Times New Roman';
            app.PoreVolume_DBJH.Position = [125 118 42 22];

            % Create MeanRadius_DBJH
            app.MeanRadius_DBJH = uieditfield(app.ParametersforPorousStructurePanel, 'text');
            app.MeanRadius_DBJH.FontName = 'Times New Roman';
            app.MeanRadius_DBJH.Position = [125 47 42 22];

            % Create Porosity_DBJH
            app.Porosity_DBJH = uieditfield(app.ParametersforPorousStructurePanel, 'text');
            app.Porosity_DBJH.FontName = 'Times New Roman';
            app.Porosity_DBJH.Position = [125 83 42 22];

            % Create PeakRadius_DBJH
            app.PeakRadius_DBJH = uieditfield(app.ParametersforPorousStructurePanel, 'text');
            app.PeakRadius_DBJH.FontName = 'Times New Roman';
            app.PeakRadius_DBJH.Position = [125 10 42 22];

            % Create Lable_BJH
            app.Lable_BJH = uilabel(app.ParametersforPorousStructurePanel);
            app.Lable_BJH.BackgroundColor = [0.902 0.902 0.902];
            app.Lable_BJH.HorizontalAlignment = 'center';
            app.Lable_BJH.FontName = 'Times New Roman';
            app.Lable_BJH.FontColor = [0.4941 0.1843 0.5569];
            app.Lable_BJH.Position = [80 183 42 22];
            app.Lable_BJH.Text = 'BJH';

            % Create Lable_DBJH
            app.Lable_DBJH = uilabel(app.ParametersforPorousStructurePanel);
            app.Lable_DBJH.BackgroundColor = [0.902 0.902 0.902];
            app.Lable_DBJH.HorizontalAlignment = 'center';
            app.Lable_DBJH.FontName = 'Times New Roman';
            app.Lable_DBJH.FontColor = [0.9294 0.6941 0.1255];
            app.Lable_DBJH.Position = [125 183 42 22];
            app.Lable_DBJH.Text = 'D-BJH';

            % Create MeanRadius_Unit
            app.MeanRadius_Unit = uilabel(app.ParametersforPorousStructurePanel);
            app.MeanRadius_Unit.BackgroundColor = [0.9412 0.9412 0.9412];
            app.MeanRadius_Unit.FontName = 'Times New Roman';
            app.MeanRadius_Unit.Position = [170 47 37 22];
            app.MeanRadius_Unit.Text = 'nm';

            % Create Porosity_Unit
            app.Porosity_Unit = uilabel(app.ParametersforPorousStructurePanel);
            app.Porosity_Unit.BackgroundColor = [0.9412 0.9412 0.9412];
            app.Porosity_Unit.FontName = 'Times New Roman';
            app.Porosity_Unit.Position = [170 82 37 22];
            app.Porosity_Unit.Text = '%';

            % Create PoreVolume_Unit
            app.PoreVolume_Unit = uilabel(app.ParametersforPorousStructurePanel);
            app.PoreVolume_Unit.BackgroundColor = [0.9412 0.9412 0.9412];
            app.PoreVolume_Unit.FontName = 'Times New Roman';
            app.PoreVolume_Unit.Position = [168 118 41 22];
            app.PoreVolume_Unit.Text = 'cm^3/g';

            % Create SurfaceArea_Unit
            app.SurfaceArea_Unit = uilabel(app.ParametersforPorousStructurePanel);
            app.SurfaceArea_Unit.BackgroundColor = [0.9412 0.9412 0.9412];
            app.SurfaceArea_Unit.FontName = 'Times New Roman';
            app.SurfaceArea_Unit.Position = [169 154 37 22];
            app.SurfaceArea_Unit.Text = 'm^2/g';

            % Create PeakRadius_Unit
            app.PeakRadius_Unit = uilabel(app.ParametersforPorousStructurePanel);
            app.PeakRadius_Unit.BackgroundColor = [0.9412 0.9412 0.9412];
            app.PeakRadius_Unit.FontName = 'Times New Roman';
            app.PeakRadius_Unit.Position = [170 10 37 22];
            app.PeakRadius_Unit.Text = 'nm';

            % Create CopyRightMingYANGinSUSTechEMailyoung_94126comLabel
            app.CopyRightMingYANGinSUSTechEMailyoung_94126comLabel = uilabel(app.DBJH_FIG);
            app.CopyRightMingYANGinSUSTechEMailyoung_94126comLabel.BackgroundColor = [1 0 0];
            app.CopyRightMingYANGinSUSTechEMailyoung_94126comLabel.FontName = 'Times New Roman';
            app.CopyRightMingYANGinSUSTechEMailyoung_94126comLabel.FontSize = 14;
            app.CopyRightMingYANGinSUSTechEMailyoung_94126comLabel.FontColor = [1 1 1];
            app.CopyRightMingYANGinSUSTechEMailyoung_94126comLabel.Enable = 'off';
            app.CopyRightMingYANGinSUSTechEMailyoung_94126comLabel.Position = [744 14 210 123];
            app.CopyRightMingYANGinSUSTechEMailyoung_94126comLabel.Text = {' CopyRight:'; ''; '    MingYANG in SUSTech'; '    E-Mail: young_94@126.com      '};

            % Show the figure after all components are created
            app.DBJH_FIG.Visible = 'on';
        end
    end

    % App creation and deletion
    methods (Access = public)

        % Construct app
        function app = DBJH_Numerical_exported

            % Create UIFigure and components
            createComponents(app)

            % Register the app with App Designer
            registerApp(app, app.DBJH_FIG)

            if nargout == 0
                clear app
            end
        end

        % Code that executes before app deletion
        function delete(app)

            % Delete UIFigure when app is deleted
            delete(app.DBJH_FIG)
        end
    end
end