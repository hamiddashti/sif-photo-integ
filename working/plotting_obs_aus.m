clear 
clc
close all
load obs
names = ["L1", "L2", "L3","L4","L5","L6","L7","L8","L9"];

aus_licor_co2 = obs.aus_data.licor_co2;
aus_spec_co2 = obs.aus_data.spec_co2;
aus_licor_lrc = obs.aus_data.licor_lrc;
aus_spec_lrc = obs.aus_data.spec_lrc;
%% Plotting co2 experiment

% Plot the gas-exchange variables
figure(1)
set(gcf,'units','inch','Position',[25 25 9 12],'color','w')
mycolor = [.7 .7 .7];
% cmap = brewermap(length(names)-1,"Dark2")
% Legend=cell(length(names)-1,1)
Fsize=11;
for i =1:length(names)-1
    % Reading data
    Ca = aus_licor_co2.Ca(:,i);
    Fs = aus_licor_co2.Fs(:,i);
    Fo = aus_licor_co2.Fo(:,i);
    Fop = aus_licor_co2.Fo_(:,i);
    Fm = aus_licor_co2.Fm(:,i);
    Fmp = aus_licor_co2.Fm_(:,i);
    An = aus_licor_co2.A(:,i);
    gsw = aus_licor_co2.gsw(:,i);
    
    % Calculate PAM indices
    Fs_Fo = Fs./Fo; 
    Fm_Fo = Fm./Fo;
    Fmp_F0 = Fmp./Fo;
    Fv_Fm = 1-(Fo./Fm);
    PhiPS2 = 1-(Fs./Fmp);
    qP = (Fmp-Fs)./(Fmp-Fop);
    qL = qP.*(Fop./Fs);
    qN = (Fm-Fmp)./(Fm-Fop);
    NPQ = (Fm./Fmp)-1;
    Y_NO = Fs./Fm;
    Y_NPQ = (Fs./Fmp)-(Fs./Fm);
    ETR = aus_licor_co2.ETR(:,i);
    qL_1 = 1-qL;
   
    s(1) = subplot(5,3,1)
    hold on
    plot(Ca,An ,".-",color=mycolor)
    box on
    ylabel("A_{n} [µmol m^{-2} s^{-1}]",FontSize=Fsize)
    
    s(2) = subplot(5,3,2)
    hold on
    plot(Ca,gsw,".-",color=mycolor)
    box on
    ylabel("gsw [mol m^{-2} s^{-1}]",FontSize=Fsize)

    s(3) = subplot(5,3,3)
    hold on
    plot(Ca,Fs_Fo,".-",color=mycolor)
    box on
    ylabel("Fs",FontSize=Fsize)

    s(4) = subplot(5,3,4)
    hold on
    plot(Ca,Fm_Fo,".-",color=mycolor)
    box on
    ylabel("Fm",FontSize=Fsize)
    ylim([4.5, 6.5])

    s(5) = subplot(5,3,5)
    hold on
    plot(Ca,Fmp_F0,".-",color=mycolor)
    box on 
    ylabel("Fm'",FontSize=Fsize)

    s(6) = subplot(5,3,6)
    hold on
    plot(Ca,Fv_Fm,".-",color=mycolor)
    box on
    ylabel("Fv/Fm",FontSize=Fsize)

    s(7) = subplot(5,3,7)
    hold on
    plot(Ca,PhiPS2,".-",color=mycolor)
    box on
    ylabel("PhiPS2",FontSize=Fsize)

    s(8) = subplot(5,3,8)
    hold on
    plot(Ca,qP,".-",color=mycolor)
    box on 
    ylabel("qP",FontSize=Fsize)
    
    s(9) = subplot(5,3,9)
    hold on
    plot(Ca,qL,".-",color=mycolor)
    box on
    ylabel("qL",FontSize=Fsize)

    s(10) = subplot(5,3,10)
    hold on
    plot(Ca,qN,".-",color=mycolor)
    box on
    ylabel("qN",FontSize=Fsize)

    s(11) = subplot(5,3,11)
    hold on
    plot(Ca,NPQ,".-",color=mycolor)
    box on 
    ylabel("NPQ",FontSize=Fsize)

    s(12) = subplot(5,3,12)
    hold on
    plot(Ca,Y_NO,".-",color=mycolor)
    box on
    ylabel("Y(NO)",FontSize=Fsize)

    s(13) = subplot(5,3,13)
    hold on
    plot(Ca,Y_NPQ,".-",color=mycolor)
    box on
    ylabel("Y(NPQ)",FontSize=Fsize)
    xlabel("Ambient CO_{2} [µmol mol^{-1}]",FontSize=Fsize)

    s(14) = subplot(5,3,14)
    hold on
    plot(Ca,ETR,".-",color=mycolor)
    box on
    ylabel("ETR [µmol e^{-1} s^{-1}]",FontSize=Fsize)
    xlabel("Ambient CO_{2} [µmol mol^{-1}]",FontSize=Fsize)
    
    s(15) = subplot(5,3,15)
    hold on
    plot(Ca,qL_1,".-",color=mycolor)
    box on
    ylabel("1-qL",FontSize=Fsize)
    xlabel("Ambient CO_{2} [µmol mol^{-1}]",FontSize=Fsize)
      
end

sgtitle("CO_{2} experiment")
% newPosition = [0.5 0.05 0.01 0.01];
% newUnits = 'normalized';
% hL = legend(Legend,"orientation","horizontal","location","southoutside")
% set(hL,'Position', newPosition,'Units', newUnits,"FontSize",12);
saveas(gcf,"./Figures/Aus_GasExc_co2.png")
% ----------------------------------------------------------------------

% Plotting refl/trans/SIF for the co2 experiment
% Smoothing refl/trans and calculate PRI 
mean_co2_ca =mean(aus_licor_co2.Ca,2);

smoothed_refl_co2 = zeros(length(aus_licor_co2.Ca),length(aus_spec_co2.wl_refl_green),...
    length(names)-1);
smoothed_trans_co2 = zeros(length(aus_licor_co2.Ca),length(aus_spec_co2.wl_trans_green),...
    length(names)-1);

smoothed_fmax_back_co2 = zeros(length(aus_licor_co2.Ca),length(aus_spec_co2.wl_back_f),...
    length(names)-1);
smoothed_fmax_forw_co2 = zeros(length(aus_licor_co2.Ca),length(aus_spec_co2.wl_forw_f),...
    length(names)-1);

smoothed_fmin_back_co2 = zeros(length(aus_licor_co2.Ca),length(aus_spec_co2.wl_back_f),...
    length(names)-1);
smoothed_fmin_forw_co2 = zeros(length(aus_licor_co2.Ca),length(aus_spec_co2.wl_forw_f),...
    length(names)-1);

smoothed_fs_back_co2 = zeros(length(aus_licor_co2.Ca),length(aus_spec_co2.wl_back_f),...
    length(names)-1);
smoothed_fs_forw_co2 = zeros(length(aus_licor_co2.Ca),length(aus_spec_co2.wl_forw_f),...
    length(names)-1);

% ------------------- Smooth refl/trans ------------------------------
for i =1:17
    for j=1:8
        smoothed_refl_co2(i,:,j) = sgolayfilt(aus_spec_co2.refl(i,:,j),2,7);
        smoothed_trans_co2(i,:,j) = sgolayfilt(aus_spec_co2.trans(i,:,j),2,7);
    end
end
% Take the mean over 8 observations
smoothed_refl_mean = mean(smoothed_refl_co2,3);
smoothed_trans_mean = mean(smoothed_trans_co2,3);
%---------------------------------------------------------------------

% ------------------- Calculate PRI ----------------------------------
wl_refl_531 = find(aus_spec_co2.wl_refl_green>=531 & aus_spec_co2.wl_refl_green<532);
wl_refl_570 = find(aus_spec_co2.wl_refl_green>=570 & aus_spec_co2.wl_refl_green<571);
wl_trans_531 = find(aus_spec_co2.wl_trans_green>=531 & aus_spec_co2.wl_trans_green<532);
wl_trans_570 = find(aus_spec_co2.wl_trans_green>=570 & aus_spec_co2.wl_trans_green<571);

% PRI from refl
a_refl = squeeze(smoothed_refl_co2(:,wl_refl_531,:));
b_refl = squeeze(smoothed_refl_co2(:,wl_refl_570,:));
pri_refl_co2 = (a_refl-b_refl)./(a_refl+b_refl);

% PRI from trans
a_trans = squeeze(smoothed_trans_co2(:,wl_trans_531,:));
b_trans = squeeze(smoothed_trans_co2(:,wl_trans_570,:));
pri_trans_co2 = (a_trans-b_trans)./(a_trans+b_trans);

% PRI from albedo (refl+trans)
a_albedo = a_refl+a_trans; 
b_albedo = b_refl+b_trans;
pri_albedo_co2 = (a_albedo-b_albedo)./(a_albedo+b_albedo);

% Take the mean for all samples
pri_refl_co2_mean = mean(pri_refl_co2,2);
pri_trans_co2_mean = mean(pri_trans_co2,2);
pri_albedo_co2_mean = mean(pri_albedo_co2,2);
%---------------------------------------------------------------------

% ------------------- Smooth SIF ------------------------------------
for i =1:length(aus_licor_co2.Ca)
    for j=1:length(names)-1
        smoothed_fmax_back_co2(i,:,j) = sgolayfilt(aus_spec_co2.fmax_back(i,:,j),2,7);
        smoothed_fmax_forw_co2(i,:,j) = sgolayfilt(aus_spec_co2.fmax_forw(i,:,j),2,7);

        smoothed_fmin_back_co2(i,:,j) = sgolayfilt(aus_spec_co2.fmin_back(i,:,j),2,7);
        smoothed_fmin_forw_co2(i,:,j) = sgolayfilt(aus_spec_co2.fmin_forw(i,:,j),2,7);

        smoothed_fs_back_co2(i,:,j) = sgolayfilt(aus_spec_co2.fs_back(i,:,j),2,7);
        smoothed_fs_forw_co2(i,:,j) = sgolayfilt(aus_spec_co2.fs_forw(i,:,j),2,7);


    end
end
% Take the mean over 8 observations
smoothed_fmax_back_co2_mean = mean(smoothed_fmax_back_co2,3);
smoothed_fmax_forw_co2_mean = mean(smoothed_fmax_forw_co2,3);

smoothed_fmin_back_co2_mean = mean(smoothed_fmin_back_co2,3);
smoothed_fmin_forw_co2_mean = mean(smoothed_fmin_forw_co2,3);

smoothed_fs_back_co2_mean = mean(smoothed_fs_back_co2,3);
smoothed_fs_forw_co2_mean = mean(smoothed_fs_forw_co2,3);
%--------------------------------------------------------------------

% ----------------------Plot refl/trans and SIF
figure(3)
set(gcf,'units','inch','Position',[1 2 8 9],'color','w')
% Legend=cell(length(aus_licor_co2.Ca),1);
cmap = jet(17);
for i=1:17 % Loop over different Ca

    subplot(3,3,1)
    hold on
    plot(aus_spec_co2.wl_refl_green,smoothed_refl_mean(i,:),color=cmap(i,:))
    box on 
    xlabel("wl [nm]",FontSize=12)
    ylabel("Reflectance",FontSize=12)

    subplot(3,3,2)
    hold on
    plot(aus_spec_co2.wl_trans_green,smoothed_trans_mean(i,:),color=cmap(i,:))
    box on
    xlabel("wl [nm]",FontSize=12)
    ylabel("Transmittance",FontSize=12)

    subplot(3,3,3)
    hold on
    plot(aus_spec_co2.wl_trans_green,smoothed_refl_mean(i,:)+smoothed_trans_mean(i,:),color=cmap(i,:))
    box on
    xlabel("wl [nm]",FontSize=12)
    ylabel("Leaf albedo",FontSize=12)

    subplot(3,3,4)
    hold on
    plot(aus_spec_co2.wl_back_f,smoothed_fmax_back_co2_mean(i,:),color=cmap(i,:))
    box on 
    xlabel("wl [nm]",FontSize=12)
    ylabel("Fmax_{b} \muW cm^{-2} nm^{-1}",FontSize=12)

    subplot(3,3,5)
    hold on
    plot(aus_spec_co2.wl_forw_f,smoothed_fmax_forw_co2_mean(i,:),color=cmap(i,:))
    box on 
    xlabel("wl [nm]",FontSize=12)
    ylabel("Fmax_{f} \muW cm^{-2} nm^{-1}",FontSize=12)

    subplot(3,3,6)
    hold on
    plot(aus_spec_co2.wl_back_f,smoothed_fmin_back_co2_mean(i,:),color=cmap(i,:))
    box on 
    xlabel("wl [nm]",FontSize=12)
    ylabel("Fmin_{b} \muW cm^{-2} nm^{-1}",FontSize=12)

    subplot(3,3,7)
    hold on
    plot(aus_spec_co2.wl_forw_f,smoothed_fmin_forw_co2_mean(i,:),color=cmap(i,:))
    box on 
    xlabel("wl [nm]",FontSize=12)
    ylabel("Fmin_{f} \muW cm^{-2} nm^{-1}",FontSize=12)

    subplot(3,3,8)
    hold on
    plot(aus_spec_co2.wl_back_f,smoothed_fs_back_co2_mean(i,:),color=cmap(i,:))
    box on 
    xlabel("wl [nm]",FontSize=12)
    ylabel("Fs_b \muW cm^{-2} nm^{-1}",FontSize=12)

    subplot(3,3,9)
    hold on
    plot(aus_spec_co2.wl_forw_f,smoothed_fs_forw_co2_mean(i,:),color=cmap(i,:))
    box on 
    xlabel("wl [nm]",FontSize=12)
    ylabel("Fs_f \muW cm^{-2} nm^{-1}",FontSize=12)

end
c= colorbar('southoutside');
colormap(c,'jet')
caxis([min(aus_licor_co2.Ca(:,1)),max(aus_licor_co2.Ca(:,1))])
c.Label.String = "Ambient CO_{2} [µmol mol^{-1}]";
newPosition = [0.35 0.05 0.4 0.015];
newUnits = 'normalized';
c.FontSize = 10;
c.Units = newUnits;
c.Position = newPosition;
sgtitle("CO2 Experiment");
saveas(gcf,"./Figures/Aus_spec_co2.png")
% -----------------------------------------------------------------

% -----------Plot Fs in 757 (Magney 2017) for different Ca----------
figure(4)
set(gcf,'units','inch','Position',[1 2 5 5],'color','w')
I = find(aus_spec_co2.wl_back_f>757 & aus_spec_co2.wl_back_f<757.5);

h(1) = plot(mean_co2_ca,smoothed_fs_back_co2_mean(:,I),'r',LineWidth=1)
hold on
h(2) = plot(mean_co2_ca,smoothed_fs_forw_co2_mean(:,I),'k',LineWidth=1)

xlabel('Ambient CO_{2} [µmol mol^{-1}]')
ylabel('Fs [\muW cm^{-2} nm^{-1}]')
hl = legend(h,'757nm backward','757nm forward')
hl.FontSize = 10;
saveas(gcf,"./Figures/Aus_fs_757_co2.png")
%---------------------------------------------------------------------

% --------------------- Plot PRI -------------------------------------
figure(5)
set(gcf,'units','inch','Position',[1 2 10 3],'color','w')
Legend=cell(length(aus_licor_co2.Ca),1);
cmap = brewermap(length(names)-1,"Dark2");
for i =1:length(names)-1
    subplot(1,3,1)
    hold on
    plot(mean_co2_ca,pri_refl_co2(:,i),".-",color=cmap(i,:))
    plot(mean_co2_ca,pri_refl_co2_mean,"-",color="black",LineWidth=2)
    box on 
    xlabel('Ambient CO_{2} [µmol mol^{-1}]')
    ylabel('PRI')


    subplot(1,3,2)
    hold on
    plot(mean_co2_ca,pri_trans_co2(:,i),".-",color=cmap(i,:))
    plot(mean_co2_ca,pri_trans_co2_mean,"-",color="black",LineWidth=2)
    box on 
    xlabel('Ambient CO_{2} [µmol mol^{-1}]')
    ylabel('PRI')

    subplot(1,3,3)
    hold on
    plot(mean_co2_ca,pri_albedo_co2(:,i),".-",color=cmap(i,:))
    plot(mean_co2_ca,pri_albedo_co2_mean,"-",color="black",LineWidth=2)
    box on 
    xlabel('Ambient CO_{2} [µmol mol^{-1}]')
    ylabel('PRI')
end
sgtitle("CO_{2} experiment")
saveas(gcf,"./Figures/Aus_PRI_co2.png")
% -------------------------------------------------------------------

%% Plotting LRC experiment
figure(6)
set(gcf,'units','inch','Position',[25 25 9 12],'color','w')
cmap = brewermap(length(names),"Dark2");
Legend=cell(length(names),1);
Fsize=11;
for i =1:length(names)
    par = aus_licor_lrc.Qin(:,i);
    Fs = aus_licor_lrc.Fs(:,i);
    Fo = aus_licor_lrc.Fo(:,i);
    Fop = aus_licor_lrc.Fo_(:,i);
    Fm = aus_licor_lrc.Fm(:,i);
    Fmp = aus_licor_lrc.Fm_(:,i);
        
    An = aus_licor_lrc.A(:,i);
    gsw = aus_licor_lrc.gsw(:,i);
    Fs_Fo = Fs./Fo; 
    Fm_Fo = Fm./Fo;
    Fmp_F0 = Fmp./Fo;
    Fv_Fm = 1-(Fo./Fm);
    PhiPS2 = 1-(Fs./Fmp);
    qP = (Fmp-Fs)./(Fmp-Fop);
    qL = qP.*(Fop./Fs);
    qN = (Fm-Fmp)./(Fm-Fop);
    NPQ = (Fm./Fmp)-1;
    Y_NO = Fs./Fm;
    Y_NPQ = (Fs./Fmp)-(Fs./Fm);
    ETR = aus_licor_lrc.ETR(:,i);
    qL_1 = 1-qL;
   
    s(1) = subplot(5,3,1)
    hold on
    plot(par,An ,".-",color=cmap(i,:))
    box on
%     xlabel("Ambient CO_{2} [µmol mol^{-1}]",FontSize=12)
    ylabel("A_{n} [µmol m^{-2} s^{-1}]",FontSize=Fsize)
    
    s(2) = subplot(5,3,2)
    hold on
    plot(par,gsw,".-",color=cmap(i,:))
    box on
%     xlabel("Ambient CO_{2} [µmol mol^{-1}]",FontSize=12)
    ylabel("gsw [mol m^{-2} s^{-1}]",FontSize=Fsize)

    s(3) = subplot(5,3,3)
    hold on
    plot(par,Fs_Fo,".-",color=cmap(i,:))
    box on
    ylabel("Fs",FontSize=Fsize)
    ylim([0.85,2])

    s(4) = subplot(5,3,4)
    hold on
    plot(par,Fm_Fo,".-",color=cmap(i,:))
    box on
    ylabel("Fm",FontSize=Fsize)
    ylim([4.5, 6.5])

    s(5) = subplot(5,3,5)
    hold on
    plot(par,Fmp_F0,".-",color=cmap(i,:))
    box on 
    ylabel("Fm",FontSize=Fsize)

    s(6) = subplot(5,3,6)
    hold on
    plot(par,Fv_Fm,".-",color=cmap(i,:))
    box on
    ylabel("Fv/Fm",FontSize=Fsize)

    s(7) = subplot(5,3,7)
    hold on
    plot(par,PhiPS2,".-",color=cmap(i,:))
    box on
    ylabel("PhiPS2",FontSize=Fsize)

    s(8) = subplot(5,3,8)
    hold on
    plot(par,qP,".-",color=cmap(i,:))
    box on 
    ylabel("qP",FontSize=Fsize)
    
    s(9) = subplot(5,3,9)
    hold on
    plot(par,qL,".-",color=cmap(i,:))
    box on
    ylabel("qL",FontSize=Fsize)

    s(10) = subplot(5,3,10)
    hold on
    plot(par,qN,".-",color=cmap(i,:))
    box on
    ylabel("qN",FontSize=Fsize)

    s(11) = subplot(5,3,11)
    hold on
    plot(par,NPQ,".-",color=cmap(i,:))
    box on 
    ylabel("NPQ",FontSize=Fsize)

    s(12) = subplot(5,3,12)
    hold on
    plot(par,Y_NO,".-",color=cmap(i,:))
    box on
    ylabel("Y(NO)",FontSize=Fsize)

    s(13) = subplot(5,3,13)
    hold on
    plot(par,Y_NPQ,".-",color=cmap(i,:))
    box on
    ylabel("Y(NPQ)",FontSize=Fsize)
    xlabel("PAR [µmol m^{-2} s^{-1}]",FontSize=Fsize)

    s(14) = subplot(5,3,14)
    hold on
    plot(par,ETR,".-",color=cmap(i,:))
    box on
    ylabel("ETR [µmol e^{-1} s^{-1}]",FontSize=Fsize)
    xlabel("PAR [µmol m^{-2} s^{-1}]",FontSize=Fsize)
    
    s(15) = subplot(5,3,15)
    hold on
    plot(par,qL_1,".-",color=cmap(i,:))
    box on
    ylabel("1-qL",FontSize=Fsize)
    xlabel("PAR [µmol m^{-2} s^{-1}]",FontSize=Fsize)
    Legend{i} = names(i);


%     
end
sgtitle("LRC experiment")
newPosition = [0.5 0.05 0.01 0.01];
newUnits = 'normalized';
hL = legend(Legend,"orientation","horizontal","location","southoutside")
set(hL,'Position', newPosition,'Units', newUnits,"FontSize",12);
saveas(gcf,"./Figures/Aus_GasExc_lrc.png")

% Plotting refl/trans/SIF for the co2 experiment
% Smoothing refl/trans and calculate PRI 
mean_lrc_qin =mean(aus_licor_lrc.Qin,2);

smoothed_refl_lrc = zeros(length(aus_licor_lrc.Qin),length(aus_spec_lrc.wl_refl_green),...
    length(names)-1);
smoothed_trans_co2 = zeros(length(aus_licor_lrc.Qin),length(aus_spec_lrc.wl_trans_green),...
    length(names)-1);

smoothed_fmax_back_lrc = zeros(length(aus_licor_lrc.Qin),length(aus_spec_lrc.wl_back_f),...
    length(names)-1);
smoothed_fmax_forw_lrc = zeros(length(aus_licor_lrc.Qin),length(aus_spec_lrc.wl_forw_f),...
    length(names)-1);

smoothed_fmin_back_lrc = zeros(length(aus_licor_lrc.Qin),length(aus_spec_lrc.wl_back_f),...
    length(names)-1);
smoothed_fmin_forw_lrc = zeros(length(aus_licor_lrc.Qin),length(aus_spec_lrc.wl_forw_f),...
    length(names)-1);

smoothed_fs_back_lrc = zeros(length(aus_licor_lrc.Qin),length(aus_spec_lrc.wl_back_f),...
    length(names)-1);
smoothed_fs_forw_lrc = zeros(length(aus_licor_lrc.Qin),length(aus_spec_lrc.wl_forw_f),...
    length(names)-1);

% ------------------- Smooth refl/trans ------------------------------
for i =1:length(aus_licor_lrc.Qin)
    for j=1:9
        smoothed_refl_lrc(i,:,j) = sgolayfilt(aus_spec_lrc.refl(i,:,j),2,7);
        smoothed_trans_lrc(i,:,j) = sgolayfilt(aus_spec_lrc.trans(i,:,j),2,7);
    end
end
% Take the mean over 8 observations
smoothed_refl_mean = mean(smoothed_refl_lrc,3);
smoothed_trans_mean = mean(smoothed_trans_lrc,3);
%---------------------------------------------------------------------

% ------------------- Calculate PRI ----------------------------------
wl_refl_531 = find(aus_spec_lrc.wl_refl_green>=531 & aus_spec_lrc.wl_refl_green<532);
wl_refl_570 = find(aus_spec_lrc.wl_refl_green>=570 & aus_spec_lrc.wl_refl_green<571);
wl_trans_531 = find(aus_spec_lrc.wl_trans_green>=531 & aus_spec_lrc.wl_trans_green<532);
wl_trans_570 = find(aus_spec_lrc.wl_trans_green>=570 & aus_spec_lrc.wl_trans_green<571);

% PRI from refl
a_refl = squeeze(smoothed_refl_lrc(:,wl_refl_531,:));
b_refl = squeeze(smoothed_refl_lrc(:,wl_refl_570,:));
pri_refl_lrc = (a_refl-b_refl)./(a_refl+b_refl);

% PRI from trans
a_trans = squeeze(smoothed_trans_lrc(:,wl_trans_531,:));
b_trans = squeeze(smoothed_trans_lrc(:,wl_trans_570,:));
pri_trans_lrc = (a_trans-b_trans)./(a_trans+b_trans);

% PRI from albedo (refl+trans)
a_albedo = a_refl+a_trans; 
b_albedo = b_refl+b_trans;
pri_albedo_lrc = (a_albedo-b_albedo)./(a_albedo+b_albedo);

% Take the mean for all samples
pri_refl_lrc_mean = mean(pri_refl_lrc,2);
pri_trans_lrc_mean = mean(pri_trans_lrc,2);
pri_albedo_lrc_mean = mean(pri_albedo_lrc,2);
%---------------------------------------------------------------------

% ------------------- Smooth SIF ------------------------------------
for i =1:length(aus_licor_lrc.Qin)
    for j=1:length(names)
        smoothed_fmax_back_lrc(i,:,j) = sgolayfilt(aus_spec_lrc.fmax_back(i,:,j),2,7);
        smoothed_fmax_forw_lrc(i,:,j) = sgolayfilt(aus_spec_lrc.fmax_forw(i,:,j),2,7);

        smoothed_fmin_back_lrc(i,:,j) = sgolayfilt(aus_spec_lrc.fmin_back(i,:,j),2,7);
        smoothed_fmin_forw_lrc(i,:,j) = sgolayfilt(aus_spec_lrc.fmin_forw(i,:,j),2,7);

        smoothed_fs_back_lrc(i,:,j) = sgolayfilt(aus_spec_lrc.fs_back(i,:,j),2,7);
        smoothed_fs_forw_lrc(i,:,j) = sgolayfilt(aus_spec_lrc.fs_forw(i,:,j),2,7);
    end
end

% Take the mean over 8 observations
smoothed_fmax_back_lrc_mean = mean(smoothed_fmax_back_lrc,3);
smoothed_fmax_forw_lrc_mean = mean(smoothed_fmax_forw_lrc,3);

smoothed_fmin_back_lrc_mean = mean(smoothed_fmin_back_lrc,3);
smoothed_fmin_forw_lrc_mean = mean(smoothed_fmin_forw_lrc,3);

smoothed_fs_back_lrc_mean = mean(smoothed_fs_back_lrc,3);
smoothed_fs_forw_lrc_mean = mean(smoothed_fs_forw_lrc,3);
%--------------------------------------------------------------------

% ----------------------Plot refl/trans and SIF
figure(7)
set(gcf,'units','inch','Position',[1 2 8 9],'color','w')
% Legend=cell(length(aus_licor_co2.Ca),1);
cmap = jet(length(aus_licor_lrc.Qin));
for i=1:length(aus_licor_lrc.Qin) % Loop over different Ca

    subplot(3,3,1)
    hold on
    plot(aus_spec_lrc.wl_refl_green,smoothed_refl_mean(i,:),color=cmap(i,:))
    box on 
    xlabel("wl [nm]",FontSize=12)
    ylabel("Reflectance",FontSize=12)

    subplot(3,3,2)
    hold on
    plot(aus_spec_lrc.wl_trans_green,smoothed_trans_mean(i,:),color=cmap(i,:))
    box on
    xlabel("wl [nm]",FontSize=12)
    ylabel("Transmittance",FontSize=12)

    subplot(3,3,3)
    hold on
    plot(aus_spec_lrc.wl_trans_green,smoothed_refl_mean(i,:)+smoothed_trans_mean(i,:),color=cmap(i,:))
    box on
    xlabel("wl [nm]",FontSize=12)
    ylabel("Leaf albedo",FontSize=12)

    subplot(3,3,4)
    hold on
    plot(aus_spec_lrc.wl_back_f,smoothed_fmax_back_lrc_mean(i,:),color=cmap(i,:))
    box on 
    xlabel("wl [nm]",FontSize=12)
    ylabel("Fmax_{b} \muW cm^{-2} nm^{-1}",FontSize=12)

    subplot(3,3,5)
    hold on
    plot(aus_spec_lrc.wl_forw_f,smoothed_fmax_forw_lrc_mean(i,:),color=cmap(i,:))
    box on 
    xlabel("wl [nm]",FontSize=12)
    ylabel("Fmax_{f} \muW cm^{-2} nm^{-1}",FontSize=12)

    subplot(3,3,6)
    hold on
    plot(aus_spec_co2.wl_back_f,smoothed_fmin_back_co2_mean(i,:),color=cmap(i,:))
    box on 
    xlabel("wl [nm]",FontSize=12)
    ylabel("Fmin_{b} \muW cm^{-2} nm^{-1}",FontSize=12)

    subplot(3,3,7)
    hold on
    plot(aus_spec_lrc.wl_forw_f,smoothed_fmin_forw_lrc_mean(i,:),color=cmap(i,:))
    box on 
    xlabel("wl [nm]",FontSize=12)
    ylabel("Fmin_{f} \muW cm^{-2} nm^{-1}",FontSize=12)

    subplot(3,3,8)
    hold on
    plot(aus_spec_lrc.wl_back_f,smoothed_fs_back_lrc_mean(i,:),color=cmap(i,:))
    box on 
    xlabel("wl [nm]",FontSize=12)
    ylabel("Fs_b \muW cm^{-2} nm^{-1}",FontSize=12)

    subplot(3,3,9)
    hold on
    plot(aus_spec_lrc.wl_forw_f,smoothed_fs_forw_lrc_mean(i,:),color=cmap(i,:))
    box on 
    xlabel("wl [nm]",FontSize=12)
    ylabel("Fs_f \muW cm^{-2} nm^{-1}",FontSize=12)

end
c= colorbar('southoutside');
colormap(c,'jet')
caxis([min(aus_licor_co2.Ca(:,1)),max(aus_licor_co2.Ca(:,1))])
c.Label.String = "PAR [µmol m^{2} s^{-1}]";
newPosition = [0.35 0.05 0.4 0.015];
newUnits = 'normalized';
c.FontSize = 10;
c.Units = newUnits;
c.Position = newPosition;
sgtitle("LRC Experiment")
saveas(gcf,"./Figures/Aus_spec_lrc.png")
%%
% -----------------------------------------------------------------

% -----------Plot Fs in 757 (Magney 2017) for different Ca----------
figure(8)
set(gcf,'units','inch','Position',[1 2 5 5],'color','w')
I = find(aus_spec_lrc.wl_back_f>757 & aus_spec_lrc.wl_back_f<757.5);

h(1) = plot(mean_lrc_qin,smoothed_fs_back_lrc_mean(:,I),'r',LineWidth=1);
hold on
h(2) = plot(mean_lrc_qin,smoothed_fs_forw_lrc_mean(:,I),'k',LineWidth=1);

xlabel("PAR [µmol m^{2} s^{-1}]")
ylabel('Fs [\muW cm^{-2} nm^{-1}]')
hl = legend(h,'757nm backward','757nm forward')
hl.FontSize = 10;
saveas(gcf,"./Figures/Aus_fs_757_lrc.png")
%---------------------------------------------------------------------
%%
% --------------------- Plot PRI -------------------------------------
figure(9)
set(gcf,'units','inch','Position',[1 2 10 3],'color','w')
Legend=cell(length(aus_licor_lrc.Qin),1);
cmap = brewermap(length(names),"Dark2");
for i =1:length(names)
    subplot(1,3,1)
    hold on
    plot(mean_lrc_qin,pri_refl_lrc(:,i),".-",color=cmap(i,:))
    plot(mean_lrc_qin,pri_refl_lrc_mean,"-",color="black",LineWidth=2)
    box on 
    xlabel("PAR [µmol m^{2} s^{-1}]")
    ylabel('PRI')

    subplot(1,3,2)
    hold on
    plot(mean_lrc_qin,pri_trans_lrc(:,i),".-",color=cmap(i,:))
    plot(mean_lrc_qin,pri_trans_lrc_mean,"-",color="black",LineWidth=2)
    box on 
    xlabel("PAR [µmol m^{2} s^{-1}]")
    ylabel('PRI')

    subplot(1,3,3)
    hold on
    plot(mean_lrc_qin,pri_albedo_lrc(:,i),".-",color=cmap(i,:))
    plot(mean_lrc_qin,pri_albedo_lrc_mean,"-",color="black",LineWidth=2)
    box on 
    xlabel("PAR [µmol m^{2} s^{-1}]")
    ylabel('PRI')
end
sgtitle("LRC experiment")
saveas(gcf,"./Figures/Aus_PRI_lrc.png")
% -------------------------------------------------------------------
