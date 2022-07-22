clear 
clc

licor_dir =  '/home/hamid/SIF/sif-photo-integ/data/william_woodgate/LI6800/'
asd_dir =  '/home/hamid/SIF/sif-photo-integ/data/william_woodgate/ASD/'
qep_dir =   '/home/hamid/SIF/sif-photo-integ/data/william_woodgate/QEP/'

% Read data for co2 experiment


% photosynthesis and PAM data

names = ["L1", "L2", "L3","L4","L5","L6","L7","L8","L9"];

% Find the definition and units here:https://www.licor.com/env/support/LI-6800/topics/symbols.html
vars_licor = {"Ca","Ci", "A","gsw","gbw","gtc","Rabs","TleafEB","TleafCnd",...
    "TleafCnd","VPDleaf","E","Qabs","Qin","Fs","Fm_","PhiPS2","PS2_1",...
    "Fo","Fm","Fv","Fv_Fm","ETR","Fv__Fm_","PhiCO2","NPQ","Fo_",...
    "Fv_","qP","qN","qP_Fo","qN_Fo","qL","x1_qL","Qin","Qabs","alpha",...
    "Tair","Tleaf",};
vars_qep = {"refl","trans","fmax_back","fmax_forw","fmin_back","fmin_forw",...
    "fs_back","fs_forw"};
green_wl_low = 500;
green_wl_high = 640;

f_wl_low = 660;
f_wl_high = 800;
%% The co2 experiment

for k =1:length(names)
    if k==9
        % There is no CO2 data for this leaf
        continue
    end
    data_co2_licor = readtable(strcat(licor_dir,names(k),"_CO2.csv"),'NumHeaderLines',1); 
    data_co2_licor(1,:) = [];    % Remove unit row

    for i =1:length(vars_licor)
        if k==6
            licor_co2.(vars_licor{i})(:,k) = data_co2_licor.(vars_licor{i})(1:17);
        else
            licor_co2.(vars_licor{i})(:,k) = data_co2_licor.(vars_licor{i})(2:18);
        end
        
    end
end

% Reflectance and Transmittance

for k = 1:length(names)
    
    if k==9
      continue
    end
    
    for i = 1:length(vars_qep)
        data_co2_spec = readtable(strcat(qep_dir,names(k),"_CO2_",(vars_qep{i}),".csv")); 
        data_co2_spec(:,1) = [];  %Remove time column
        
        data_co2_licor = readtable(strcat(licor_dir,names(k),"_CO2.csv"),'NumHeaderLines',1); 
        data_co2_licor(1,:) = [];    % Remove unit row
        
        wl = table2array(data_co2_spec(1,:));
        wl_green = find(wl>=green_wl_low & wl<=green_wl_high);
        wl_f = find(wl>=f_wl_low & wl<=f_wl_high);
        data_co2_spec(1,:) = [];  %Remove wl rows
        
        if size(data_co2_spec,1) ~= length(data_co2_licor.Ca)
            disp("ERROR Check co2 levels")
            continue
        end
        
        if (i==1) | (i==2)

        % Check we pick right vars corresponding to diff Ca
            
            if k==6
                spec_co2.(vars_qep{i})(:,:,k) = table2array(data_co2_spec(1:17,wl_green));
    
            else
                
                spec_co2.(vars_qep{i})(:,:,k) = table2array(data_co2_spec(2:18,wl_green));
            end
            
            
        else
            if k==6
                spec_co2.(vars_qep{i})(:,:,k) = table2array(data_co2_spec(1:17,wl_f));
    
            else
                
                spec_co2.(vars_qep{i})(:,:,k) = table2array(data_co2_spec(2:18,wl_f));
                
            end
        end 
    end 

end
% For some reasons the wl recorded for the refl and trans are not the same
data_co2_spec = readtable(strcat(qep_dir,names(1),"_CO2_refl.csv")); 
data_co2_spec(:,1) = [];  %Remove time column        
wl_reflectance = table2array(data_co2_spec(1,:));
wl_refl_green = wl_reflectance(find(wl_reflectance>=500 & wl_reflectance<=640));
wl_back_f = wl_reflectance(find(wl_reflectance>=660 & wl_reflectance<=800));

data_co2_spec = readtable(strcat(qep_dir,names(1),"_CO2_trans.csv")); 
data_co2_spec(:,1) = [];  %Remove time column        
wl_transmittance = table2array(data_co2_spec(1,:));
wl_trans_green= wl_transmittance(find(wl_transmittance>=500 & wl_transmittance<=640));
wl_forw_f = wl_transmittance(find(wl_transmittance>=660 & wl_transmittance<=800));

spec_co2.wl_refl_green = wl_refl_green;
spec_co2.wl_back_f = wl_back_f;

spec_co2.wl_trans_green = wl_trans_green;
spec_co2.wl_forw_f = wl_forw_f;

%% The light curve data

for k =1:length(names)

    data_lrc_licor = readtable(strcat(licor_dir,names(k),"_LRC.csv"),'NumHeaderLines',1); 
    data_lrc_licor(1,:) = [];    % Remove unit row
    
    for i =1:length(vars_licor)
        if (k<=3)
            licor_lrc.(vars_licor{i})(:,k) = data_lrc_licor.(vars_licor{i})(2:17);
        else
            licor_lrc.(vars_licor{i})(:,k) = data_lrc_licor.(vars_licor{i})(1:16);
        end
        
    end
end



for k = 1:length(names)
    
    for i = 1:length(vars_qep)
        data_lrc_spec = readtable(strcat(qep_dir,names(k),"_LRC_",(vars_qep{i}),".csv")); 
        data_lrc_spec(:,1) = [];  %Remove time column
        
        data_lrc_licor = readtable(strcat(licor_dir,names(k),"_LRC.csv"),'NumHeaderLines',1); 
        data_lrc_licor(1,:) = [];    % Remove unit row
        
        wl = table2array(data_lrc_spec(1,:));
        wl_green = find(wl>=green_wl_low & wl<=green_wl_high);
        wl_f = find(wl>=f_wl_low & wl<=f_wl_high);
        data_lrc_spec(1,:) = [];  %Remove wl rows
        
        if size(data_lrc_spec,1) ~= length(data_lrc_licor.Qin)
            disp("ERROR Check co2 levels")
            continue
        end
        
        if (i==1) | (i==2)

        % Check we pick right vars corresponding to diff Ca
            
            if k<=3
                spec_lrc.(vars_qep{i})(:,:,k) = table2array(data_lrc_spec(2:17,wl_green));
    
            else
                
                spec_lrc.(vars_qep{i})(:,:,k) = table2array(data_lrc_spec(1:16,wl_green));
            end
            
            
        else
            if k<=3
                spec_lrc.(vars_qep{i})(:,:,k) = table2array(data_lrc_spec(2:17,wl_f));
    
            else
                
                spec_lrc.(vars_qep{i})(:,:,k) = table2array(data_lrc_spec(1:16,wl_f));
                
            end
        end 
    end 

end

spec_lrc.wl_refl_green = wl_refl_green;
spec_lrc.wl_back_f = wl_back_f;

spec_lrc.wl_trans_green = wl_trans_green;
spec_lrc.wl_forw_f = wl_forw_f;

%% Plotting co2 experiment
figure(1)
set(gcf,'units','inch','Position',[25 25 12 5],'color','w')
cmap = brewermap(length(names)-1,"Set1")
Legend=cell(length(names)-1,1)
for i =1:length(names)-1
    
    s(1) = subplot(1,2,1)
    hold on
    plot(licor_co2.Ca(:,i),licor_co2.A(:,i),".-",color=cmap(i,:))
    box on
    xlabel("Ambient CO_{2} [µmol mol^{-1}]",FontSize=12)
    ylabel("A_{n} [µmol m^{-2} s^{-1}]",FontSize=12)
    
    s(2) = subplot(1,2,2)
    hold on
    plot(licor_co2.Ca(:,i),licor_co2.gsw(:,i),".-",color=cmap(i,:))
    box on
    xlabel("Ambient CO_{2} [µmol mol^{-1}]",FontSize=12)
    ylabel("gsw [mol m^{-2} s^{-1}]",FontSize=12)

    Legend{i} = names(i);
end
sgtitle("CO_{2} experiment")
legend(s(1),Legend,"location","southeast")
saveas(gcf,"./figures/GasExc_co2.png")
%% Plotting refl/trans for the co2 experiment
% Smoothing refl/trans and calculate PRI 
smoothed_refl_co2 = zeros(length(licor_co2.Ca),length(spec_co2.wl_refl_green),...
    length(names)-1);
smoothed_trans_co2 = zeros(length(licor_co2.Ca),length(spec_co2.wl_trans_green),...
    length(names)-1);

for i =1:17
    for j=1:8
        smoothed_refl_co2(i,:,j) = sgolayfilt(spec_co2.refl(i,:,j),2,7);
        smoothed_trans_co2(i,:,j) = sgolayfilt(spec_co2.trans(i,:,j),2,7);
    end
end

wl_refl_531 = find(spec_co2.wl_refl_green>=531 & spec_co2.wl_refl_green<532);
wl_refl_570 = find(spec_co2.wl_refl_green>=570 & spec_co2.wl_refl_green<571);
wl_trans_531 = find(spec_co2.wl_trans_green>=531 & spec_co2.wl_trans_green<532);
wl_trans_570 = find(spec_co2.wl_trans_green>=570 & spec_co2.wl_trans_green<571);

a_refl = squeeze(spec_co2.refl(:,wl_refl_531,:));
b_refl = squeeze(spec_co2.refl(:,wl_refl_570,:));
pri_refl_co2 = (a_refl-b_refl)./(a_refl+b_refl);

a_trans = squeeze(spec_co2.trans(:,wl_trans_531,:));
b_trans = squeeze(spec_co2.trans(:,wl_trans_570,:));
pri_trans_co2 = (a_trans-b_trans)./(a_trans+b_trans);

% Take the mean over 8 observations
smoothed_refl_mean = mean(smoothed_refl_co2,3);
smoothed_trans_mean = mean(smoothed_trans_co2,3);

figure(2)
set(gcf,'units','inch','Position',[1 2 18 4],'color','w')
Legend=cell(length(licor_co2.Ca),1)
cmap = jet(17);
for i=1:17
    subplot(1,3,1)
    hold on
    plot(spec_co2.wl_refl_green,smoothed_refl_mean(i,:),color=cmap(i,:))
    box on 
    xlabel("wl",FontSize=12)
    ylabel("Reflectance",FontSize=12)

    subplot(1,3,2)
    hold on
    plot(spec_co2.wl_trans_green,smoothed_trans_mean(i,:),color=cmap(i,:))
    box on
    xlabel("wl",FontSize=12)
    ylabel("Transmittance",FontSize=12)

    subplot(1,3,3)
    hold on
    plot(spec_co2.wl_trans_green,smoothed_refl_mean(i,:)+smoothed_trans_mean(i,:),color=cmap(i,:))
    box on
    xlabel("wl",FontSize=12)
    ylabel("Leaf albedo",FontSize=12)
end
c= colorbar
colormap(c,'jet')
caxis([min(licor_co2.Ca(:,1)),max(licor_co2.Ca(:,1))])
c.Label.String = "Ambient CO_{2} [µmol mol^{-1}]";
c.FontSize = 12
sgtitle("CO_{2} experiment")
saveas(gcf,"./figures/spec_co2.png")

%% Plotting PAR experiment
figure(3)
set(gcf,'units','inch','Position',[25 25 12 5],'color','w')
cmap = brewermap(length(names),"Set1")
Legend=cell(length(names),1)
for i =1:length(names)
    
    s(1) = subplot(1,2,1)
    hold on
    plot(licor_lrc.Qin(:,i),licor_lrc.A(:,i),".-",color=cmap(i,:))
    box on
    xlabel("Incoming PAR [µmol m^{-2} s^{-1}]",FontSize=12)
    ylabel("A_{n} [µmol m^{-2} s^{-1}]",FontSize=12)
    
    s(2) = subplot(1,2,2)
    hold on
    plot(licor_lrc.Qin(:,i),licor_lrc.gsw(:,i),".-",color=cmap(i,:))
    box on
    xlabel("Incoming PAR [µmol m^{-2} s^{-1}]",FontSize=12)
    ylabel("gsw [mol m^{-2} s^{-1}]",FontSize=12)

    Legend{i} = names(i);
end
sgtitle("LRC experiment")
legend(s(1),Legend,"location","southeast")
sgtitle("LRC experiment")
saveas(gcf,"./figures/GasExc_lrc.png")

%% Plotting refl/trans for the co2 experiment
% Smoothing refl/trans and calculate PRI 
smoothed_refl_lrc = zeros(length(licor_lrc.Qin),length(spec_lrc.wl_refl_green),...
    length(names));
smoothed_trans_lrc = zeros(length(licor_lrc.Qin),length(spec_lrc.wl_trans_green),...
    length(names));
for i =1:length(licor_lrc.Qin)
    for j=1:length(names)
        smoothed_refl_lrc(i,:,j) = sgolayfilt(spec_lrc.refl(i,:,j),2,7);
        smoothed_trans_lrc(i,:,j) = sgolayfilt(spec_lrc.trans(i,:,j),2,7);
    end
end


wl_refl_531 = find(spec_lrc.wl_refl_green>=531 & spec_lrc.wl_refl_green<532);
wl_refl_570 = find(spec_lrc.wl_refl_green>=570 & spec_lrc.wl_refl_green<571);
wl_trans_531 = find(spec_lrc.wl_trans_green>=531 & spec_lrc.wl_trans_green<532);
wl_trans_570 = find(spec_lrc.wl_trans_green>=570 & spec_lrc.wl_trans_green<571);

a_refl = squeeze(spec_lrc.refl(:,wl_refl_531,:));
b_refl = squeeze(spec_lrc.refl(:,wl_refl_570,:));
pri_refl_lrc = (a_refl-b_refl)./(a_refl+b_refl);

a_trans = squeeze(spec_lrc.trans(:,wl_trans_531,:));
b_trans = squeeze(spec_lrc.trans(:,wl_trans_570,:));
pri_trans_lrc = (a_trans-b_trans)./(a_trans+b_trans);

% Take the mean over 8 observations
smoothed_refl_lrc_mean = mean(smoothed_refl_lrc,3);
smoothed_trans_lrc_mean = mean(smoothed_trans_lrc,3);

figure(4)
set(gcf,'units','inch','Position',[1 2 18 4],'color','w')
Legend=cell(length(licor_lrc.Qin),1)
cmap = jet(length(licor_lrc.Qin));
for i=1:length(licor_lrc.Qin)
    subplot(1,3,1)
    hold on
    plot(spec_lrc.wl_refl_green,smoothed_refl_lrc_mean(i,:),color=cmap(i,:))
    box on 
    xlabel("wl",FontSize=12)
    ylabel("Reflectance",FontSize=12)

    subplot(1,3,2)
    hold on
    plot(spec_lrc.wl_trans_green,smoothed_trans_lrc_mean(i,:),color=cmap(i,:))
    box on
    xlabel("wl",FontSize=12)
    ylabel("Transmittance",FontSize=12)

    subplot(1,3,3)
    hold on
    plot(spec_lrc.wl_trans_green,smoothed_refl_lrc_mean(i,:)+smoothed_trans_lrc_mean(i,:),color=cmap(i,:))
    box on
    xlabel("wl",FontSize=12)
    ylabel("Leaf albedo",FontSize=12)
end
c= colorbar
colormap(c,'jet')
caxis([min(licor_lrc.Qin(:,1)),max(licor_lrc.Qin(:,1))])
c.Label.String = "Ambient CO_{2} [µmol mol^{-1}]";
c.FontSize = 12;
sgtitle("LRC experiment");
saveas(gcf,"./figures/spec_lrc.png")

%% Plot PRI 
mean_co2_ca =mean(licor_co2.Ca,2);
mean_lrc_qin = mean(licor_lrc.Qin,2);
pri_refl_co2_mean = mean(pri_refl_co2,2);
pri_trans_co2_mean = mean(pri_trans_co2,2);

pri_refl_lrc_mean = mean(pri_refl_lrc,2);
pri_trans_lrc_mean = mean(pri_trans_lrc,2);


figure(5)
set(gcf,'units','inch','Position',[25 25 12 5],'color','w')
cmap = brewermap(length(names),"Set1")
Legend=cell(length(names),1)

for i =1:length(names)
    if i<9
        subplot(2,2,1)
        hold on
        plot(mean_co2_ca,pri_refl_co2(:,i),".-",color=cmap(i,:))
        plot(mean_co2_ca,pri_refl_co2_mean,"-",color="black",LineWidth=2)
        title('PRI Co2 experiment (refl)')
        xlabel('wl')
        ylabel('PRI')
        
        subplot(2,2,2)
        hold on
        plot(mean_lrc_qin,pri_refl_lrc(:,i),".-",color=cmap(i,:))
        plot(mean_lrc_qin,pri_refl_lrc_mean,"-",color="black",LineWidth=2)
        title('PRI LRC experiment (refl)')
        xlabel('wl')
        ylabel('PRI')

        subplot(2,2,3)
        hold on
        plot(mean_co2_ca,pri_trans_co2(:,i),".-",color=cmap(i,:))
        plot(mean_co2_ca,pri_trans_co2_mean,"-",color="black",LineWidth=2)
        title('PRI Co2 experiment (trans)')
        xlabel('wl')
        ylabel('PRI')

        subplot(2,2,4)
        hold on
        plot(mean_lrc_qin,pri_trans_lrc(:,i),".-",color=cmap(i,:))
        plot(mean_lrc_qin,pri_trans_lrc_mean,"-",color="black",LineWidth=2)
        title('PRI LRC experiment (trans)')
        xlabel('wl')
        ylabel('PRI')

    end
end

sgtitle("PRI");
saveas(gcf,"./figures/PRI.png")

%% Plotting SIF for the co2 experiment
% Smoothing SIF 

smoothed_fmax_back_co2 = zeros(length(licor_co2.Ca),length(spec_co2.wl_back_f),...
    length(names)-1);
smoothed_fmax_forw_co2 = zeros(length(licor_co2.Ca),length(spec_co2.wl_forw_f),...
    length(names)-1);

smoothed_fmin_back_co2 = zeros(length(licor_co2.Ca),length(spec_co2.wl_back_f),...
    length(names)-1);
smoothed_fmin_forw_co2 = zeros(length(licor_co2.Ca),length(spec_co2.wl_forw_f),...
    length(names)-1);

smoothed_fs_back_co2 = zeros(length(licor_co2.Ca),length(spec_co2.wl_back_f),...
    length(names)-1);
smoothed_fs_forw_co2 = zeros(length(licor_co2.Ca),length(spec_co2.wl_forw_f),...
    length(names)-1);

for i =1:length(licor_co2.Ca)
    for j=1:length(names)-1
        smoothed_fmax_back_co2(i,:,j) = sgolayfilt(spec_co2.fmax_back(i,:,j),2,7);
        smoothed_fmax_forw_co2(i,:,j) = sgolayfilt(spec_co2.fmax_forw(i,:,j),2,7);

        smoothed_fmin_back_co2(i,:,j) = sgolayfilt(spec_co2.fmin_back(i,:,j),2,7);
        smoothed_fmin_forw_co2(i,:,j) = sgolayfilt(spec_co2.fmin_forw(i,:,j),2,7);

        smoothed_fs_back_co2(i,:,j) = sgolayfilt(spec_co2.fs_back(i,:,j),2,7);
        smoothed_fs_forw_co2(i,:,j) = sgolayfilt(spec_co2.fs_forw(i,:,j),2,7);


    end
end

% Take the mean over 8 observations
smoothed_fmax_back_co2_mean = mean(smoothed_fmax_back_co2,3);
smoothed_fmax_forw_co2_mean = mean(smoothed_fmax_forw_co2,3);

smoothed_fmin_back_co2_mean = mean(smoothed_fmin_back_co2,3);
smoothed_fmin_forw_co2_mean = mean(smoothed_fmin_forw_co2,3);

smoothed_fs_back_co2_mean = mean(smoothed_fs_back_co2,3);
smoothed_fs_forw_co2_mean = mean(smoothed_fs_forw_co2,3);

figure(2)
set(gcf,'units','inch','Position',[1 2 18 8],'color','w')
Legend=cell(length(licor_co2.Ca),1)
cmap = jet(length(licor_co2.Ca));
for i=1:17
    subplot(3,2,1)
    hold on
    plot(spec_co2.wl_back_f,smoothed_fmax_back_co2_mean(i,:),color=cmap(i,:))
    box on 
    xlabel("wl",FontSize=12)
    ylabel("Fmax bs \muW cm^{-2} nm^{-1}",FontSize=12)

    subplot(3,2,2)
    hold on
    plot(spec_co2.wl_forw_f,smoothed_fmax_forw_co2_mean(i,:),color=cmap(i,:))
    box on 
    xlabel("wl",FontSize=12)
    ylabel("Fmax fs \muW cm^{-2} nm^{-1}",FontSize=12)

    subplot(3,2,3)
    hold on
    plot(spec_co2.wl_back_f,smoothed_fmin_back_co2_mean(i,:),color=cmap(i,:))
    box on 
    xlabel("wl",FontSize=12)
    ylabel("Fmin bs \muW cm^{-2} nm^{-1}",FontSize=12)

    subplot(3,2,4)
    hold on
    plot(spec_co2.wl_forw_f,smoothed_fmin_forw_co2_mean(i,:),color=cmap(i,:))
    box on 
    xlabel("wl",FontSize=12)
    ylabel("Fmin fs \muW cm^{-2} nm^{-1}",FontSize=12)

    subplot(3,2,5)
    hold on
    plot(spec_co2.wl_back_f,smoothed_fs_back_co2_mean(i,:),color=cmap(i,:))
    box on 
    xlabel("wl",FontSize=12)
    ylabel("Fs bs \muW cm^{-2} nm^{-1}",FontSize=12)

    subplot(3,2,6)
    hold on
    plot(spec_co2.wl_forw_f,smoothed_fs_forw_co2_mean(i,:),color=cmap(i,:))
    box on 
    xlabel("wl",FontSize=12)
    ylabel("Fs fs \muW cm^{-2} nm^{-1}",FontSize=12)

end
c= colorbar
colormap(c,'jet')
caxis([min(licor_co2.Ca(:,1)),max(licor_co2.Ca(:,1))])
c.Label.String = "Ambient CO_{2} [µmol mol^{-1}]";
c.FontSize = 12
sgtitle("CO_{2} experiment")
saveas(gcf,"./figures/F_co2.png")

%% Plotting SIF for the co2 experiment
% Smoothing SIF 

smoothed_fmax_back_lrc = zeros(length(licor_lrc.Qin),length(spec_lrc.wl_back_f),...
    length(names));
smoothed_fmax_forw_lrc = zeros(length(licor_lrc.Qin),length(spec_lrc.wl_forw_f),...
    length(names));

smoothed_fmin_back_lrc = zeros(length(licor_lrc.Qin),length(spec_lrc.wl_back_f),...
    length(names));
smoothed_fmin_forw_lrc = zeros(length(licor_lrc.Qin),length(spec_co2.wl_forw_f),...
    length(names));

smoothed_fs_back_lrc = zeros(length(licor_lrc.Qin),length(spec_lrc.wl_back_f),...
    length(names));
smoothed_fs_forw_lrc = zeros(length(licor_lrc.Qin),length(spec_lrc.wl_forw_f),...
    length(names));

for i =1:length(licor_lrc.Qin)
    for j=1:length(names)
        smoothed_fmax_back_lrc(i,:,j) = sgolayfilt(spec_lrc.fmax_back(i,:,j),2,7);
        smoothed_fmax_forw_lrc(i,:,j) = sgolayfilt(spec_lrc.fmax_forw(i,:,j),2,7);

        smoothed_fmin_back_lrc(i,:,j) = sgolayfilt(spec_lrc.fmin_back(i,:,j),2,7);
        smoothed_fmin_forw_lrc(i,:,j) = sgolayfilt(spec_lrc.fmin_forw(i,:,j),2,7);

        smoothed_fs_back_lrc(i,:,j) = sgolayfilt(spec_lrc.fs_back(i,:,j),2,7);
        smoothed_fs_forw_lrc(i,:,j) = sgolayfilt(spec_lrc.fs_forw(i,:,j),2,7);


    end
end

% Take the mean over 8 observations
smoothed_fmax_back_lrc_mean = mean(smoothed_fmax_back_lrc,3);
smoothed_fmax_forw_lrc_mean = mean(smoothed_fmax_forw_lrc,3);

smoothed_fmin_back_lrc_mean = mean(smoothed_fmin_back_lrc,3);
smoothed_fmin_forw_lrc_mean = mean(smoothed_fmin_forw_lrc,3);

smoothed_fs_back_lrc_mean = mean(smoothed_fs_back_lrc,3);
smoothed_fs_forw_lrc_mean = mean(smoothed_fs_forw_lrc,3);

figure(2)
set(gcf,'units','inch','Position',[1 2 18 8],'color','w')
Legend=cell(length(licor_lrc.Qin),1)
cmap = jet(length(licor_lrc.Qin));
for i=1:length(licor_lrc.Qin)
    subplot(3,2,1)
    hold on
    plot(spec_lrc.wl_back_f,smoothed_fmax_back_lrc_mean(i,:),color=cmap(i,:))
    box on 
    xlabel("wl",FontSize=12)
    ylabel("Fmax bs \muW cm^{-2} nm^{-1}",FontSize=12)

    subplot(3,2,2)
    hold on
    plot(spec_lrc.wl_forw_f,smoothed_fmax_forw_lrc_mean(i,:),color=cmap(i,:))
    box on 
    xlabel("wl",FontSize=12)
    ylabel("Fmax fs \muW cm^{-2} nm^{-1}",FontSize=12)

    subplot(3,2,3)
    hold on
    plot(spec_lrc.wl_back_f,smoothed_fmin_back_lrc_mean(i,:),color=cmap(i,:))
    box on 
    xlabel("wl",FontSize=12)
    ylabel("Fmin bs \muW cm^{-2} nm^{-1}",FontSize=12)

    subplot(3,2,4)
    hold on
    plot(spec_lrc.wl_forw_f,smoothed_fmin_forw_lrc_mean(i,:),color=cmap(i,:))
    box on 
    xlabel("wl",FontSize=12)
    ylabel("Fmin fs \muW cm^{-2} nm^{-1}",FontSize=12)

    subplot(3,2,5)
    hold on
    plot(spec_lrc.wl_back_f,smoothed_fs_back_lrc_mean(i,:),color=cmap(i,:))
    box on 
    xlabel("wl",FontSize=12)
    ylabel("Fs bs \muW cm^{-2} nm^{-1}",FontSize=12)

    subplot(3,2,6)
    hold on
    plot(spec_lrc.wl_forw_f,smoothed_fs_forw_lrc_mean(i,:),color=cmap(i,:))
    box on 
    xlabel("wl",FontSize=12)
    ylabel("Fs fs \muW cm^{-2} nm^{-1}",FontSize=12)

end
c= colorbar
colormap(c,'jet')
caxis([min(licor_co2.Ca(:,1)),max(licor_co2.Ca(:,1))])
c.Label.String = "Ambient CO_{2} [µmol mol^{-1}]";
c.FontSize = 12
sgtitle("LRC experiment")
saveas(gcf,"./figures/f_LRC.png")





