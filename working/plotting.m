% load fluspect data
clear all 
clc
my_dir = '/home/hamid/SIF/sif-photo-integ/working/';
cd(my_dir)
datdir = './../data/output/fluspect_output/2015/2019-06-19-1452/';
addpath("/home/hamid/SIF/sif-photo-integ/")
out_dir = "./Figures/"

load([datdir 'spectral.mat']);
load([datdir 'measured.mat']);
load([datdir 'leafopt.mat']);
load([datdir 'leafbio.mat']);
load([datdir 'optipar.mat']);
I = csvread(['./../data/measured/Fluowat/2015/I.csv'],1); % Reading Irradiance data
I_mean = mean(I,2);   % Mean of measured E looks high! 
leafbio_std.Cab = 40;
leafbio_std.Cca = 10;
leafbio_std.Cdm = 0.012;
leafbio_std.Cw = 0.009;
leafbio_std.Cs = 0;
leafbio_std.N = 1.4;
leafbio_std.Cx = 0;
leafbio_std.Cant = 0;
n_rand = 100;

%Change E it in way PAR become ~0-2500[umol s-1 m-2]
E_range = linspace(0.001,0.80,100);
%% Run the model 

biovars = {'Cab','Cca','Cdm','Cw','Cs','N','Cx','Cant'};
n=length(biovars);

% generate n numbers for each fluspect parameter 
tmp.Cab = linspace(20,100,n_rand);
tmp.Cca = linspace(0,30,n_rand);
tmp.Cdm = linspace(0,0.5,n_rand);
tmp.Cw = linspace(0,0.4,n_rand);
tmp.Cs = linspace(0,0.6,n_rand);
tmp.N = linspace(1,4,n_rand);
tmp.Cx = linspace(0,1.5,n_rand);
tmp.Cant = linspace(0,10,n_rand);
% tmp.fqe1 = linspace(0.0001,0.2,n_rand);
% tmp.fqe = linspace(0.0001,0.2,n_rand);

for k =1:length(biovars)
    for i = 1:n_rand
        for j =1:length(E_range)
            E_scaled = I_mean*E_range(j); % This will replace the default par (Qin)
            E_par = calc_par(E_scaled,spectral.wlM);
            data.Qin = E_par;   % Its just a temporary PAR that will be replaced later
            data.Tin = 25;                 % Leaf temperature, C
            data.Cin = 200;                % Mesophyll CO2, ubar
            data.Oin = 209;                 % Atmospheric O2, mbar
            v = configure_fun(data);
            leafbio_std.(biovars{k}) = tmp.(biovars{k})(i);
            v.vcmax = (1.30 *leafbio_std.Cab+3.72)*1e-6;
            [out_opt out_photo] = main_fun_v1(spectral,leafbio_std,optipar,E_scaled,v);
            an(k,i,j) = out_photo.m.An_a;
            npq(k,i,j) = out_photo.m.PAM9_a;
            q(k,i,j) = out_photo.m.Q;
            a1(k,i,j) = out_photo.m.a1;
            a2(k,i,j) = out_photo.m.a2;
            cbf6_a(k,i,j) = out_photo.m.CB6F_a;
            jp680_1(k,i,j) = out_photo.m.JP680_a;
            jp700_1(k,i,j) = out_photo.m.JP700_a;
            cx(k,i,j) = out_photo.Cx;
            lue(k,i,j) = out_photo.m.An_a/out_photo.m.Q;
        end
    end
end
%%
labels = ['Chlorophyll (\mug cm^{2})',"Carotenoid (\mug cm^{2})",...
    'Dry matter content(g cm^{2})','water content (cm)','Brown pigment',...
    'Leaf mesophyll structure parameter','Photochemical reflectance parameter'...
    'Anthocyanin content (\mug cm^{2})'];
map = jet(n_rand);
for k=1:length(biovars)
    fig = figure (1)
    set(gcf,'units','inch','Position',[50 50 15 12],'color','w')
    for i = 1:n_rand
   
        s(1) = subplot(2,4,1);
        title('An')
        hold on
        plot(squeeze(q(k,i,:))*1e6, squeeze(an(k,i,:))*1e6,'color',map(i,:))
        box on
        xlabel('PAR (umol PPFD m-2 s-1)','FontSize', 11);
        ylabel('photosynthesis (umol CO2 m-2 s-1)','FontSize', 11);
            
        s(2) = subplot(2,4,2);
        title('NPQ')
        hold on
        plot(squeeze(q(k,i,:))*1e6, squeeze(npq(k,i,:)),'color',map(i,:))
        box on
        xlabel('PAR (umol PPFD m-2 s-1)','FontSize', 11);
        ylabel('NPQ','FontSize', 11);

        s(3) = subplot(2,4,3);
        title('PsI cross section')
        hold on
        plot(squeeze(q(k,i,:))*1e6, squeeze(a1(k,i,:)),'color',map(i,:))
        box on
        xlabel('PAR (umol PPFD m-2 s-1)','FontSize', 11);
        ylabel('PSI cross section','FontSize', 11);

        s(4) = subplot(2,4,4);
        title('PsII cross section')
        hold on
        plot(squeeze(q(k,i,:))*1e6, squeeze(a2(k,i,:)),'color',map(i,:))
        box on 
        xlabel('PAR (umol PPFD m-2 s-1)','FontSize', 11);
        ylabel('PSII cross section','FontSize', 11);

        s(5) =subplot(2,4,5);
        title('Closed PSI')
        hold on
        plot(squeeze(q(k,i,:))*1e6, squeeze(cbf6_a(k,i,:)),'color',map(i,:))
        box on 
        xlabel('PAR (umol PPFD m-2 s-1)','FontSize', 11);
        ylabel('closed PSI','FontSize', 11);

        s(6) =subplot(2,4,6);
        title('ETR PSI')
        hold on
        plot(squeeze(q(k,i,:))*1e6, squeeze(jp680_1(k,i,:)),'color',map(i,:))
        box on
        xlabel('PAR (umol PPFD m-2 s-1)','FontSize', 11);
        ylabel('ETR PSI','FontSize', 11);

        s(7) =subplot(2,4,7);
        title('ETR PSII')
        hold on
        plot(squeeze(q(k,i,:))*1e6, squeeze(jp700_1(k,i,:)),'color',map(i,:))
        box on 
        xlabel('PAR (umol PPFD m-2 s-1)','FontSize', 11);
        ylabel('ETR PSII','FontSize', 11);
        
        s(8) =subplot(2,4,8);
        title('LUE (An/PAR)')
        hold on
        plot(squeeze(q(k,i,:))*1e6, squeeze(lue(k,i,:)),'color',map(i,:))
        box on
        xlabel('PAR (umol PPFD m-2 s-1)','FontSize', 11);
        ylabel('LUE','FontSize', 11);

    end
h = axes(fig,'visible','off'); 
c = colorbar(h,'Position',[0.93 0.168 0.022 0.7]);  % attach colorbar to h
colormap(c,'jet')
caxis(h,[min(tmp.(biovars{k})) max(tmp.(biovars{k}))])
c.Label.String = labels(k);
c.FontSize = 12
saveas(gcf,strcat(out_dir,biovars(k),"_fig3.png") )
close all
end 

%%
% How photosynthesis changes 
leafbio_std.Cab = 40;
leafbio_std.Cca = 10;
leafbio_std.Cdm = 0.012;
leafbio_std.Cw = 0.009;
leafbio_std.Cs = 0;
leafbio_std.N = 1.4;
leafbio_std.Cx = 0;
leafbio_std.Cant = 0;

data.Qin = 1000;   % Its just a temporary PAR that will be replaced later
data.Tin = 25;                 % Leaf temperature, C
data.Cin = 200;                % Mesophyll CO2, ubar
data.Oin = 209;                 % Atmospheric O2, mbar
v = configure_fun(data);
E_range = linspace(0.001,0.80,10);


X_Labels = {'CB6F','RUB','Kf','Kd','Kp1','Kp2','Kn1'} ;
x_center = [v.CB6F,v.RUB,v.Kf,v.Kd,v.Kp1,v.Kp2,v.Kn1];
x_min = x_center-(x_center*0.1);
x_max = x_center+(x_center*0.1);

n_range = 10;
tmp.CB6F = linspace(x_min(1),x_max(1),n_range);  % Rubisco density
tmp.RUB = linspace(x_min(2),x_max(2),n_range);  % Rate constant for fluoresence at PSII and PSI
tmp.Kf = linspace(x_min(3),x_max(3),n_range);  % Rate constant for constitutive heat loss at PSII and PSI
tmp.Kd = linspace(x_min(4),x_max(4),n_range); % Rate constant for photochemistry at PSI
tmp.Kp1 = linspace(x_min(5),x_max(5),n_range); % Rate constant for photochemistry at PSI
tmp.Kp2 = linspace(x_min(6),x_max(6),n_range);   % Rate constant for regulated heat loss at PSI
tmp.Kn1 = linspace(x_min(7),x_max(7),n_range);  % Rate constant for regulated heat loss at PSI
% tmp.beta = linspace(x_min(8),x_max(8),n_range);   % PSII fraction of total leaf absorptance
% tmp.Ku = linspace(x_min(9),x_max(9),n_range);   % Rate constant for exciton sharing at PSII


I_531 = find(spectral.wlP==531);
I_570 = find(spectral.wlP==570);

for k = 1:length(X_Labels)
    for i =1:length(E_range)  %Loop over par
        for j = 1:n_range % Loop over variable ranges
            E_scaled = I_mean*E_range(i); 
            E_par = calc_par(E_scaled,spectral.wlM);
            data.Qin = E_par;   % Its just a temporary PAR that will be replaced later
            data.Tin = 25;                 % Leaf temperature, C
            data.Cin = 200;                % Mesophyll CO2, ubar
            data.Oin = 209;                 % Atmospheric O2, mbar
            v = configure_fun(data);
            v.vcmax = (1.30 *leafbio_std.Cab+3.72)*1e-6;
            v.(X_Labels{k}) =  tmp.(X_Labels{k})(j);
            [out_opt out_photo] = main_fun_v1(spectral,leafbio_std,optipar,E_scaled,v);
            out_fu(k,i,j,:) = out_opt.Fu;
            out_fd(k,i,j,:) = out_opt.Fd;
            out_refl(k,i,j,:) = out_opt.refl;
            out_tran(k,i,j,:) = out_opt.tran;
            b_531 = out_opt.refl(I_531);
            b_570 = leafopt_all.refl(I_570);
            PRI(k,i,j) = (b_531-b_570)/(b_531+b_570);
            
            x_wlf(k,i,j,:) = spectral.wlF;
            x_ref(k,i,j,:) = spectral.wlP;
            y_f(k,i,j,:) = E_par*ones(211,1);
            y_r(k,i,j,:) = E_par*ones(2001,1);
        end
    end
end


%%
% size(x(:))
% size(y)
% size(out_fu)
%plot3(x(:),y(:),out_fu(:))
%x = spectral.wlF;
%y = [1,2,3,4];
I_vis = find(spectral.wlP>=510&spectral.wlP<=570)
Labels = {'Cyt b6f density (mol sites m^{-2})','Rubisco density (mol sites m^{-2})',...
    'Rate constant for fluoresence at PSII and PSI (s^{-1})',...
    'Rate constant for constitutive heat loss at PSII and PSI (s^{-1})',...
    'Rate constant for photochemistry at PSI (s^{-1})',...
    'Rate constant for photochemistry at PSII (s^{-1})',...
    'Rate constant for regulated heat loss at PSI (s^{-1})'}

% {'CB6F','RUB','Kf','Kd','Kp1','Kp2','Kn1'} ;

for k = 1:length(X_Labels)
    fig = figure (2)
    set(gcf,'units','inch','Position',[50 50 15 12],'color','w')
    map = jet(n_range)
    s(1) = subplot(2,2,1)
    hold on
    for n = 1:length(E_range)
        for b = 1:n_range
            xx = squeeze(x_ref(k,n,b,I_vis));
            yy = squeeze(y_r(k,n,b,I_vis))';
            zz = squeeze(out_refl(k,n,b,I_vis));
            plot3(xx,yy,zz,'color',map(b,:))
        end
    end
    view(45,35)
    xlabel('wl (nm)','FontSize', 12)
    ylabel("PAR (umol PPFD m-2 s-1)",'FontSize', 12)
    zlabel("Reflectance",'FontSize', 12)
    
    s(2) = subplot(2,2,2)
    hold on
    for n = 1:length(E_range)
        for b = 1:n_range
            xx = squeeze(x_ref(k,n,b,I_vis));
            yy = squeeze(y_r(k,n,b,I_vis))';
            zz = squeeze(out_tran(k,n,b,I_vis));
            plot3(xx,yy,zz,'color',map(b,:))
        end
    end
    view(45,35)
    xlabel('wl (nm)','FontSize', 12)
    ylabel("PAR (umol PPFD m-2 s-1)",'FontSize', 12)
    zlabel("Tranmitance",'FontSize', 12)
    
    
    s(3) = subplot(2,2,3)
    hold on
    for n = 1:length(E_range)
        for b = 1:n_range
            xx = squeeze(x_wlf(k,n,b,:));
            yy = squeeze(y_f(k,n,b,:))';
            zz = squeeze(out_fu(k,n,b,:));
            plot3(xx,yy,zz,'color',map(b,:))
        end
    end
    view(45,35)
    xlabel('wl (nm)','FontSize', 12)
    ylabel("PAR (umol PPFD m-2 s-1)",'FontSize', 12)
    zlabel('L_Fu (W m^{-2}\mum^{-1}sr^{-1})','FontSize', 12)
    
    s(4) = subplot(2,2,4)
    hold on
    for n = 1:length(E_range)
        for b = 1:n_range
            xx = squeeze(x_wlf(k,n,b,:));
            yy = squeeze(y_f(k,n,b,:))';
            zz = squeeze(out_fd(k,n,b,:));
            plot3(xx,yy,zz,'color',map(b,:))
        end
    end
    view(45,35)
    xlabel('wl (nm)','FontSize', 12)
    ylabel("PAR (umol PPFD m-2 s-1)",'FontSize', 12)
    zlabel('L_Fd (W m^{-2}\mum^{-1}sr^{-1})','FontSize', 12)
    h = axes(fig,'visible','off'); 
    c = colorbar(h,'Position',[0.93 0.168 0.022 0.7]);  % attach colorbar to h
    colormap(c,'jet')
    caxis(h,[min(tmp.(X_Labels{k})) max(tmp.(X_Labels{k}))])
    c.Label.String = Labels(k);
    c.FontSize = 12
    saveas(gcf,strcat(out_dir,X_Labels(k),"_fig.png") )
    close all
end

%% How does PAM compare with Fu and Fd


% generate n numbers for each fluspect parameter 
tmp.Cab = linspace(0,100,n_rand);
E_range = linspace(0.001,0.80,10);

for i = 1:length(E_range)
    for j =1:n_rand
            E_scaled = I_mean*E_range(i); % This will replace the default par (Qin)
            E_par = calc_par(E_scaled,spectral.wlM);
            data.Qin = E_par;   % Its just a temporary PAR that will be replaced later
            data.Tin = 25;                 % Leaf temperature, C
            data.Cin = 200;                % Mesophyll CO2, ubar
            data.Oin = 209;                 % Atmospheric O2, mbar
            v = configure_fun(data);
            
            leafbio_std.Cab = tmp.Cab(j);
            [out_opt out_photo] = main_fun_v1(spectral,leafbio_std,optipar,E_scaled,v);
            an(i,j) = out_photo.m.An_a;
            npq(i,j) = out_photo.m.PAM9_a;
            q(i,j) = out_photo.m.Q;
            fm_a(i,j) = out_photo.m.Fm_a;
            fmp_a(i,j) = out_photo.m.Fmp_a;
            fs_a(i,j) = out_photo.m.Fs_a;
            fu(i,j,:) = out_opt.Fu;
            fd(i,j,:) = out_opt.Fd;
            ref(i,j,:) = out_opt.refl;
            tran(i,j,:) = out_opt.tran;
            cx(i,j) = out_photo.Cx;
            lue(i,j) = out_photo.m.An_a/out_photo.m.Q;
    end
end

%% plotting

fig = figure (3)
set(gcf,'units','inch','Position',[50 50 15 12],'color','w')
map = jet(length(E_range));
I_686 = find(spectral.wlF == 686);
I_740 = find(spectral.wlF == 740);
I_757 = find(spectral.wlF == 757);
I_771 = find(spectral.wlF == 771);
I_800 = find(spectral.wlF == 800);
for i = 1:5
    s(1) = subplot(2,3,1)
    hold on
    for k = 1:length(E_range)
        plot(squeeze(fmp_a(k,:)),fu(k,:,I_686),'color',map(k,:))
    end
    title('686 (nm)')
    xlabel('F^{`}_m','FontSize',12)
    ylabel('Fu (W m^{-2}\mum^{-1}sr^{-1})','FontSize', 12)

    s(2) = subplot(2,3,2)
    hold on
    for k = 1:length(E_range)
        plot(squeeze(fmp_a(k,:)),fu(k,:,I_740),'color',map(k,:))
    end
    title('740 (nm)')
    xlabel('F^{`}_m','FontSize',12)
    ylabel('Fu (W m^{-2}\mum^{-1}sr^{-1})','FontSize', 12)

    s(3) = subplot(2,3,3)
    hold on
    for k = 1:length(E_range)
        plot(squeeze(fmp_a(k,:)),fu(k,:,I_757),'color',map(k,:))
    end
    title('757 (nm)')
    xlabel('F^{`}_m','FontSize',12)
    ylabel('Fu (W m^{-2}\mum^{-1}sr^{-1})','FontSize', 12)
    
    s(4) = subplot(2,3,4)
    hold on
    for k = 1:length(E_range)
        plot(squeeze(fmp_a(k,:)),fu(k,:,I_771),'color',map(k,:))
    end
    title('771 (nm)')
    xlabel('F^{`}_m','FontSize',12)
    ylabel('Fu (W m^{-2}\mum^{-1}sr^{-1})','FontSize', 12)

    s(5) = subplot(2,3,5)
    hold on
    for k = 1:length(E_range)
        plot(squeeze(fmp_a(k,:)),fu(k,:,I_800),'color',map(k,:))
    end
    title('800 (nm)')
    xlabel('F^{`}_m','FontSize',12)
    ylabel('Fu (W m^{-2}\mum^{-1}sr^{-1})','FontSize', 12)

end
h = axes(fig,'visible','off'); 
c = colorbar(h,'Position',[0.93 0.168 0.022 0.7]);  % attach colorbar to h
colormap(c,'jet')
caxis(h,[0 2500])
c.Label.String = 'PAR (umol PPFD m-2 s-1)';
c.FontSize = 12
saveas(gcf,strcat(out_dir,"fmp_fu_fig.png") )
close all

fig = figure (3)
set(gcf,'units','inch','Position',[50 50 15 12],'color','w')
map = jet(length(E_range));
for i = 1:5
    s(1) = subplot(2,3,1)
    hold on
    for k = 1:length(E_range)
        plot(squeeze(fm_a(k,:)),fu(k,:,I_686),'color',map(k,:))
    end
    title('686 (nm)')
    xlabel('F_m','FontSize',12)
    ylabel('Fu (W m^{-2}\mum^{-1}sr^{-1})','FontSize', 12)

    s(2) = subplot(2,3,2)
    hold on
    for k = 1:length(E_range)
        plot(squeeze(fm_a(k,:)),fu(k,:,I_740),'color',map(k,:))
    end
    title('740 (nm)')
    xlabel('F_m','FontSize',12)
    ylabel('Fu (W m^{-2}\mum^{-1}sr^{-1})','FontSize', 12)

    s(3) = subplot(2,3,3)
    hold on
    for k = 1:length(E_range)
        plot(squeeze(fm_a(k,:)),fu(k,:,I_757),'color',map(k,:))
    end
    title('757 (nm)')
    xlabel('F_m','FontSize',12)
    ylabel('Fu (W m^{-2}\mum^{-1}sr^{-1})','FontSize', 12)
    
    s(4) = subplot(2,3,4)
    hold on
    for k = 1:length(E_range)
        plot(squeeze(fm_a(k,:)),fu(k,:,I_771),'color',map(k,:))
    end
    title('771 (nm)')
    xlabel('F_m','FontSize',12)
    ylabel('Fu (W m^{-2}\mum^{-1}sr^{-1})','FontSize', 12)

    s(5) = subplot(2,3,5)
    hold on
    for k = 1:length(E_range)
        plot(squeeze(fm_a(k,:)),fu(k,:,I_800),'color',map(k,:))
    end
    title('800 (nm)')
    xlabel('F_m','FontSize',12)
    ylabel('Fu (W m^{-2}\mum^{-1}sr^{-1})','FontSize', 12)

end
h = axes(fig,'visible','off'); 
c = colorbar(h,'Position',[0.93 0.168 0.022 0.7]);  % attach colorbar to h
colormap(c,'jet')
caxis(h,[0 2500])
c.Label.String = 'PAR (umol PPFD m-2 s-1)';
c.FontSize = 12
saveas(gcf,strcat(out_dir,"fm_fu_fig.png") )
close all


fig = figure (4)
set(gcf,'units','inch','Position',[50 50 15 12],'color','w')
map = jet(length(E_range));
for i = 1:5
    s(1) = subplot(2,3,1)
    hold on
    for k = 1:length(E_range)
        plot(squeeze(fs_a(k,:)),fu(k,:,I_686),'color',map(k,:))
    end
    title('686 (nm)')
    xlabel('F_s','FontSize',12)
    ylabel('Fu (W m^{-2}\mum^{-1}sr^{-1})','FontSize', 12)

    s(2) = subplot(2,3,2)
    hold on
    for k = 1:length(E_range)
        plot(squeeze(fs_a(k,:)),fu(k,:,I_740),'color',map(k,:))
    end
    title('740 (nm)')
    xlabel('F_s','FontSize',12)
    ylabel('Fu (W m^{-2}\mum^{-1}sr^{-1})','FontSize', 12)

    s(3) = subplot(2,3,3)
    hold on
    for k = 1:length(E_range)
        plot(squeeze(fs_a(k,:)),fu(k,:,I_757),'color',map(k,:))
    end
    title('757 (nm)')
    xlabel('F_s','FontSize',12)
    ylabel('Fu (W m^{-2}\mum^{-1}sr^{-1})','FontSize', 12)
    
    s(4) = subplot(2,3,4)
    hold on
    for k = 1:length(E_range)
        plot(squeeze(fs_a(k,:)),fu(k,:,I_771),'color',map(k,:))
    end
    title('771 (nm)')
    xlabel('F_s','FontSize',12)
    ylabel('Fu (W m^{-2}\mum^{-1}sr^{-1})','FontSize', 12)

    s(5) = subplot(2,3,5)
    hold on
    for k = 1:length(E_range)
        plot(squeeze(fs_a(k,:)),fu(k,:,I_800),'color',map(k,:))
    end
    title('800 (nm)')
    xlabel('F_s','FontSize',12)
    ylabel('Fu (W m^{-2}\mum^{-1}sr^{-1})','FontSize', 12)

end
h = axes(fig,'visible','off'); 
c = colorbar(h,'Position',[0.93 0.168 0.022 0.7]);  % attach colorbar to h
colormap(c,'jet')
caxis(h,[0 2500])
c.Label.String = 'PAR (umol PPFD m-2 s-1)';
c.FontSize = 12
saveas(gcf,strcat(out_dir,"fs_fu_fig.png") )
close all

fig = figure (5)
set(gcf,'units','inch','Position',[50 50 15 12],'color','w')
hold on
for i = 1:length(E_range)
plot(tmp.Cab, an(i,:),'color',map(i,:))
end
xlabel('Chlorophyll (\mug cm^{2})')
ylabel('An (umol CO2 m-2 s-1)')
h = axes(fig,'visible','off'); 
c = colorbar(h,'Position',[0.93 0.168 0.022 0.7]);  % attach colorbar to h
colormap(c,'jet')
caxis(h,[0 2500])
c.Label.String = 'PAR (umol PPFD m-2 s-1)';
c.FontSize = 12
saveas(gcf,strcat(out_dir,"chl_an_fig.png") )
close all

fig = figure (6)
set(gcf,'units','inch','Position',[50 50 15 12],'color','w')
hold on
for i = 1:length(E_range)
plot(tmp.Cab, lue(i,:),'color',map(i,:))
end
xlabel('Chlorophyll (\mug cm^{2})')
ylabel('LUE')
h = axes(fig,'visible','off'); 
c = colorbar(h,'Position',[0.93 0.168 0.022 0.7]);  % attach colorbar to h
colormap(c,'jet')
caxis(h,[0 2500])
c.Label.String = 'PAR (umol PPFD m-2 s-1)';
c.FontSize = 12
saveas(gcf,strcat(out_dir,"chl_lue_fig.png") )
close all