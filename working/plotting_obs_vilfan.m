clear
clc

load obs
data = obs.vilfan_data;
data_co2 = data.data_co2;
data_lrc = data.data_lrc;
leaf_id_co2 = data_co2.leaf_id;
leaf_id_lrc = data_lrc.leaf_id;
%% plot CO2 experiment, licor data
figure(1)
set(gcf,'units','inch','Position',[5 5 9 12],'color','w')
licor = data_co2.licor_pam;
names = licor.SampleName;
leaf_IDs = unique(data_co2.leaf_id); 
n_sample = length(leaf_IDs);
% cmap = brewermap(n_sample,"Dark2")
Legend=cell(n_sample,1)
Fsize=11;
mycolor = [.7 .7 .7];

for i =1:n_sample
    I = data_co2.leaf_id==i;
    Ca = licor.CO2(I);
    An = licor.Photo(I);
    Fs = licor.Ft(I); % Note in CSV file it is called Ft
    Fmp = licor.Fmp(I);
    Fm = licor.Fm(I);
    Fo = licor.Fo(I);
    pari = licor.PARi(I);

    % Calculate PAM indices
    Fs_Fo = Fs./Fo;
    Fm_Fo = Fm./Fo;
    Fop = 1./((1./Fo)-(1./Fm)+(1./Fmp));
    Fmp_F0 = Fmp./Fo;
    Fv_Fm = 1-(Fo./Fm);
    PhiPS2 = 1-(Fs./Fmp);
    qP = (Fmp-Fs)./(Fmp-Fop);
    qL = qP.*(Fop./Fs);
    qN = (Fm-Fmp)./(Fm-Fop);
    NPQ = (Fm./Fmp)-1;
    Y_NO = Fs./Fm;
    Y_NPQ = (Fs./Fmp)-(Fs./Fm);
    % ETR calculation is based on miniPAM manual which is used by Vilfan
    ETR = pari.*0.84.*0.5.*PhiPS2;
    qL_1 = 1-qL;


    s(1) = subplot(5,3,1)
    hold on
    plot(Ca,An ,".-",color=mycolor)
    box on
    ylabel("A_{n} [µmol m^{-2} s^{-1}]",FontSize=Fsize)
    ylim([-10 33])

    s(2) = subplot(5,3,2)
    hold on
    plot(Ca,Fs_Fo,".-",color=mycolor)
    box on
    ylabel("Fs",FontSize=Fsize)
    ylim([0.5 2.1])

    s(3) = subplot(5,3,3)
    hold on
    plot(Ca,Fm_Fo,".-",color=mycolor)
    box on
    ylabel("Fm",FontSize=Fsize)
    ylim([4.5 6.6])
        
    s(4) = subplot(5,3,4)
    hold on
    plot(Ca,Fmp_F0,".-",color=mycolor)
    box on
    ylabel("Fm'",FontSize=Fsize)

    s(5) = subplot(5,3,5)
    hold on
    plot(Ca,Fv_Fm,".-",color=mycolor)
    box on
    ylabel("Fv/Fm",FontSize=Fsize)
    ylim([0.79 0.85]) 

    s(6) = subplot(5,3,6)
    hold on
    plot(Ca,PhiPS2,".-",color=mycolor)
    box on
    ylabel("PhiPS2",FontSize=Fsize)
    ylim([0.1 0.7])

    s(7) = subplot(5,3,7)
    hold on
    plot(Ca,qP,".-",color=mycolor)
    box on
    ylabel("qP",FontSize=Fsize)

    s(8) = subplot(5,3,8)
    hold on
    plot(Ca,qL,".-",color=mycolor)
    box on
    ylabel("qL",FontSize=Fsize)

    s(9) = subplot(5,3,9)
    hold on
    plot(Ca,qN,".-",color=mycolor)
    box on
    ylabel("qN",FontSize=Fsize)

    s(10) = subplot(5,3,10)
    hold on
    plot(Ca,NPQ,".-",color=mycolor)
    box on
    ylabel("NPQ",FontSize=Fsize)

    s(11) = subplot(5,3,11)
    hold on
    plot(Ca,Y_NO,".-",color=mycolor)
    box on
    ylabel("Y(NO)",FontSize=Fsize)

    s(12) = subplot(5,3,12)
    hold on
    plot(Ca,Y_NPQ,".-",color=mycolor)
    box on
    ylabel("Y(NPQ)",FontSize=Fsize)
    xlabel("Ambient CO_{2} [µmol mol^{-1}]",FontSize=Fsize)

    s(13) = subplot(5,3,13)
    hold on
    plot(Ca,ETR,".-",color=mycolor)
    box on
    ylabel("ETR [µmol e^{-1} s^{-1}]",FontSize=Fsize)
    xlabel("Ambient CO_{2} [µmol mol^{-1}]",FontSize=Fsize)

    s(14) = subplot(5,3,14)
    hold on
    plot(Ca,qL_1,".-",color=mycolor)
    box on
    ylabel("1-qL",FontSize=Fsize)
    xlabel("Ambient CO_{2} [µmol mol^{-1}]",FontSize=Fsize)
    Legend{i} = names(i);  
end

sgtitle("CO_{2} experiment")
% newPosition = [0.5 0.05 0.01 0.01];
% newUnits = 'normalized';
% hL = legend(Legend,"orientation","horizontal","location","southoutside")
% set(hL,'Position', newPosition,'Units', newUnits,"FontSize",12);
saveas(gcf,"./Figures/Vilfan_GasExc_co2.png")