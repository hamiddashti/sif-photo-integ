clear
clc


spectral = define_bands;

iwle = (400:750)-349;
iwlF = (640:850)-399;
datdir = './data/output/fluspect_output/2015/2019-06-19-1452/';
load([datdir 'spectral.mat']);
load([datdir 'measured.mat']);
load([datdir 'leafopt.mat']);
load([datdir 'Tsimulated.mat']);
load([datdir 'Rsimulated.mat']);
load([datdir 'leafbio.mat']);
load([datdir 'optipar.mat']);

%% Figure A.11, old and new Optipar fluorescence spectra, plotted for leaf 28
k = 28;
leafbio2.Cab= leafbio_all.Cab(k);
leafbio2.Cdm= leafbio_all.Cdm(k);
leafbio2.Cw= leafbio_all.Cw(k);
leafbio2.Cs= leafbio_all.Cs(k);
leafbio2.Cca= leafbio_all.Cca(k);
leafbio2.N= leafbio_all.N(k);
leafbio2.Cx= leafbio_all.Cx(k);
leafbio2.Cant= leafbio_all.Cant(k);
leafbio2.fqe(1) = leafbio_all.fqe_opt(k)*0.2;
leafbio2.fqe(2) = leafbio_all.fqe_opt(k)*0.8;

leafbio3 = leafbio2;
leafbio3.fqe = leafbio_all.fqe_opt(k);
[leafopt_old] = fluspect_B_CX(spectral,leafbio2,optipar);
[leafopt_new] = fluspect_B_CX_PSI_PSII_combined(spectral,leafbio3,optipar);

leafopt_old.Fu = (leafopt_old.MbI+leafopt_old.MbII)* interp1(spectral.wlM,measured.E(:,k),spectral.wlE)';
leafopt_old.Fd = (leafopt_old.MfI+leafopt_old.MfII)* interp1(spectral.wlM,measured.E(:,k),spectral.wlE)';
leafopt_new.Fu = (leafopt_new.Mb)* interp1(spectral.wlM,measured.E(:,k),spectral.wlE)';
leafopt_new.Fd = (leafopt_new.Mf)* interp1(spectral.wlM,measured.E(:,k),spectral.wlE)';
iwlE = find(spectral.wlM>399 & spectral.wlM<750);
Etot = 1E-3*Sint(measured.E(iwle,k),spectral.wlM(iwle));

Fum = interp1(spectral.wlM,measured.Fu(:,k),spectral.wlF)
Fum(spectral.wlF<660) = NaN;
Fdm = interp1(spectral.wlM,measured.Fd(:,k),spectral.wlF);
Fdm(spectral.wlF<660) = NaN;

chi2 = chi_fun(spectral, leafbio2,optipar,"not_combined",measured.E(:,k),Fum,Fdm)


figure(11), clf
s(1) = subplot(131);
z=plot(spectral.wlP,[ optipar.phiI optipar.phiII optipar.phi],'k');
set(z(1),'color',[.3 .3 .3])
set(z(3),'linewidth',2)
legend('PSI','PSII','new')
set(gca,'xlim',[640 850],'ylim',[0 0.03])
xlabel('wl (nm)')
ylabel('\phi')
text(710,0.027,'a','FontSize',12)
set(gca,'FontSize',12)

s(2) = subplot(132);
z=plot(spectral.wlF,[ Fum'/Etot  leafopt_old.Fu/Etot leafopt_new.Fu/Etot],'k');
set(z(1),'lineStyle','none','Marker','x','MarkerSize',3)
set(z(2),'color',[.3 .3 .3])
set(z(3),'linewidth',2)
set(gca,'xlim',[640 850],'ylim',[0 0.03])
ylabel('F_b/E_{tot} (\mum^{-1})')
xlabel('wl (nm)')
legend('measured','old','new')
text(710,0.027,'b','FontSize',12)
set(gca,'FontSize',12)

s(3) = subplot(133);
z=plot(spectral.wlF,[Fdm'/Etot leafopt_old.Fd/Etot leafopt_new.Fd/Etot  ],'k');
set(z(1),'lineStyle','none','Marker','x','MarkerSize',3)
set(z(2),'color',[.3 .3 .3])
set(z(3),'linewidth',2)
set(gca,'xlim',[640 850],'ylim',[0 0.03])
ylabel('F_f/E_{tot} (\mum^{-1})')
legend('measured','old','new')
xlabel('wl (nm)')
text(710,0.027,'c','FontSize',12)
set(gca,'FontSize',12)

resizefigure(s,3,1,.09,.22,.12,.07)
