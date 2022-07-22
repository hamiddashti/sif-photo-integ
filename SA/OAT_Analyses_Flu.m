%% The purpose of this 
%% Step 1: set paths
clear all 
clc
my_dir = '/home/hamid/SIF/sif-photo-integ/SA/safe_R1.1/';
% my_dir = './SA/safe_R1.1'; % use the 'pwd' command if you have already setup the Matlab
% current directory to the SAFE directory. Otherwise, you may define
% 'my_dir' manually by giving the path to the SAFE directory, e.g.:
% my_dir = '/Users/francescapianosi/Documents/safe_R1.0';
cd(my_dir)
datdir = './../../data/output/fluspect_output/2015/2019-06-19-1452/';
% Set current directory to 'my_dir' and add path to sub-folders:

addpath(genpath(my_dir));
out_dir = strcat('/home/hamid/SIF/sif-photo-integ/working/Figures/');

%% Step 2: setup the model and define input ranges

% load fluspect data
iwle = (400:750)-349;
load([datdir 'spectral.mat']);
load([datdir 'measured.mat']);
load([datdir 'leafopt.mat']);
load([datdir 'Tsimulated.mat']);
load([datdir 'Rsimulated.mat']);
load([datdir 'leafbio.mat']);
load([datdir 'optipar.mat']);

%%
% Randomly select a sample from Vilfan dataset
k=17;
biovars = {'Cab','Cca','Cdm','Cw','Cs','N','Cx','Cant','fqe1','fqe2'};
n=length(biovars);
n_rand = 50;
% generate n numbers for each fluspect parameter 
tmp.Cab = linspace(0,100,n_rand);
tmp.Cca = linspace(0,30,n_rand);
tmp.Cdm = linspace(0,0.5,n_rand);
tmp.Cw = linspace(0,0.4,n_rand);
tmp.Cs = linspace(0,0.6,n_rand);
tmp.N = linspace(1,4,n_rand);
tmp.Cx = linspace(0,1.5,n_rand);
tmp.Cant = linspace(0,10,n_rand);
tmp.fqe1 = linspace(0.0001,0.2,n_rand);
tmp.fqe2 = linspace(0.0001,0.2,n_rand);

% iwlE = find(spectral.wlM>399 & spectral.wlM<751);
E = (measured.E(:,k)); % Measured PAR 

%initialize to the middle value and change each in the loops later
leafbio.Cab = tmp.Cab(ceil(n_rand/2));
leafbio.Cca = tmp.Cca(ceil(n_rand/2));
leafbio.Cdm = tmp.Cdm(ceil(n_rand/2));
leafbio.Cw = tmp.Cw(ceil(n_rand/2));
leafbio.Cs = tmp.Cs(ceil(n_rand/2));
leafbio.N = tmp.N(ceil(n_rand/2));
leafbio.Cx = tmp.Cx(ceil(n_rand/2));
leafbio.Cant = tmp.Cant(ceil(n_rand/2));
leafbio.fqe(1) = tmp.fqe1(ceil(n_rand/2));
leafbio.fqe(2) = tmp.fqe2(ceil(n_rand/2));

map = jet(n_rand)
figure (1)
set(gcf,'units','inch','Position',[50 50 15 12],'color','w')
% tiledlayout(10,4,'TileSpacing','tight','Padding','compact')
spfig = zeros(10,4);
counter = 1;
for i = 1:n
    for j = 1:n_rand
        if i==9
            leafbio_loop = leafbio;
            leafbio_loop.fqe(1) = tmp.(biovars{i})(j);
            leafopt = main_flu(spectral,leafbio_loop,optipar,...
                "not_combined",E);
            hold on
            spfig3(1,i) = subplot(10,4,1+(i-1)*4)
            plot(spectral.wlP,leafopt.refl,'color',map(j,:))
            
            hold on
            spfig3(2,i) = subplot(10,4,2+(i-1)*4)
            plot(spectral.wlP,1-leafopt.tran,'color',map(j,:))
            
            hold on 
            spfig3(3,i) = subplot(10,4,3+(i-1)*4)
            plot(spectral.wlF,leafopt.Fu,'color',map(j,:))
            hold on
            spfig3(4,i) = subplot(10,4,4+(i-1)*4)
            plot(spectral.wlF,leafopt.Fd,'color',map(j,:))
        elseif i==10
            leafbio_loop = leafbio;
            leafbio_loop.fqe(2) = tmp.(biovars{i})(j);
            leafopt = main_flu(spectral,leafbio_loop,optipar,...
                "not_combined",E);
            
            hold on
            spfig3(1,i) = subplot(10,4,1+(i-1)*4)
            plot(spectral.wlP,leafopt.refl,'color',map(j,:))
            
            hold on
            spfig3(2,i) = subplot(10,4,2+(i-1)*4)
            plot(spectral.wlP,1-leafopt.tran,'color',map(j,:))
            
            hold on 
            spfig3(3,i) = subplot(10,4,3+(i-1)*4)
            plot(spectral.wlF,leafopt.Fu,'color',map(j,:))
            hold on
            spfig3(4,i) = subplot(10,4,4+(i-1)*4)
            plot(spectral.wlF,leafopt.Fd,'color',map(j,:))
        else
            
            leafbio_loop = leafbio;
            leafbio_loop.(biovars{i}) = tmp.(biovars{i})(j);
            leafopt = main_flu(spectral,leafbio_loop,optipar,...
                "not_combined",E);
            hold on
            spfig3(1,i) = subplot(10,4,1+(i-1)*4)
            plot(spectral.wlP,leafopt.refl,'color',map(j,:))
            
            hold on
            spfig3(2,i) = subplot(10,4,2+(i-1)*4)
            plot(spectral.wlP,1-leafopt.tran,'color',map(j,:))
            
            hold on 
            spfig3(3,i) = subplot(10,4,3+(i-1)*4)
            plot(spectral.wlF,leafopt.Fu,'color',map(j,:))
            hold on
            spfig3(4,i) = subplot(10,4,4+(i-1)*4)
            plot(spectral.wlF,leafopt.Fd,'color',map(j,:))
            counter = counter+1;
        end

    end

end
spfig3 = spfig3';
ylabel(spfig3(1,1),'Cab','FontSize',12,'FontWeight','bold')
ylabel(spfig3(2,1),'Cca','FontWeight','bold')
ylabel(spfig3(3,1),'Cdm','FontWeight','bold')
ylabel(spfig3(4,1),'Cw','FontWeight','bold')
ylabel(spfig3(5,1),'Cs','FontWeight','bold')
ylabel(spfig3(6,1),'N','FontWeight','bold')
ylabel(spfig3(7,1),'Cx','FontWeight','bold')
ylabel(spfig3(8,1),'Cant','FontWeight','bold')
ylabel(spfig3(9,1),'fqe1','FontWeight','bold')
ylabel(spfig3(10,1),'fqe2','FontWeight','bold')
for k = 1:4
    xlabel(spfig3(10,k),'wl (nm)')
end
title(spfig3(1,1),'Refl')
title(spfig3(1,2),'Tran')
title(spfig3(1,3),'Fu')
title(spfig3(1,4),'Fd') 
saveas(gcf,strcat(out_dir,"/OAT_Analyses.png"))

