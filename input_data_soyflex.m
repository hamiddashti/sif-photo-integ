function [wl,measured,include,target,c,wlrange,leafbio,outdirname] = input_data_soyflex(year,day) %#ok<STOUT>

% The data of the campaigns (KCA (2015) and Udine (2016)) are not in the
% same format. I decided to keep them in the format they were provided to 
% to me by the people who did the measurements (Nastiassia Vilfan in 2015 
% and MaPi Cendrero in 2016), and load them separately rather than first 
% changing the data format (CvdT)

switch year
    case 2015
        %2015 KCA data
        load('../data/measured/Fluowat/2015/CKA2015_Diurnal_ITCvFIN.mat');
        wl = G.wl{:,:};
        measured.tran = G.T{:,:};
        measured.refl = G.R{:,:};
        measured.Euf = G.If{:,:};
        measured.E = G.If{:,:};
        measured.Fu = G.Fu{:,:};
        measured.Fd = G.Fd{:,:};  
        measured.variety = G.IDs(:,3);
        measured.variety = measured.variety{:,:};
    case 2016
        if day==203
            % %2016 Udine data 22 July
            
            direc = '../data/measured/Fluowat/2016/';
            wl = dlmread([direc '2016_07_22_Reflectance_Top_Sun.txt'],'',1,0);
            wl = wl(:,1);
            measured.refl = dlmread([direc '2016_07_22_Reflectance_Top_Sun.txt'],'',1,1);
            measured.tran = dlmread([direc '2016_07_22_transmittance_Top_Sun.txt'],'',1,1);
            measured.E = 1E4*dlmread([direc '2016_07_22_Incoming_Radiance_with_Filter_Top_Sun.txt'],'',1,1);
            measured.Euf = 1E4*dlmread([direc '2016_07_22_Incoming_Radiance_Top_Sun.txt'],'',1,1);
            measured.Fu = 1E4*dlmread([direc '2016_07_22_Fup_Top_Sun.txt'],'',1,1);
            measured.Fd = 1E4*dlmread([direc '2016_07_22_Fdw_Top_Sun.txt'],'',1,1);
  
        elseif day==207
            % %2016 Udine data 26 July
            
            direc = '/../data/measured/Fluowat/2016/';
            wl = dlmread([direc '2016_07_26_Reflectance_Top_Sun.txt'],'',1,0);
            wl = wl(:,1);
            measured.refl = dlmread([direc '2016_07_26_Reflectance_Top_Sun.txt'],'',1,1);
            measured.tran = dlmread([direc '2016_07_26_transmittance_Top_Sun.txt'],'',1,1);
            measured.E = 1E4*dlmread([direc '2016_07_26_Incoming_Radiance_with_Filter_Top_Sun.txt'],'',1,1);
            measured.Euf = 1E4*dlmread([direc '2016_07_26_Incoming_Radiance_Top_Sun.txt'],'',1,1);
            measured.Fu = 1E4*dlmread([direc '2016_07_26_Fup_Top_Sun.txt'],'',1,1);
            measured.Fd = 1E4*dlmread([direc '2016_07_26_Fdw_Top_Sun.txt'],'',1,1);
  
            
        elseif day==210
            %2016 Udine data 29 July       
            direc = '/../data/measured/Fluowat/2016/';
            wl = dlmread([direc '2016_07_29_Sun_Reflectance_Top_Md_Bt.txt'],'',1,0);
            wl = wl(:,1);
            measured.refl = dlmread([direc '2016_07_29_Sun_Reflectance_Top_Md_Bt.txt'],'',1,1);
            measured.tran = dlmread([direc '2016_07_29_Sun_Transmittance_Top_Md_Bt.txt'],'',1,1);
            measured.Euf = 1E3*dlmread([direc '2016_07_29_Sun_Incoming_Radiance_Top_Md_Bt.txt'],'',1,1);
        end
end
measured.std = .03*ones(length(wl),size(measured.refl,2));

%%
measured.refl(measured.Euf<2) = NaN;
measured.tran(measured.Euf<2) = NaN;

%% Tell me, where should I write the output
outdirname              = ['../data/output/fluspect_output/' num2str(year)];
if year == 2016, outdirname  = [outdirname '/' num2str(day)]; end

%% Which columns in the reflectance/transmittance should I use ?
% in other words: how many spectra do you want to tune?
c               = -999; % 1: first column ; 2: second column, [1,2]: first and second, ect
% -999: all columns in the file

include.Cab     = 1;
include.Cdm     = 1;
include.Cw      = 1;
include.Cs      = 1;
include.Cca     = 1;
include.Cx      = 1;
include.Cant    = 1;
include.N       = 1;

%% which outputs should I calibrate?
target          = 0; %#ok<*NASGU> %0: calibrate T&R, 1: calibrate only R; 2: calibrate only T

%% which spectral region should I calibrate?
wlmin           = 400;          % starting wavelength (nm)
wlmax           = 2400;         % ending wavelength (nm)

%% initialize parameters for retrieval
% these will be calibrated to your reflectance data if you said so above
leafbio.Cab     = 30;           % chlorophyll content               [ug cm-2]
leafbio.Cdm     = 0.001;        % dry matter content                [g cm-2]
leafbio.Cw      = 0.002;        % leaf water thickness equivalent   [cm]
leafbio.Cs      = 0.15;          % senescent material                [fraction]
leafbio.Cca     = 4.3;            % carotenoids                       [?]
leafbio.Cant    = 0;
leafbio.N       = 1.5;          % leaf structure parameter (affects the ratio of refl: transmittance) []
leafbio.Cx      = 0;
leafbio.fqe     = 0.01;                     % quantum yield
leafbio.V2Z     = 0;

wlrange.wlmin   = wlmin;
wlrange.wlmax   = wlmax;

%load Optipar2017_ProspectD.mat
%optipar.phi = optipar.phiII; %#ok<NODEF>
%load optiparnew
load ../data/optical_coefficients/optipar2019
