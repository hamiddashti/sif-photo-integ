clear all
%% options
retrieve_Prospect = 1;
calibrate_Phi = 1;
% data_dir = '/home/hamid/SIF/data/original/data' 

%% optical coefficient file for Fluspect
% load(strcat(data_dir, '/optical_coefficients/optipar2019.mat'));
load('../data/optical_coefficients/optipar2019.mat');
%% this part can be commented out once it the optical parameters have been retrieved and the output has been saved
datdir = '../data/output/fluspect_output/2015/2019-06-19-1452/'; %this is the directory where the retrieved data are saved
load([datdir 'spectral.mat']);
load([datdir 'measured.mat']);
load([datdir 'leafopt.mat']);
load([datdir 'Tsimulated.mat']);
load([datdir 'Rsimulated.mat']);
load([datdir 'leafbio.mat']);
load([datdir 'optipar.mat']);
% and these are the files in this directory.

% this is the spectral range from 640 to 850 nm where SIF takes place.
spectral.iwlF = (640:850)-399;

%% retrieve PROSPECT-D parameters
% this part of the code retrieves optical properties from the measured
% transmittance and reflectance spectra of 66 leaves.

% if retrieve_Prospect
%     [wl,measured,include,target,c,wlrange,leafbio,outdirname] = input_data_soyflex(2015); %#ok<*UNRCH>
%     if c==-999; c = (1:size(measured.refl,2));
%     else
%         if min(c)<0, fprintf('%s \r', 'negative incides are not allowed (except -999 for all spectra), program ends'), return, end
%         if max(c)>size(measured.refl,2), fprintf('%s \r', ['you only have ' num2str(size(measured.refl,2)) ' spectra, not ' num2str(max(c)) ' in your measurement file, program ends']), return, end
%     end
%     [er,~,~,spectral,reflFluspect,tranFluspect,leafopt,leafbio,leafopt_all, leafbio_all]= Fluspect_retrievals(leafbio,c,include,target,measured,wl,wlrange,optipar);
%     string          = clock;
%     Output_dir      =fullfile(outdirname, sprintf('%4.0f-%02.0f-%02.0f-%02.0f%02.0f/',[string(1) string(2) string(3) string(4) string(5)]));
%     %Output_dir      = Output_dir{1};
%     mkdir(Output_dir)
%     
%     save([Output_dir 'optipar'],'optipar');
%     save([Output_dir 'leafbio'],'leafbio_all');
%     save([Output_dir 'Rsimulated'],'reflFluspect')
%     save([Output_dir 'Tsimulated'],'tranFluspect')
%     save([Output_dir 'leafopt'],'leafopt_all')
%     save([Output_dir 'spectral'],'spectral')
%     save([Output_dir 'measured'],'measured')
%     
% end
% disp('ALL Done!')
%% calibrate Phi
% this is the code that has been used to calibrate the specific
% fluorescence emission spectrum of Chlorophyll in the leaf. It has been
% calibrated to a subset of the leaves with Chlorophyll between 20 and 50
% \mug cm^{-2}.

if calibrate_Phi
    leafbio.fqe = 0.01;
    %I = find(leafbio_all.Cab>20 & leafbio_all.Cab<50);
    J = [  2     7    10    18    22    24    26    29];
    % A selection of leaves of different Cab. The calibration is too time
    % consuming of all leaves are used!
    
    sample = J;
    leafbio_all2.Cab = leafbio_all.Cab(sample);
    leafbio_all2.Cdm = leafbio_all.Cdm(sample);
    leafbio_all2.Cs = leafbio_all.Cs(sample);
    leafbio_all2.Cca = leafbio_all.Cca(sample);
    leafbio_all2.N = leafbio_all.N(sample);
    leafbio_all2.Cx = leafbio_all.Cx(sample);
    leafbio_all2.Cant = leafbio_all.Cant(sample);
    leafbio_all2.Cw = leafbio_all.Cw(sample);
    
    [optipar]= Fluspect_retrievals_phi(leafbio_all2,target,measured,wl,optipar);
    save('../data/optical_coefficients/optipar_test.mat', 'optipar')
    disp("All Done!")
end