function aus_data = prepare_aus_data(licor_dir,asd_dir,qep_dir)
%
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
spectral = define_bands();

f_wl_low = 660;
f_wl_high = 800;

%% The ASD data
% asd_dir = '/home/hamid/SIF/sif-photo-integ/data/william_woodgate/ASD/'

data_asd = readtable(strcat(asd_dir,'tumba_spectra_LRC_CO2_comb_smooth.csv'));

for i = 1:9
    Rname = strcat(names(i),'_R');
    asd_R(:,i) = data_asd.(Rname);

    Tname = strcat(names(i),'_T');
    asd_T(:,i) = data_asd.(Tname);

    Aname = strcat(names(i),'_A');
    asd_A(:,i) = data_asd.(Aname);

end
wl_asd = spectral.wl_ASD';

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

for i = 1:length(names)
    if i==9
        % There is no CO2 data for this leaf
        continue
    end
    Qin_wave_co2(:,i) = par2wave(spectral,mean(licor_co2.Qin(:,i)));
end
licor_co2.Qin_wave = Qin_wave_co2;


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

% Dostributing the par values recorded by licor over PAR wl
Qin_wave = zeros(length(spectral.wlPAR),size(licor_lrc.Qin,1),size(licor_lrc.Qin,2));
for i = 1:length(names)
    Qin_wave(:,:,i) = par2wave(spectral,licor_lrc.Qin(:,i));
end
licor_lrc.Qin_wave = Qin_wave;


%%
% Gather all the data
aus_data.asd.refl = asd_R;
aus_data.asd.trans = asd_T;
aus_data.asd.abs = asd_A;
aus_data.asd.wl_asd = wl_asd;

aus_data.spec_co2 = spec_co2;
aus_data.licor_co2 = licor_co2;

aus_data.spec_lrc = spec_lrc;
aus_data.licor_lrc = licor_lrc;

aus_data.names = names;
aus_data.vars_licor = vars_licor;
aus_data.vars_qep = vars_qep;
end





