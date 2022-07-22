% Prepare Nastassia Vilfan's data
function vilfan_data = prepare_vilfan_data(data_dir)

% data_dir = '/home/hamid/SIF/sif-photo-integ/data/Nastassia_Vilfan/';

info = readtable(strcat(data_dir,'Data_ACiCurves.xlsx'),...
    'Sheet','Info','NumHeaderLines',1);
wl = 350:1:2500;
%% The Co2 experiment

% Licor/PAM
licor_pam = readtable(strcat(data_dir,'Data_ACiCurves.xlsx'),...
    'Sheet','Li6400+PAM');
Vars = licor_pam.Properties.VariableNames;
leaf_id = licor_pam.Leaf;

empty=  cell2table(cell(0,length(Vars)),'VariableNames',Vars);
for i = 1:20
    b = licor_pam(leaf_id==i,:);
    c = sortrows(b,'CO2');
    empty = [empty;c];
end

for i = 1:length(Vars)
    data_co2.licor_pam.(Vars{i}) = empty.(Vars{i});
end 

data_co2.leaf_id = leaf_id;

% Trans
data_co2.Ti = table2array(readtable(strcat(data_dir,'Data_ACiCurves.xlsx'),...
    'Sheet','Ti'));

data_co2.TnormRatio565 = table2array(readtable(strcat(data_dir,'Data_ACiCurves.xlsx'),...
    'Sheet','TnormRatio565'));
data_co2.If = table2array(readtable(strcat(data_dir,'Data_ACiCurves.xlsx'),...
    'Sheet','If')); 

data_co2.tau = table2array(readtable(strcat(data_dir,'Data_ACiCurves.xlsx'),...
    'Sheet','tau')); 

data_co2.Fd = table2array(readtable(strcat(data_dir,'Data_ACiCurves.xlsx'),...
    'Sheet','Fd')); 
%% The LRC experiment
% Licor/PAM
licor_pam = readtable(strcat(data_dir,'Data_LightCurves.xlsx'),...
    'Sheet','Li6400+PAM');
Vars = licor_pam.Properties.VariableNames;
Leaf_id = licor_pam.Leaf;

for i = 1:length(Vars)
    data_lrc.licor_pam.(Vars{i}) = licor_pam.(Vars{i})
end 
data_lrc.leaf_id = Leaf_id;
% Trans
data_lrc.Ti = table2array(readtable(strcat(data_dir,'Data_LightCurves.xlsx'),...
    'Sheet','Ti'));

data_lrc.TnormRatio565 = table2array(readtable(strcat(data_dir,'Data_LightCurves.xlsx'),...
    'Sheet','TnormRatio565'));
data_lrc.If = table2array(readtable(strcat(data_dir,'Data_LightCurves.xlsx'),...
    'Sheet','If')); 

data_lrc.tau = table2array(readtable(strcat(data_dir,'Data_LightCurves.xlsx'),...
    'Sheet','tau')); 

data_lrc.Fd = table2array(readtable(strcat(data_dir,'Data_LightCurves.xlsx'),...
    'Sheet','Fd')); 

vilfan_data.info = info;

vilfan_data.data_co2 = data_co2;
vilfan_data.data_lrc = data_lrc;

end
