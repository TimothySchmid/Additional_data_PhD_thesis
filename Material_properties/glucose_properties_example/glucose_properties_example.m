% creates figures for estimating viscosity of different glucose syrup types
% Modified after Schmid, T., Zwaan, F., Corbi, F., Funiciello, F., and
% Schreurs, G., 2022, Rheology of glucose syrup from the Tectonic Modelling
% Lab (TecLab) of the University of Bern (CH)

close all
clear
clc

% colors
c_map = [0      0.4470 0.7410;...
         0.8500 0.3250 0.0980;...
         0.9290 0.6940 0.1250;...
         0.4940 0.1840 0.5560;...
         0.4660 0.6740 0.1880;...
         0.3010 0.7450 0.9330;...
         0.6350 0.0780 0.1840];

% linear fit
ft = fittype('m * x + q','dependent',{'y'},'independent',{'x'},...
             'coefficients',{'m','q'});

% import data
data_set = ["1_GlucoWheat_45","2_GlucoWheat_60",...
    "3_GlucoSweet_44","4_GlucoSweet_62"];
            
            
for iset = 1:length(data_set)
    set_name   = data_set(iset);
    db         = importdata(set_name + '.txt');
    dbv        = importdata(set_name + '_t.txt');
    DATA(iset) = struct('set_number',data_set{iset}(1),...
                        'name',data_set{iset}(isletter(set_name)),...
                        'DE',data_set{iset}(end-1:end),...
                        'shear_rate',db.data(:,2),...
                        'shear_stress',db.data(:,3),...
                        'shear_viscosity',db.data(:,4)./1e3,...
                        'temperature',dbv.data(:,5),...
                        'temp_viscosity',dbv.data(:,4)./1e3,...
                        'exponent',diff(log(db.data(:,2)./db.data(:,4)))./diff(log(db.data(:,3))));
                    
    clearvars db*
end

tiledlayout(1,3);
set(gcf,'Units','normalized','Position',[.1 .1 .7 .6])

% CALCULATE STRESS EXPONENT
% ======================================================================= %
nexttile
sample_slope = zeros(1,length(data_set));
col_num      = 1;
range        = 5;

for isample = 1:length(data_set)
    
   shear_rate   = log10(DATA(isample).shear_rate);
   shear_stress = log10(DATA(isample).shear_stress);
   
   f = fit(shear_rate(1:range),shear_stress(1:range),ft);
    
   loglog(10.^shear_rate(1:range),10.^(f.m*shear_rate(1:range) + f.q),...
        'color',c_map(col_num,:),'LineWidth',2,'handlevisibility','off')
   hold on
   loglog(10.^shear_rate,10.^shear_stress,'.','MarkerSize',20,'color',c_map(col_num,:))

   sample_slope(1,col_num) = f.m;
   col_num = col_num + 1;
end

title('Shear rate vs Shear stress')
xlabel('Shear rate (s^{-1})')
ylabel('Shear stress (Pa)')
legend(strrep(data_set, '_', ' '),'location','SouthEast')

hAx=gca;
hAx.LineWidth=2;
hAx.FontSize = 14;

axis square
axis([5e-4 2e1 1e-3 1e3])

box on
set(gca, 'Layer', 'Top')

disp(1/mean(sample_slope))

% VISCOSITY UNCORRECTED
% ======================================================================= %
nexttile
col_num = 1;
for isample = 1:length(data_set)
    
   shear_rate   = log10(DATA(isample).shear_rate);
   shear_stress = log10(DATA(isample).shear_stress);
   viscosity    = log10(DATA(isample).shear_viscosity);
    
    loglog(10.^shear_rate,10.^viscosity,...
        'color',c_map(col_num,:),'LineWidth',2,'handlevisibility','off')
    hold on
    loglog(10.^shear_rate,10.^viscosity,...
        '.','MarkerSize',20,'color',c_map(col_num,:))
    
    col_num = col_num + 1;
end

title('Viscosity from n = 1')
xlabel('Shear rate (s^{-1})')
ylabel('Viscosity (Pa s)')
legend(strrep(data_set, '_', ' '),'location','SouthEast')

hAx=gca;
hAx.LineWidth=2;
hAx.FontSize = 14;

axis square
axis([5e-4 2e1 1e0 5e2])

box on
set(gca, 'Layer', 'Top')


% VISCOSITY CORRECTED
% ======================================================================= %
nexttile
col_num = 1;
for isample = 1:length(data_set)
    
   shear_rate   = log10(DATA(isample).shear_rate);
   shear_stress = log10(DATA(isample).shear_stress);
   viscosity    = log10(DATA(isample).shear_viscosity);
    
    loglog(10.^shear_rate,10.^shear_stress ./ (10.^shear_rate).^mean(sample_slope),...
        'color',c_map(col_num,:),'LineWidth',2,'handlevisibility','off')
    hold on
    loglog(10.^shear_rate,10.^shear_stress ./ (10.^shear_rate).^mean(sample_slope),...
        '.','MarkerSize',20,'color',c_map(col_num,:))
    
    col_num = col_num + 1;
end

title(['Viscosity from n = ', num2str(1/mean(sample_slope),'%2.2f')])
xlabel('Shear rate (s^{-1})')
ylabel('Viscosity (Pa s)')
legend(strrep(data_set, '_', ' '),'location','SouthEast')

hAx=gca;
hAx.LineWidth=2;
hAx.FontSize = 14;

axis square
axis([5e-4 2e1 1e0 5e2])

box on
set(gca, 'Layer', 'Top')
drawnow

% SAVING FIGURE
% ======================================================================= %
% h = gcf;
% set(h, 'PaperPositionMode','auto')
% print('-dpng','-r600','-noui','fig_viscosity_syrup.png')
% print(gcf,'-depsc','-r300','fig_viscosity_syrup.eps')
