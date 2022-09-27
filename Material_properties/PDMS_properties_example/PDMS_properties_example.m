% creates figures for estimating viscosity of different PDMS/corundum sand
% mixtures Modified after Zwaan, F., Schreurs, G., Ritter, M., and Rosenau,
% M. 2018, Rheology of PDMS-corundum sand mixtures from the Tectonic
% Modelling Lab (TecLab) of the University of Bern (CH)

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
data = load('Results.txt');
all_data = data(:,:);
clearvars data

% assign data
shear_rate   = log10(all_data(:,1));

tiledlayout(1,3);
set(gcf,'Units','normalized','Position',[.1 .1 .7 .6])

% CALCULATE STRESS EXPONENT
% ======================================================================= %
nexttile
[m,n]        = size(all_data);
sample_slope = zeros(1,(n-1)/2);
col_num      = 1;
range        = 10;

for isample = 2:2:n
    
    shear_stress = log10(all_data(:,isample));
    viscosity    = all_data(:,isample+1);
    
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
legend('PDMS','PDMS-crds rho=1.2','PDMS-crds rho=1.3','PDMS-crds rho=1.4',...
       'PDMS-crds rho=1.5','PDMS-crds rho=1.6','location','SouthEast')

hAx=gca;
hAx.LineWidth=2;
hAx.FontSize = 14;

axis square
axis([5e-5 2e-1 1e0 1e4])

box on
set(gca, 'Layer', 'Top')

disp(1/mean(sample_slope))

% VISCOSITY UNCORRECTED
% ======================================================================= %
nexttile
col_num = 1;
for isample = 2:2:n
    
    shear_stress = log10(all_data(:,isample));
    viscosity    = log10(all_data(:,isample+1));
    
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
legend('PDMS','PDMS-crds rho=1.2','PDMS-crds rho=1.3','PDMS-crds rho=1.4',...
       'PDMS-crds rho=1.5','PDMS-crds rho=1.6','location','SouthEast')

hAx=gca;
hAx.LineWidth=2;
hAx.FontSize = 14;

axis square
axis([5e-5 2e-1 1e4 2e5])

box on
set(gca, 'Layer', 'Top')

% VISCOSITY CORRECTED
% ======================================================================= %
nexttile
col_num = 1;
for isample = 2:2:n
    
    shear_stress = log10(all_data(:,isample));
    viscosity    = log10(all_data(:,isample+1));
    
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
legend('PDMS','PDMS-crds rho=1.2','PDMS-crds rho=1.3','PDMS-crds rho=1.4',...
       'PDMS-crds rho=1.5','PDMS-crds rho=1.6','location','SouthEast')

hAx=gca;
hAx.LineWidth=2;
hAx.FontSize = 14;

axis square
axis([5e-5 2e-1 1e4 2e5])

box on
set(gca, 'Layer', 'Top')
drawnow

% SAVING FIGURE
% ======================================================================= %
h = gcf;
set(h, 'PaperPositionMode','auto')
print('-dpng','-r600','-noui','fig_viscosity.png')
print(gcf,'-depsc','-r300','fig_viscosity.eps')
