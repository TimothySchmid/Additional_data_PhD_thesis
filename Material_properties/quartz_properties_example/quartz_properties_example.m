% creates figures for estimating mechanical properties of quartz sand and
% corundum sand Modified after Schmid, T., Schreurs, G., Warsitzka, M., and
% Rosenau, M., 2020, Effect of sieving height on density and friction of
% brittle analogue material: Ring-shear test data of quarz sand used for
% analogue experiments in the Tectonic Modelling Lab of the University of
% Bern

clear
close all
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
     
% parameters
A  = 0.022619; % area of shear zone (= surface area of lid) (m^2)
li = 0.0776;   % inner lever (center cell to center of shear zone)
lo = 0.1250;   % outer lever (center cell to hinge point)
v  =      3;   % shear velocity (mm/min)
N  =      2;   % Number of repetitions 


% FRICTION DATA
% ======================================================================= %
filename = '443-03_UB_quartzsand_30_ts.txt';
data     = readtable(filename);

[m,n] = size(data);

% get arrays
time  = data(2:end,1);

shear_displacement = table2array(time) * v/60;

col_num      = [1,1,2,2,3,3,4,4,5,5];

tiledlayout(1,3);
set(gcf,'Units','normalized','Position',[.3 .2 .7 .6])

% plotting
nexttile
for idata = 1:n-1
   
force = data(2:end,idata+1);

shear_stress = (table2array(force) * lo) / (li * A);

plot(shear_displacement,smoothdata(shear_stress,'gaussian',1),'-','LineWidth',2,'color',c_map(col_num(idata),:))
hold on
end

title('Shear displacement vs Shear stress')
xlabel('Shear displacement (mm)')
ylabel('Shear stress (Pa)')

hAx=gca;
hAx.LineWidth=2;
hAx.FontSize = 14;

axis square
axis([0 30 0 1600])

box on
set(gca, 'Layer', 'Top')

clearvars data A li lo v N


% REGRESSION
% ======================================================================= %
nexttile
files = {'443-03_UB_quartzsand_30_peak.txt',...
         '443-03_UB_quartzsand_30_dynamic.txt',...
         '443-03_UB_quartzsand_30_reactivation.txt'};
   
col_num      = 1;

for ifile = 1:length(files)
filename = files{ifile};
data     = readtable(filename);
data     = sortrows(data,'NormalStress_Pa_');

[m,n] = size(data);

normal_stress = table2array(data(:,1));
shear_stress  = table2array(data(:,2));

f = fit(normal_stress,shear_stress,ft);

P1 = m * sum(normal_stress .* shear_stress) - sum(normal_stress) * sum(shear_stress);
P2 = m * sum(normal_stress .^ 2) - sum(normal_stress)^2;

slope = P1/P2;
res   = shear_stress - f.q - f.m * normal_stress;
s     = sum(res.^2) / (m-2);

std_m = sqrt(m * s / P2);
std_q = sqrt(sum(normal_stress.^2) * s / P2);

plot(normal_stress,f.m*normal_stress + f.q,...
    'color',c_map(col_num,:),'LineWidth',2,'handlevisibility','off')
hold on
plot(normal_stress,shear_stress,'.','MarkerSize',20,'color',c_map(col_num,:))


col_num = col_num + 1;
end

title('Linear regression')
xlabel('Normal stress (Pa s)')
ylabel('Shear stress (Pa s)')
legend(strrep(files, '_', ' '),'location','SouthEast')

hAx=gca;
hAx.LineWidth=2;
hAx.FontSize = 14;

axis square
axis([0 2500 0 1600])

box on
set(gca, 'Layer', 'Top')


% MUTUAL TWO POINT REGRESSION
% ======================================================================= %
nexttile
for ifile = 1:1%length(files)
filename = files{ifile};
data     = readtable(filename);
data     = sortrows(data,'NormalStress_Pa_');

[m,n] = size(data);

num_reg = factorial(10)/(factorial(2) * factorial(10-2));
MU = zeros(num_reg,1);
C  = zeros(num_reg,1);
ipos = 1;

normal_stress = table2array(data(:,1));
shear_stress  = table2array(data(:,2));

for ipoint = 1:m
   
    x_fix = normal_stress(ipoint,1);
    y_fix = shear_stress(ipoint,1);
    
    count_array = 1:m;
    excl_array  = count_array(find(count_array~=ipoint));
    
    for iflex = excl_array
        x_flex = normal_stress(iflex,1);
        y_flex = shear_stress(iflex,1);
        
        
        dx = x_fix - x_flex;
        dy = y_fix - y_flex;

        
        MU(ipos,1) = dy/dx;
        C(ipos,1)  = y_fix - dy/dx * x_fix;
        
        ipos = ipos + 1;
    end
end

MU = MU(MU~=0 & isfinite(MU));
C  = C(C~=0 & isfinite(C));

end

pd = fitdist(MU,'normal');
x_gauss = 0.6:0.001:0.8;
f = exp(-(x_gauss-pd.mu).^2./(2*pd.sigma^2))./(pd.sigma*sqrt(2*pi));

gauss_curve = normpdf(MU,pd.mu,pd.sigma);

h = histogram(MU,9);
h.EdgeColor = 'w';
h.FaceColor = c_map(1,:);
h.FaceAlpha = 1;

hold on
plot(x_gauss,f,'LineWidth',2)

title('Mutual two-point regression')
xlabel('Friction coefficient µ')
ylabel('Counts')
% legend(strrep(files, '_', ' '),'location','NorthWest')

hAx=gca;
hAx.LineWidth=2;
hAx.FontSize = 14;

axis square
axis([0.65 0.8 0 30])

box on
set(gca, 'Layer', 'Top')


% SAVING FIGURE
% ======================================================================= %
h = gcf;
set(h, 'PaperPositionMode','auto')
print('-dpng','-r600','-noui','brittle_material.png')
print(gcf,'-depsc','-r300','brittle_material.eps')
