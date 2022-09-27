% script for scaling analogue material properties and natural rock
% materials

clear
close all
clc

yes_diary = 1;

% PHYSICAL PARAMETERS

    g  = 9.81;      % m/s^2
    th = 4;         % runtime in h
    
% BRITTLE PARAMETER

    % model
    rho_b_m = 1560;        % kg/m^3
    C_m     = 50;          % Pa;
    mu_m    = 0.72;        % unitless
    ang_m   = atand(mu_m); % friction angle in degree

    % nature
    rho_b_n = 2700;        % kg/m^3
    C_n     = 50e6;        % Pa, Byerlee, 1978
    mu_n    = 0.6;         % unitless, Byerlee, 1978
    ang_n   = atand(mu_n); % friction angle in degree

% VISCOUS PARAMETER

    % model
    rho_v_m = 1600;        % kg/m^3
    eta_m   = 1e5;         % Pa s
    n_m     = 1.05;        % stress exponent

    % nature
    rho_v_n = 2900;        % kg/m^3
    eta_n   = 5e20;        % Pa s
    
% GEOMETRY

    % model
    h_b_m = 4.0e-2;        % m
    h_v_m = 6e-2-h_b_m;    % m
    v_m   = 10e-3;         % m/h
    v_m_s = v_m/60/60;     % m/s
    
    % nature
    h_b_n   = 15.0e3;      % m
    h_v_n   = 30e3-h_b_n;  % m
     
% SCALING

  % stress Pa
  rho_r    = rho_b_m/rho_b_n;
  g_r      = 1;
  h_r      = h_b_m/h_b_n;
  
  stress_m = rho_b_m * g * h_b_m;
  stress_r = rho_r * g_r * h_r;
  stress_n = stress_m/stress_r;
 
  % strain rate 1/s
  eta_r = eta_m/eta_n;
  
  strain_rate_m = v_m_s/(h_b_m+h_v_m);
  strain_rate_r = stress_r/eta_r;
  strain_rate_n = strain_rate_m/strain_rate_r;
  
  % time s
  t_m = th*60*60;
  t_r = 1/strain_rate_r;
  t_n = t_m/t_r;
  
  % velocity m/s
  
  v_r   = strain_rate_r*h_r;
  v_n   = v_m_s/v_r;
  
  % DYNAMIC SIMILARITY
  %  Smoluchowski number for brittle similarity
  sm = (rho_b_m*g*h_b_m)/(C_m+mu_m*rho_b_m*g*h_b_m);
  sn = (rho_b_n*g*h_b_n)/(C_n+mu_n*rho_b_n*g*h_b_n);
  
  % Ramberg number for viscous similarity
  rm = (rho_v_m*g*h_v_m^2)/(eta_m*v_m_s);
  rn = (rho_v_n*g*h_v_n^2)/(eta_n*v_n);
  
  % Reynold number for ratios between inertial forces and viscous forces
  rem = (rho_v_m*v_m_s*h_v_m)/eta_m;
  ren = (rho_v_n*v_n*h_v_n)/eta_n;
  
  % STRENGTH
  
%   % After Zwaan et al. 2019
%   % --------------------------------------------------------------------- %
%   strength_bm = rho_b_m * g * h_b_m^2 * sind(ang_m) + 2*C_m * h_b_m * cosd(ang_m);
%   strength_vm = 4 * eta_m * 1/2*v_m_s * h_v_m / 0.15;

%   % After Brun 2002
%   % --------------------------------------------------------------------- %
%   profile_bm = 2/3 * rho_b_m * g * h_b_m;
%   profile_vm = 1 * eta_m * strain_rate_m;
%   
%   % integrate over layer thickness
%   strength_bm = 1/3 * rho_b_m * g * h_b_m^2;
%   strength_vm = eta_m * strain_rate_m * h_v_m;
  
  % After Cruden
  % --------------------------------------------------------------------- %
  profile_bm = C_m + mu_m * rho_b_m * g * h_b_m;
  profile_vm = (eta_m * strain_rate_m)^(1/n_m);
  
  % integrate over layer thickness
  strength_bm = (C_m * h_b_m) + (1/2 * mu_m * rho_b_m * g * h_b_m^2);
  strength_vm = (eta_m * strain_rate_m)^(1/n_m) * h_v_m;
  
  strength_ratio = round(strength_bm/strength_vm);
  
  
  zm = [0 -h_b_m -h_b_m -(h_b_m+h_v_m)];
  xm = [C_m profile_bm profile_vm profile_vm];
  
  % PLOT
  figure(1)
  clf
  set(gcf,'Units','Normalized','Position',[.2 .2 .2 .7],'PaperPositionMode','auto')

  plot(xm,zm,'-','LineWidth',3)
  hold on
  
  xlabel('\sigma_{diff} [Pa]','FontSize',14)
  ylabel('Depth','FontSize',14)
  yticks(-0.06:0.01:0)
  yticklabels({'6 cm','5 cm','4 cm','3 cm','2 cm','1 cm','0 cm'})
  set(gca,'xaxisLocation','top')
  hAx=gca;
  hAx.LineWidth=2.5;
  hAx.FontSize = 14;
  
  print('-depsc','-r300','-noui',['Brittle thickness_',num2str(1e2*h_b_m),'cm.eps'])
%   print('-depsc','-r300','-noui','3_1_reference.eps')
  
  % PRINT
  if yes_diary
      diary(['Brittle thickness_',num2str(1e2*h_b_m),'cm.log'])
  end
  
  fprintf(['VALUES FOR SAND THICKNESS OF ', num2str(1e2*h_b_m),'cm'])
  fprintf('\n\n')
  fprintf('SCALING RATIOS\n')
  fprintf('================================ \n \n')
  fprintf('stress scaling :%1.E\n\n',stress_r)
  fprintf('strain rate scaling :%1.E\n\n',strain_rate_r)
  fprintf('time scaling :%1.E\n\n',t_r)
  fprintf('velocity scaling :%1.E\n\n\n',v_r)
  
  fprintf('MODEL VALUES ~ VALUES IN NATURE\n')
  fprintf('================================ \n \n')
  fprintf('1 hour ~ %.2f Ma\n\n',t_n/60/60/24/365/1e6/th)
  fprintf('1 cm ~ %.2f km\n\n',1e-2/h_r*1e-3)
  fprintf([num2str(v_m*1e3), 'mm/h ~ %1.f mm/a\n\n'],1e3*v_n*60*60*24*365)
  fprintf([num2str(strain_rate_m,'%1.E'), ' 1/s ~ ',num2str(strain_rate_n,'%1.E'),' 1/s'])
  fprintf('\n\n\n')
  
  fprintf('DYNAMIC SIMILARITIES\n')
  fprintf('================================ \n \n')
  fprintf('Smoluchowski number S_m\n\n')
  fprintf(['Model: ', num2str(sm,'%1.f'),'      Nature: ', num2str(sn,'%1.f')])
  fprintf('\n\n')
  fprintf('Ramberg number R_m\n\n')
  fprintf(['Model: ', num2str(rm,'%1.f'),'      Nature: ', num2str(rn,'%1.f')])
  fprintf('\n\n')
  fprintf('Reynold number R_e\n\n')
  fprintf(['Model: ', num2str(rem,'%1.E'),'      Nature: ', num2str(ren,'%1.e')])
  fprintf('\n\n\n')
  
  fprintf('STRENGTH RATIO MODELS\n')
  fprintf('================================ \n \n')
  fprintf(['brittle strength: ', num2str(strength_bm, '%1.1f'), '\n\n'])
  fprintf(['viscous strength: ', num2str(strength_vm, '%1.1f'), '\n\n'])
  fprintf(['ratio: ', num2str(strength_ratio, '%1.f'), '\n\n'])
  
  if yes_diary
      diary off
  end
  
