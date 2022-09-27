function fct_plot_xslice(INPUT, CELL)

% fct_plot_zslice
%
% PLOTTING SLICES NORMAL TO X AXIS
%
% Input:
%   - INPUT      --> Structure containing initial input parameters
%   - CELL       --> Structure containing  displacements and coordinates
%
% ======================================================================= %

pos  = INPUT.slice_position;
incr = INPUT.vector_incr;

l = size(CELL.disp.X,2);

max_pos = size(CELL.disp.X{1},1);
if pos > max_pos
    error(['Slice position exceeds maximum position of ', num2str(max_pos), '!'])
end

ycoords = double(CELL.coords.ycoords);
zcoords = double(CELL.coords.zcoords);

Y = double(CELL.coords.Y);
Z = double(CELL.coords.Z);

% Plotting 2d in the xz plane --> perpendicular to rift-axis
for itime = 1:1:l
    
    % get coordinates 2d field for pcolor field
    plt_hcoord_2d = squeeze(Z(pos,:,:));
    plt_vcoord_2d = squeeze(Y(pos,:,:));
    
    % get coordinates for vectors
    plt_hcoord_vec_2d = squeeze(Z(pos,1:incr:end,1:incr:end));
    plt_vcoord_vec_2d = squeeze(Y(pos,1:incr:end,1:incr:end));
    
    % get x vector 2d field
    plt_vec_h2d = double(squeeze(CELL.disp.Z{itime}(pos,:,:)));
    plt_vec_h2d = plt_vec_h2d(1:incr:end,1:incr:end);
    
    % get y vector 2d field
    plt_vec_v2d = double(squeeze(CELL.disp.Y{itime}(pos,:,:)));
    plt_vec_v2d = plt_vec_v2d(1:incr:end,1:incr:end);
    
    switch INPUT.displacement
        case 'x_displacement'
            plot_what_2d = double(squeeze(CELL.disp.X{itime}(pos,:,:)));
            plot_what_2d(plot_what_2d==0) = NaN;
            crange = [min(min(plot_what_2d)), max(max(plot_what_2d))];
            collim = round([-max(abs(crange)) max(abs(crange))],0);
            
        case 'y_displacement'
            plot_what_2d = double(squeeze(CELL.disp.Z{itime}(pos,:,:)));
            plot_what_2d(plot_what_2d==0) = NaN;
            crange = [min(min(plot_what_2d)), max(max(plot_what_2d))];
            collim = round([-max(abs(crange)) max(abs(crange))],0);
            
        case 'z_displacement'
            plot_what_2d = double(squeeze(CELL.disp.Y{itime}(pos,:,:)));
            plot_what_2d(plot_what_2d==0) = NaN;
            crange = [min(min(plot_what_2d)), max(max(plot_what_2d))];
            collim = round([-max(abs(crange)) max(abs(crange))],0);
                
        case '2D_tot_displacement'
            x_now = double(squeeze(CELL.disp.X{itime}(pos,:,:)));
            y_now = double(squeeze(CELL.disp.Y{itime}(pos,:,:)));
            plot_what_2d = sqrt(x_now.^2 + y_now.^2);
            plot_what_2d(plot_what_2d==0) = NaN;
            maxc = max(max(plot_what_2d));
            collim = [0, round(maxc,0)];
            
        case '3D_tot_displacement'
            plot_what_2d = double(squeeze(CELL.disp.TOT3D{itime}(pos,:,:)));
            plot_what_2d(plot_what_2d==0) = NaN;
            maxc = max(max(plot_what_2d));
            collim = [0, round(maxc,0)];
            
        otherwise
            error('displacement not defined')
    end
    
    figure(1)
    clf
    set(gcf,'Units','normalized','Position',[.15 .2 .7 .5])
    c_map = fct_colormap(INPUT);
    colormap(flipud(c_map))
    
    pcolor(plt_hcoord_2d,plt_vcoord_2d,plot_what_2d)
    shading interp
    hold on
    
    % Contour lines -------------------------------------------------------
    switch INPUT.contour_lines
        case 'yes'
            v_cont = [min(collim):2:max(collim) 0];
            [C,h]  = contour(plt_hcoord_2d,plt_vcoord_2d,plot_what_2d,...
                v_cont,'k','LineWidth',2);
            
            v_label = v_cont;
            clabel(C,h,v_label,'color',INPUT.contour_color,'FontSize',10)
            h.LineColor = INPUT.contour_color; h.LineWidth = INPUT.contour_width;
            
        case 'no'
        otherwise
            error('contour lines undefined')
    end
    % ---------------------------------------------------------------------
    
    % Vector arrows -------------------------------------------------------
    switch INPUT.vector_arrows
        case 'yes'
            quiver(plt_hcoord_vec_2d,plt_vcoord_vec_2d,plt_vec_h2d,...
                plt_vec_v2d,'color',INPUT.vector_color,'LineWidth',INPUT.vector_width,'AutoScale','off')
        case 'no'
        otherwise
            error('vector arrows undefined')
    end
    % ---------------------------------------------------------------------
    
    cb = colorbar;
    cb.Label.String = INPUT.displacement;
    cb.Label.Interpreter = 'none';
    caxis(collim)
    
    title(['slide position ' num2str(pos), ' time: ' num2str(itime * CELL.exp.dt)])
    
    hAx=gca;
    hAx.LineWidth=2;
    hAx.FontSize = 14;
    axis equal
    axis([zcoords(1)-10 zcoords(end)+10 ycoords(1)-2 ycoords(end)+2])
    
    box on
    set(gca, 'Layer', 'Top')
    drawnow
    
    hold off
    
    % SAVE FIGURE
    % =================================================================== %
    fct_save_plots(INPUT, CELL, itime)
end
end

