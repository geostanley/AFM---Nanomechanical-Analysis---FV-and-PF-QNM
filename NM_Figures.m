%%%=== NM_Figures ===%%%
 
% This is script 5/5 for the nanomechanical analysis procedure.

% This script generates publication quality figures.

%%  Enter load and save directories and file names etc

clear variables
close all
clc

% Make nice colours for plotting
N = 4;
C = linspecer(N);
c1 = C(1,:);
c2 = C(2,:);
c3 = C(3,:);
c4 = C(4,:);

display('NM_Figures')

%%%%%%=== File to be loaded === %%%%%%%%%%%%%%%%%

GenericFileName = '2kHz_test_cyto';
GenericSaveName = '2kHz_test_cyto';
shadecolour = c1;


%%%%%%%=== Data structure to be loaded
LoadFolder = ['Z:\Users\George\Documents\PhD\Data\'...
    'Nanomechanical_Outputs_Hertz_YM_CP\Test'];

%%%%%%%=== Output folder
OutputFolder = ['Z:\Users\George\Documents\PhD\Data\'...
    'Nanomechanical_Outputs_Hertz_YM_CP\Test'];
%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%=== Colour limits for stiffness heatmap (pN/nm)
cLow_hmap  = 0;
cHigh_hmap = 35;
%%%%%%%%%%%%%%%%%%%%%%%%

% per pore save name (pp) and averaged result save name
pp_figsavename  = ' - IndividualPoresSubPlot_';
all_figsavename = ' - AllHeatMapYMSubPlot_';

% Plots for individual pores or just averaged?
individual_subplots_on = 0;
individual_plots_on    = 0;

%% Enter plotting limits

% if want to use custom_cmap, use colormapeditor, decide on the map, and
% save it out as thiscolourmap.mat (you can save it out by using the code
% at the end of the script).
custom_cmap = 0;

% select shade of grey for plotting
grey = 0.65;

fc_x_min = -40;
fc_x_max =  80;
fc_y_min = -50;
fc_y_max = 400;
sc_y_min = -15;
sc_y_max =  40;

x_radial_min = -90;
x_radial_max =  90;
height_min   = -50;
height_max   =  50;

% E_eff limits for averaged pores
ym_min = 0;
ym_max = 5;

% E_eff limits for individual pores
ym_pp_min = 0;
ym_pp_max = 5;

% indiviudal pores, height image limits
height_mat_clim_min = -10;
height_mat_clim_max = 60;
%% Load the data structure

LoadFileName = strcat(GenericFileName, ' - NanomechanicalProcessedData', '.mat');
FullFileName = fullfile(LoadFolder, LoadFileName);

display('Loading data...')
load(FullFileName)

% pull-out the data for overall averaged pores
stiffhmap_pN_interp = RB_Processed.hmapplot.stiffhmap_pN_interp;
stiff_pN_hmap       = RB_Processed.hmapplot.stiff_pN_hmap;
th_plot_0nm         = RB_Processed.hmapplot.th_plot_0nm;
hd_plot_0nm         = RB_Processed.hmapplot.hd_plot_0nm;
rv_plot             = RB_Processed.rv_plot;
ym_plot             = RB_Processed.ym_plot;
hd_plot             = RB_Processed.hd_plot;
ym_mi_plot          = RB_Processed.ym_mi_plot;
fcs_pN_ave          = RB_Processed.fcs_pN_ave;                  
scs_pN_ave          = RB_Processed.scs_pN_ave; 

% and standard deviation (halved) for error bar plots
ym_rb_halfstd_plot = RB_Processed.ym_rb_halfstd_plot;
hd_rb_halfstd_plot = RB_Processed.hd_rb_halfstd_plot;

% and for force curves averaged by maximum indentation
stiff_pN_hmap_mi       = RB_Processed.hmapplot.stiff_pN_hmap_mi;
stiffhmap_pN_interp_mi = RB_Processed.hmapplot.stiffhmap_pN_interp_mi;
fcs_pN_ave_mi          = RB_Processed.fcs_pN_ave_mi;
scs_pN_ave_mi          = RB_Processed.scs_pN_ave_mi;

% and for the pores rotated individually
stiffhmap_pp_pN_interp_cell = PP_Processed.hmapplot.stiffhmap_pp_pN_interp_cell;
stiff_pp_pN_hmap_cell       = PP_Processed.hmapplot.stiff_pp_pN_hmap_cell;
th_pp_plot_0nm_cell         = PP_Processed.hmapplot.th_pp_plot_0nm_cell;
hd_pp_plot_0nm_cell         = PP_Processed.hmapplot.hd_pp_plot_0nm_cell;
ym_pp_plot_cell             = PP_Processed.ym_pp_plot;
hd_pp_plot                  = PP_Processed.hd_pp_plot;

% pull out the images (height and Hertz) for cropped pores
heightdata_nm_cell_cat      = matrices.heightdata_nm_cell_cat;
YM_MPa_cell_cat             = matrices.YM_MPa_cell_cat;
centres_cropped_xy_cat_cell = matrices.centres_cropped_xy_cat_cell;

%% Save out E_eff data on its own for plotting in another script if want

FullFileOutput = fullfile(OutputFolder, strcat(GenericSaveName, ' - YMRotationallyAveraged', '.mat'));

% save data out for the pores individually
YM_Plot.rv_plot            = rv_plot;
YM_Plot.ym_plot            = ym_plot;
YM_Plot.ym_rb_halfstd_plot = ym_rb_halfstd_plot;

save(FullFileOutput, 'YM_Plot');

%% plot for pores individually - subplot

numbpores = length(ym_pp_plot_cell);

if individual_subplots_on == 1

    for n = 1:numbpores

        stiffhmap_pp_pN_interp = stiffhmap_pp_pN_interp_cell{n};

        th_pp_plot_0nm         = th_pp_plot_0nm_cell{n};
        hd_pp_plot_0nm         = hd_pp_plot_0nm_cell{n};
        ym_pp_plot             = ym_pp_plot_cell{n};

        hd_cropped = heightdata_nm_cell_cat{n};
        ym_cropped = YM_MPa_cell_cat{n};
        cc_pp_xy   = centres_cropped_xy_cat_cell{n};

        % My custom colour map 'Blue_Yellow_Red' is saved to the Matlab structure
        % 'MyColourmaps'. First must load the structure, and then assign a variable
        % to the desired colour map saved in the structure.
        % MyClrMps = load('MyClrMps');
        % BYR_Colourmap = MyClrMps.BYR_ClrMp; 

        if custom_cmap == 1
            Structure = load(fullfile(OutputFolder, 'thiscolourmap.mat'));
            BYR_Colourmap = Structure.thiscmap;
        else
            Structure = load('MyColourMaps');
            BYR_Colourmap = Structure.MyColourMaps.BYR;
        end

        figure(),
        ax1 = subplot(221); imagesc(stiffhmap_pp_pN_interp, 'XData', [round(min(rv_plot)) round(max(rv_plot))], 'YData', [round(min(rv_plot)) round(max(rv_plot))]);
        hold on
        set(gca, 'FontSize', 14)
        set(gca, 'CLim', [cLow_hmap, cHigh_hmap]);
        colormap(ax1, BYR_Colourmap); 
        g = colorbar;
        ylabel(g,'Stiffness (pN/nm)', 'FontSize', 14)
        set(gca,'ydir','normal');
        hold on
        xlabel('Radial distance (nm)', 'FontSize', 14)
        ylabel('Tip-Sample Separation (nm)', 'FontSize', 14)
        title('Heat Map (pN/nm)', 'FontSize', 14)
        hold on
        plot(rv_plot, th_pp_plot_0nm, '--','Color', [0 0 0], 'LineWidth', 3)
        hold on
        h = area(rv_plot, hd_pp_plot_0nm, height_min); % the -80 specifies the baseline value
        h.FaceColor = [0.8 0.8 0.8];
        h.EdgeColor = [0.8 0.8 0.8];
        axis([x_radial_min x_radial_max height_min height_max])
        ax1.XColor = 'k';
        ax1.YColor = 'k';
        pbaspect([2 1 1])
        hold off

        ax2 = subplot(222); plot(rv_plot, ym_pp_plot, 'k', 'LineWidth', 2);
        axis([x_radial_min x_radial_max ym_pp_min ym_pp_max])
        xlabel('Radial distance (nm)', 'FontSize', 14)
        ylabel('E (MPa)', 'FontSize', 14)
        set(gca, 'FontSize', 14)
        title('Hertz model (MPa)', 'FontSize', 14);
        ax2.XColor = 'k';
        ax2.YColor = 'k';
        pbaspect([2 1 1])

        [rr, cc] = size(hd_cropped);

        ax3 = subplot(223); imagesc(hd_cropped);
        caxis([height_mat_clim_min height_mat_clim_max]);
        colormap(ax3, 'hot');
        colorbar;
        g = colorbar;
        ylabel(g,'Height (nm)', 'FontSize', 14)
        set(gca,'ydir','normal');
        set(gca, 'FontSize', 14)
        title('Height data (nm)', 'FontSize', 14)
        pbaspect([cc rr 1])
        ax3.XColor = 'k';
        ax3.YColor = 'k';
        hold on
        viscircles(cc_pp_xy, 1,'EdgeColor','w', 'LineWidth', 2);

        ax4 = subplot(224); imagesc(ym_cropped);
        caxis([ym_pp_min ym_pp_max])
        colormap(ax4, 'parula');
        colorbar;
        g = colorbar;
        ylabel(g,'E (MPa)', 'FontSize', 14)
        set(gca,'ydir','normal');
        set(gca, 'FontSize', 14)
        title('Hertz model (MPa)', 'FontSize', 14)
        pbaspect([cc rr 1])
        ax4.XColor = 'k';
        ax4.YColor = 'k';
        hold on
        viscircles(cc_pp_xy, 1,'EdgeColor','w', 'LineWidth', 2);
        set(gcf, 'Units', 'Normalized', 'OuterPosition', [0.2, 0.04, 0.8, 0.96]);
        hold off

        porenumb = n;
        
        fullfig_ppsavename = strcat(GenericSaveName, pp_figsavename, num2str(ym_pp_max), 'MPa_n=', num2str(porenumb));
        saveas(gca, fullfile(OutputFolder, fullfig_ppsavename), 'jpg');
        close(gcf);

    end

end

%% plot for pores individually - separate plots

if individual_plots_on == 1
    
    display('Plotting for individual pores...')

    for n = 1:numbpores

        stiffhmap_pp_pN_interp = stiffhmap_pp_pN_interp_cell{n};

        th_pp_plot_0nm         = th_pp_plot_0nm_cell{n};
        hd_pp_plot_0nm         = hd_pp_plot_0nm_cell{n};
        ym_pp_plot             = ym_pp_plot_cell{n};

        hd_cropped = heightdata_nm_cell_cat{n};
        ym_cropped = YM_MPa_cell_cat{n};
        cc_pp_xy   = centres_cropped_xy_cat_cell{n};

        % My custom colour map 'Blue_Yellow_Red' is saved to the Matlab structure
        % 'MyColourmaps'. First must load the structure, and then assign a variable
        % to the desired colour map saved in the structure.
        % MyClrMps = load('MyClrMps');
        % BYR_Colourmap = MyClrMps.BYR_ClrMp; 

        if custom_cmap == 1
            Structure = load(fullfile(OutputFolder, 'thiscolourmap.mat'));
            BYR_Colourmap = Structure.thiscmap;
        else
            Structure = load('MyColourMaps');
            BYR_Colourmap = Structure.MyColourMaps.BYR;
        end
    
        porenumb = n;
    
        % plot stiffness heat map
        set(groot,{'DefaultAxesXColor','DefaultAxesYColor','DefaultAxesZColor'},{'k','k','k'})

        figure(),
        imagesc(stiffhmap_pp_pN_interp, 'XData', [round(min(rv_plot)) round(max(rv_plot))], 'YData', [round(min(rv_plot)) round(max(rv_plot))]);
        set(gca,'linewidth',2)
        set(gca, 'FontSize', 12, 'Color', 'w')
        hold on
        set(gca, 'FontSize', 14)
        set(gca, 'CLim', [cLow_hmap, cHigh_hmap]);
        colormap(BYR_Colourmap); 
        g = colorbar;
        ylabel(g,'Stiffness (pN/nm)', 'FontSize', 14)
        set(gca,'ydir','normal');
        hold on
        xlabel('Radial distance (nm)', 'FontSize', 14)
        ylabel('Tip-Sample Separation (nm)', 'FontSize', 14)
        title('Heat Map (pN/nm)', 'FontSize', 14)
        hold on
        plot(rv_plot, th_pp_plot_0nm, '--','Color', [0 0 0], 'LineWidth', 3)
        hold on
        plot(rv_plot, hd_pp_plot_0nm, 'LineWidth', 2, 'Color', [0.8 0.8 0.8]); % the -80 specifies the baseline value
        axis([x_radial_min x_radial_max height_min height_max])
        pbaspect([2 1 1])
        set(gcf, 'Units', 'Normalized', 'OuterPosition', [0.2, 0.4, 0.5, 0.6]);
        
        fullfig_ppsavename_hm = strcat(GenericSaveName, pp_figsavename, 'heatmap_', num2str(porenumb));
        saveas(gca, fullfile(OutputFolder, fullfig_ppsavename_hm), 'pdf');
        close(gcf);


        % plot ym rotationally averaged
        set(groot,{'DefaultAxesXColor','DefaultAxesYColor','DefaultAxesZColor'},{'k','k','k'})
        
        figure(),
        plot(rv_plot, ym_pp_plot, 'Color', c1, 'LineWidth', 2);
        set(gca,'linewidth',2)
        set(gca, 'FontSize', 12, 'Color', 'w')
        axis([x_radial_min x_radial_max ym_pp_min ym_pp_max])
        xlabel('Radial distance (nm)', 'FontSize', 14)
        ylabel('E (MPa)', 'FontSize', 14)
        set(gca, 'FontSize', 14)
        title('Hertz model (MPa)', 'FontSize', 14);
        pbaspect([2 1 1])
        set(gcf, 'Units', 'Normalized', 'OuterPosition', [0.2, 0.4, 0.5, 0.6]);
        
        fullfig_ppsavename_ym = strcat(GenericSaveName, pp_figsavename, 'ym_', num2str(ym_pp_max), 'MPa_n=', num2str(porenumb));
        saveas(gca, fullfile(OutputFolder, fullfig_ppsavename_ym), 'pdf');
        close(gcf);

        % plot height img
        [rr, cc] = size(hd_cropped);
        set(groot,{'DefaultAxesXColor','DefaultAxesYColor','DefaultAxesZColor'},{'k','k','k'})
        
        figure(),
        imagesc(hd_cropped);
        set(gca,'linewidth',2)
        set(gca, 'FontSize', 12, 'Color', 'w')
        caxis([height_mat_clim_min height_mat_clim_max]);
        colormap('hot');
        colorbar;
        g = colorbar;
        ylabel(g,'Height (nm)', 'FontSize', 14)
        set(gca,'ydir','normal');
        set(gca, 'FontSize', 14)
        title('Height data (nm)', 'FontSize', 14)
        pbaspect([cc rr 1])
        hold on
        viscircles(cc_pp_xy, 1,'EdgeColor','w', 'LineWidth', 2);
        set(gcf, 'Units', 'Normalized', 'OuterPosition', [0.2, 0.4, 0.5, 0.6]);
        
        fullfig_ppsavename_img = strcat(GenericSaveName, pp_figsavename, 'img_', num2str(porenumb));
        saveas(gca, fullfile(OutputFolder, fullfig_ppsavename_img), 'pdf');
        close(gcf);
        
        % plot ym img
        set(groot,{'DefaultAxesXColor','DefaultAxesYColor','DefaultAxesZColor'},{'k','k','k'})
        
        figure(),
        imagesc(ym_cropped);
        set(gca,'linewidth',2)
        set(gca, 'FontSize', 12, 'Color', 'w')
        caxis([ym_pp_min ym_pp_max])
        colormap('bone');
        colorbar;
        g = colorbar;
        ylabel(g,'E_{eff} (MPa)', 'FontSize', 14)
        set(gca,'ydir','normal');
        set(gca, 'FontSize', 14)
        title('Hertz model (MPa)', 'FontSize', 14)
        pbaspect([cc rr 1])
        hold on
        set(gcf, 'Units', 'Normalized', 'OuterPosition', [0.2, 0.4, 0.5, 0.6]);
        
        fullfig_ppsavename_ymimg = strcat(GenericSaveName, pp_figsavename, 'ymimg_',num2str(ym_pp_max), 'MPa_n=',  num2str(porenumb));
        saveas(gca, fullfile(OutputFolder, fullfig_ppsavename_ymimg), 'pdf');
        close(gcf);
       
        
    end

end




%% Plot force curves and stiffness curves - aligned by contact point

% find index values for force/stiffness curves 20, 45, and 60 nm from
% centre

rv_half = rv_plot(round(length(rv_plot)/2):end);

rv_plot_abs = abs(rv_half - 20);
idx_20nm = find(rv_plot_abs == min(rv_plot_abs));

rv_plot_abs_45 = abs(rv_half - 45);
idx_45nm = find(rv_plot_abs_45 == min(rv_plot_abs_45));

rv_plot_abs_60 = abs(rv_half - 60);
idx_60nm = find(rv_plot_abs_60 == min(rv_plot_abs_60));


% aligned by contact point
figure(), 
ax1 = subplot(211);
plot(fcs_pN_ave{idx_60nm}(:,1), fcs_pN_ave{idx_60nm}(:,2), 'Color', c4, 'LineWidth', 2)
hold on
plot(fcs_pN_ave{idx_45nm}(:,1), fcs_pN_ave{idx_45nm}(:,2), 'Color', c3, 'LineWidth', 2)
hold on
plot(fcs_pN_ave{idx_20nm}(:,1), fcs_pN_ave{idx_20nm}(:,2), 'Color', c1, 'LineWidth', 2)
plot(fcs_pN_ave{1}(:,1), fcs_pN_ave{1}(:,2), 'Color', c2, 'LineWidth', 2)
legend('60 nm from centre', '45 nm from centre', '20 nm from centre', 'Centre', 'Location', 'northeast')
xlabel('Tip-sample separation (nm)', 'FontSize', 14)
ylabel('Force (pN)', 'FontSize', 14);
set(gca, 'FontSize', 14)
title(strcat('Averaged force curves (aligned by contact point) (pN) (n=', num2str(numbpores), ')'), 'FontSize', 16);
axis([fc_x_min fc_x_max fc_y_min fc_y_max])
ax1.XColor = 'k';
ax1.YColor = 'k';
pbaspect([2 1 1])


ax2 = subplot(212);
plot(scs_pN_ave{idx_60nm}(:,1), scs_pN_ave{idx_60nm}(:,2),'Color', c4, 'LineWidth', 2)
hold on
plot(scs_pN_ave{idx_45nm}(:,1), scs_pN_ave{idx_45nm}(:,2), 'Color', c3, 'LineWidth', 2)
hold on
plot(scs_pN_ave{idx_20nm}(:,1), scs_pN_ave{idx_20nm}(:,2), 'Color', c1, 'LineWidth', 2)
plot(scs_pN_ave{1}(:,1), scs_pN_ave{1}(:,2), 'Color', c2, 'LineWidth', 2)
legend('60 nm from centre', '45 nm from centre', '20 nm from centre', 'Centre', 'Location', 'northeast')
xlabel('Tip-sample separation (nm)', 'FontSize', 14)
ylabel('Stiffness (pN/nm)', 'FontSize', 14);
set(gca, 'FontSize', 14)
title(strcat('Averaged stiffness curves (aligned by contact point) (pN/nm) (n=', num2str(numbpores), ')'), 'FontSize', 16);
axis([fc_x_min fc_x_max sc_y_min sc_y_max])
ax2.XColor = 'k';
ax2.YColor = 'k';
pbaspect([2 1 1])
set(gcf, 'Units', 'Normalized', 'OuterPosition', [0.2, 0.04, 0.5, 0.96]);

fullfig_SCsavename = strcat(GenericSaveName, '_fc_sc_cp_', num2str(numbpores));
saveas(gca, fullfile(OutputFolder, fullfig_SCsavename), 'pdf');


%% Subplot averaged height and E_eff values - shaded error bars

figure(); 
set(groot,{'DefaultAxesXColor','DefaultAxesYColor','DefaultAxesZColor'},{'k','k','k'})
ax1 = subplot(211);
hold on
for i = 1:numbpores
    hd_pp_plot_array = hd_pp_plot{i};
    
    scatter(rv_plot, hd_pp_plot_array, 'MarkerEdgeColor', grey*[1 1 1], 'LineWidth', 0.5)
    
end
shadedErrorBar(rv_plot, hd_plot, hd_rb_halfstd_plot, 'lineprops', {'-','color', c2, 'linewidth', 1}, 'patchSaturation', 0.3)
xlabel('Radial distance (nm)', 'FontSize', 14)
ylabel('Height (nm)', 'FontSize', 14);
set(gca, 'FontSize', 14)
title(strcat('Height (nm) (n=', num2str(numbpores), ')'), 'FontSize', 16);
axis([-90 90 height_mat_clim_min height_mat_clim_max])
pbaspect([2 1 1])
ax1.XColor = 'k';
ax1.YColor = 'k';

ax2 = subplot(212);
hold on
for i = 1:numbpores
    ym_pp_plot = ym_pp_plot_cell{i};
    
    scatter(rv_plot, ym_pp_plot, 'MarkerEdgeColor', grey*[1 1 1], 'LineWidth', 0.5)
    
end
shadedErrorBar(rv_plot, ym_plot, ym_rb_halfstd_plot, 'lineprops', {'-','color',shadecolour, 'linewidth', 1}, 'patchSaturation', 0.3)
xlabel('Radial distance (nm)', 'FontSize', 14)
ylabel('E_{eff} (MPa)', 'FontSize', 14);
set(gca, 'FontSize', 14)
set(gca,'linewidth',1)
title(strcat('E_{eff} (MPa) (n=', num2str(numbpores), ')'), 'FontSize', 16);
axis([-90 90 ym_min ym_max])
pbaspect([2 1 1])
ax2.XColor = 'k';
ax2.YColor = 'k';

set(gcf, 'Units', 'Normalized', 'OuterPosition', [0.2, 0.04, 0.5, 0.96]);

fullfig_HertzModel_thinline_shaded_savename = strcat(GenericSaveName, '_subplot_height_E_eff_thinline_shaded_std_', num2str(ym_max), 'MPa_n=', num2str(numbpores));
saveas(gca, fullfile(OutputFolder, fullfig_HertzModel_thinline_shaded_savename), 'pdf');

%% Plot stiffness heatmaps

% get colourmap
if custom_cmap == 1
    Structure = load(fullfile(OutputFolder, 'thiscolourmap.mat'));
    BYR_Colourmap = Structure.thiscmap;
else
    Structure = load('MyColourMaps');
    BYR_Colourmap = Structure.MyColourMaps.BYR;
end

% Subplot stiffness heatmap and E_eff
figure(), 
ax1 = subplot(211); imagesc(stiffhmap_pN_interp, 'XData', [round(min(rv_plot)) round(max(rv_plot))], 'YData', [round(min(rv_plot)) round(max(rv_plot))]);
set(gca, 'FontSize', 14)
set(gca, 'CLim', [cLow_hmap, cHigh_hmap]);
colormap(ax1, BYR_Colourmap); 
set(gca,'ydir','normal');
hold on
xlabel('Radial distance (nm)', 'FontSize', 14)
ylabel('Tip-Sample Separation (nm)', 'FontSize', 14)
title(strcat('Stiffness heatmap (pN/nm) (n=', num2str(numbpores), ')'), 'FontSize', 16);
ax1.XColor = 'k';
ax1.YColor = 'k';
hold on
plot(rv_plot, th_plot_0nm, '--','Color', [0 0 0], 'LineWidth', 3)
hold on
h = area(rv_plot, hd_plot_0nm, height_min);
h.FaceColor = [0.8 0.8 0.8];
h.EdgeColor = [0.8 0.8 0.8];
axis([x_radial_min x_radial_max height_min height_max])
pbaspect([2 1 1])

ax2 = subplot(212);
hold on
for i = 1:numbpores
    ym_pp_plot = ym_pp_plot_cell{i};
    
    scatter(rv_plot, ym_pp_plot, 'MarkerEdgeColor', grey*[1 1 1], 'LineWidth', 0.5)
    
end
shadedErrorBar(rv_plot, ym_plot, ym_rb_halfstd_plot, 'lineprops', {'-','color',shadecolour, 'linewidth', 1}, 'patchSaturation', 0.3)
xlabel('Radial distance (nm)', 'FontSize', 14)
ylabel('E_{eff} (MPa)', 'FontSize', 14);
set(gca, 'FontSize', 14)
set(gca,'linewidth',1)
title(strcat('E_{eff} (MPa) (n=', num2str(numbpores), ')'), 'FontSize', 16);
axis([-90 90 ym_min ym_max])
pbaspect([2 1 1])
ax2.XColor = 'k';
ax2.YColor = 'k';

set(gcf, 'Units', 'Normalized', 'OuterPosition', [0.2, 0.04, 0.5, 0.96]);

linkaxes([ax1,ax2],'x')

fullfig_SCsavename = strcat(GenericSaveName, '_hm_ym_cp_', num2str(ym_max), 'MPa_n=', num2str(numbpores));
saveas(gca, fullfile(OutputFolder, fullfig_SCsavename), 'pdf');


%% Stiffness heatmap - contact point - fill - colourbar

set(groot,{'DefaultAxesXColor','DefaultAxesYColor','DefaultAxesZColor'},{'k','k','k'})

figure(), 
imagesc(stiffhmap_pN_interp, 'XData', [round(min(rv_plot)) round(max(rv_plot))], 'YData', [round(min(rv_plot)) round(max(rv_plot))]);
set(gca, 'FontSize', 14, 'Color', 'k')
set(gca, 'CLim', [cLow_hmap, cHigh_hmap]);
colormap(BYR_Colourmap); 
g = colorbar;
ylabel(g,'Stiffness (pN/nm)', 'FontSize', 14, 'Color', 'k')
set(gca,'ydir','normal');
hold on
xlabel('Radial distance (nm)', 'FontSize', 14, 'Color', 'k')
ylabel('Tip-Sample Separation (nm)', 'FontSize', 14, 'Color', 'k')
title(strcat('Stiffness heatmap (pN/nm) (n=', num2str(numbpores), ')'), 'FontSize', 16);
hold on
plot(rv_plot, th_plot_0nm, '--','Color', [0 0 0], 'LineWidth', 3)
hold on
h = area(rv_plot, hd_plot_0nm, height_min);
h.FaceColor = [0.8 0.8 0.8];
h.EdgeColor = [0.8 0.8 0.8];
set(gca,'Layer','top')
h.LineStyle = 'none';
axis([x_radial_min x_radial_max height_min height_max])
pbaspect([2 1 1])

fullfig_heatmap_cp_savename_norm = strcat(GenericSaveName, '_heatmap_cp_colourbar_matlabsave_', num2str(numbpores));
saveas(gca, fullfile(OutputFolder, fullfig_heatmap_cp_savename_norm), 'pdf');

%% Stiffness heatmap - contact point - fill - no colourbar

set(groot,{'DefaultAxesXColor','DefaultAxesYColor','DefaultAxesZColor'},{'k','k','k'})

figure(), 
imagesc(stiffhmap_pN_interp, 'XData', [round(min(rv_plot)) round(max(rv_plot))], 'YData', [round(min(rv_plot)) round(max(rv_plot))]);
set(gca, 'FontSize', 14, 'Color', 'k')
set(gca, 'CLim', [cLow_hmap, cHigh_hmap]);
colormap(BYR_Colourmap); 
set(gca,'ydir','normal');
hold on
xlabel('Radial distance (nm)', 'FontSize', 14, 'Color', 'k')
ylabel('Tip-Sample Separation (nm)', 'FontSize', 14, 'Color', 'k')
title(strcat('Stiffness heatmap (pN/nm) (n=', num2str(numbpores), ')'), 'FontSize', 16);
hold on
plot(rv_plot, th_plot_0nm, '--','Color', [0 0 0], 'LineWidth', 3)
hold on
h = area(rv_plot, hd_plot_0nm, height_min);
h.FaceColor = [0.8 0.8 0.8];
h.EdgeColor = [0.8 0.8 0.8];
set(gca,'Layer','top')
h.LineStyle = 'none';
axis([x_radial_min x_radial_max height_min height_max])
pbaspect([2 1 1])

fullfig_heatmap_nocb_cp_savename = strcat(GenericSaveName, '_heatmap_nocb_cp_colourbar_', num2str(numbpores));
saveas(gca, fullfile(OutputFolder, fullfig_heatmap_nocb_cp_savename), 'pdf');


%% Stiffness heatmap - contact point - no fill - colourbar

set(groot,{'DefaultAxesXColor','DefaultAxesYColor','DefaultAxesZColor'},{'k','k','k'})

figure(), 
imagesc(stiffhmap_pN_interp, 'XData', [round(min(rv_plot)) round(max(rv_plot))], 'YData', [round(min(rv_plot)) round(max(rv_plot))]);
set(gca, 'FontSize', 14, 'Color', 'k')
set(gca, 'CLim', [cLow_hmap, cHigh_hmap]);
colormap(BYR_Colourmap); 
g = colorbar;
ylabel(g,'Stiffness (pN/nm)', 'FontSize', 14, 'Color', 'k')
set(gca,'ydir','normal');
hold on
xlabel('Radial distance (nm)', 'FontSize', 14, 'Color', 'k')
ylabel('Tip-Sample Separation (nm)', 'FontSize', 14, 'Color', 'k')
title(strcat('Stiffness heatmap (pN/nm) (n=', num2str(numbpores), ')'), 'FontSize', 16);
hold on
plot(rv_plot, th_plot_0nm, '--','Color', [0 0 0], 'LineWidth', 3)
hold on
plot(rv_plot, hd_plot_0nm, 'LineWidth', 2, 'Color', [0.8 0.8 0.8]);
axis([x_radial_min x_radial_max height_min height_max])
pbaspect([2 1 1])

fullfig_heatmap_cp_noarea_norm = strcat(GenericSaveName, '_heatmap_cp_noarea_colourbar_', num2str(numbpores));
saveas(gca, fullfile(OutputFolder, fullfig_heatmap_cp_noarea_norm), 'pdf');


%% Stiffness heatmap - contact point - no fill - no colourbar

set(groot,{'DefaultAxesXColor','DefaultAxesYColor','DefaultAxesZColor'},{'k','k','k'})

figure(), 
imagesc(stiffhmap_pN_interp, 'XData', [round(min(rv_plot)) round(max(rv_plot))], 'YData', [round(min(rv_plot)) round(max(rv_plot))]);
set(gca, 'FontSize', 14, 'Color', 'k')
set(gca, 'CLim', [cLow_hmap, cHigh_hmap]);
colormap(BYR_Colourmap); 
set(gca,'ydir','normal');
hold on
xlabel('Radial distance (nm)', 'FontSize', 14, 'Color', 'k')
ylabel('Tip-Sample Separation (nm)', 'FontSize', 14, 'Color', 'k')
title(strcat('Stiffness heatmap (pN/nm) (n=', num2str(numbpores), ')'), 'FontSize', 16);
hold on
plot(rv_plot, th_plot_0nm, '--','Color', [0 0 0], 'LineWidth', 3)
hold on
plot(rv_plot, hd_plot_0nm, 'LineWidth', 2, 'Color', [0.8 0.8 0.8]);
axis([x_radial_min x_radial_max height_min height_max])
pbaspect([2 1 1])

fullfig_heatmap_cp_noarea_norm = strcat(GenericSaveName, '_heatmap_cp_noarea_nocolourbar_', num2str(numbpores));
saveas(gca, fullfile(OutputFolder, fullfig_heatmap_cp_noarea_norm), 'pdf');


%% this is for making the custom colourmap

% thiscmap = colormap(gca);
% 
% save(fullfile(OutputFolder, 'thiscolourmap.mat'), 'thiscmap');

