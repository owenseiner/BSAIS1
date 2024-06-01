%This script will be used to generate a figure with model predicted basal
%temp from all 9 models
clear all
close all
k = 8.7e-4;
BM3_mask_ais = load("BM3_mask_ais.mat").BM3_mask_ais;


%Borehole verification for model & agreement
lat = [-66.033333; -80.016666; -75.1;
    -70.5; -66.76; -81.658; -90; -79.46765;
    -75.1];
lon = [-64.066666; -119.516666; 39.70;
    -65; 112.80; -148.81; 0; -112.08562; 123.4];
%1 = melted, 0 = frozen
basal_therm = [0; 1; 0; 0; 0; 0; 1; 1; 1];

x_coords = zeros(height(lon), 1);
y_coords = zeros(height(lat), 1);

%Assuming ll2xy function takes latitude and longitude and returns x and y
%Replace 'hemisphere' with appropriate hemisphere indicator if necessary
for i = 1:height(lat)
    [x_coords(i), y_coords(i)] = ll2xy(lat(i), lon(i), -1);
end



litempbotgr = ["litempbotgr_AIS_AWI_PISM1_hist_std.nc", ...
    "litempbotgr_AIS_DOE_MALI_ctrl_proj_std.nc", ...
    "litempbotgr_AIS_ILTS_PIK_SICOPOLIS_hist_std.nc", ...
    "litempbotgr_AIS_LSCE_GRISLI2_hist_std.nc", ...
    "litempbotgr_AIS_NCAR_CISM_hist_std.nc", ...
    "litempbotgr_AIS_PIK_PISM1_hist_open.nc", ...
    "litempbotgr_AIS_ULB_fETISh_ctrl_proj_std.nc", ...
    "litempbotgr_AIS_VUB_AISMPALEO_hist_std.nc", ...
    "litempbotgr_AIS_VUW_PISM_hist_open.nc"];

lithk = ["lithk_AIS_AWI_PISM1_hist_std.nc",...
    "lithk_AIS_DOE_MALI_ctrl_proj_std.nc",...
    "lithk_AIS_ILTS_PIK_SICOPOLIS_hist_std.nc",...
    "lithk_AIS_LSCE_GRISLI2_hist_std.nc",...
    "lithk_AIS_NCAR_CISM_hist_std.nc",...
    "lithk_AIS_PIK_PISM1_hist_open.nc",...
    "lithk_AIS_ULB_fETISh_ctrl_proj_std.nc",...
    "lithk_AIS_VUB_AISMPALEO_hist_std.nc",...
    "lithk_AIS_VUW_PISM_hist_open.nc"]; 

compileData = zeros(761);
numSteps = 256;


% Define the temperature color map breakpoints
temperature_breakpoints = [-10, -8, -6, -4, -3, -1.5, -1, -0.5];

% Define the corresponding colors
temperature_colors = [0 0.2 0.6;    % Darker Blue (for values - 10 to -8)
                      0 0 0.5;      % Dark Blue (for values -8 to -6)
                      0 0 0.5;      % Dark Blue (for -6 to -4)
                      0 0.4 0.8;    % Medium Dark Blue (for values -4 to -3)
                      0 0.5 1;      % Light Blue (for values -3 to -1.5)
                      1 1 0.8;      % Pale Yellow (for values -1.5 to -1)
                      1 0.7 0;      % Light Red (for values -1 to -0.5)
                      1 0 0];        % Dark Red (for values above -0.5)

% Generate a finer set of breakpoints for interpolation
fine_breakpoints = linspace(min(temperature_breakpoints), max(temperature_breakpoints), 500);

% Interpolate colors at the finer breakpoints
cmap = interp1(temperature_breakpoints, temperature_colors, fine_breakpoints);

    
model_names = ["AWI\_PISM" "DOE\_MALI" "ILTS\_PIK" "LSCE\_Grisli" "NCAR\_CISM" "PIK\_PISM" "ULB\_f.ETISH" "VUB\_AISMPALEO" "VUW\_PISM"];
figure(1);
set(gcf,'position',[314  222 726 535]);
for i = 1:length(litempbotgr)
    
    % Saves basal temp data
    temps = rot90(ncread(litempbotgr{i},'litempbotgr'));
    % Saves time data
    time = ncread(litempbotgr{i},'time');
    % Saves thickness data
    thicknesses = rot90(ncread(lithk{i},'lithk'));
    
    if (contains(litempbotgr(i),'hist'))
        temp = temps(:,:,length(time));
        temp(temp == 0) = NaN;
        thickness = thicknesses(:,:,length(time));
        celsius = temp - 273.15;
        pressureMeltingPoint = -(k * thickness);
        aboveBelowMelt = celsius - pressureMeltingPoint;
        aboveBelowMelt(BM3_mask_ais == 2) = NaN;
    
    
    elseif (contains(litempbotgr(i),'ctrl'))
        temp = temps(:,:,1);
        temp(temp == 0) = NaN;
        thickness = thicknesses(:,:,1);
        celsius = temp - 273.15;
        pressureMeltingPoint = -(k * thickness);
        aboveBelowMelt = celsius - pressureMeltingPoint;
        aboveBelowMelt(BM3_mask_ais == 2) = NaN;
    end
    
    cmin = min(aboveBelowMelt(~isnan(aboveBelowMelt(:))));
    cmax = max(aboveBelowMelt(~isnan(aboveBelowMelt(:))));

    x_8km = -3040:8:3040; %in km
    y_8km = -3040:8:3040; %in km
    
    %Create a subplot in a 3x3 grid and activate it
    col = mod(i-1,3);
    row = double( idivide(i-1,int16(3)));
    %Define margins
    margin_left = 0.0;  % Left margin
    margin_bottom = 0.0;  % Bottom margin
    spacing = 0.0;  % Space between plots
    %Calculate positions
    width = (1 ) / 3.2;
    height = (1 ) / 3.2;
    left =  col * width;
    bottom = 0.97 - (row + 1) * height;
    subplot(3,3,i)
    set(gca,'position',[left, bottom, width, height]);
    %set(gca,'position',[0.31*(col) 0.31*(row) 0.3 0.3]);
    colormap(cmap);
    set(gcf,'color','w');
    %aboveBelowMelt(no_ice_mask_subglacial) = NaN;
    nanIdx = isnan(aboveBelowMelt);
    h = imagesc(x_8km * 1e3, y_8km *1e3 , aboveBelowMelt);
    set(h, 'AlphaData', ~nanIdx); 
    hold on
    for ind = 1:length(basal_therm)
        if basal_therm(ind) == 1
            scatter(x_coords(ind), y_coords(ind) * -1, 50, ... % Increased size
                'filled', 'MarkerEdgeColor', [0 0 0], ...
                'MarkerFaceColor', [1 0 0], 'LineWidth', 2); % Increased outline width
        elseif basal_therm(ind) == 0
                 scatter(x_coords(ind), y_coords(ind) * -1, 50, ... % Increased size
                'filled', 'MarkerEdgeColor', [0 0 0], ...
                'MarkerFaceColor', [30/255, 144/255, 255/255], 'LineWidth', 2); % Dodger Blue color and increased outline width
        end
    end
    clim([-10, 0])
    nanIdx = isnan(aboveBelowMelt);
    aboveBelowMelt(nanIdx) = cmin - (cmax - cmin); %Set NaN values to a value outside the color range
    text(0,-2.4*10^6, model_names(i), 'HorizontalAlignment', 'center', 'FontSize', 15, 'FontWeight', 'bold', 'FontName', 'Sans Serif')
    hold on
    axis off
    x0 = -2.7E6;
    y0 = 2.1E6;
    lengthscale = 1E6;
    widthscale = 50000;
    fontcolor = 'black';
    if i == 9
        %Draw the scale bar
         p = patch([x0, x0 + lengthscale, x0 + lengthscale, x0], ...
          [y0, y0, y0 + widthscale, y0 + widthscale], fontcolor, ...
          'EdgeColor', fontcolor);
        %Add the scale bar label
         text(x0 + lengthscale/2, 1.90E6, sprintf('%.1f km', lengthscale/1E3), ...
        'HorizontalAlignment', 'center', 'Color', fontcolor, 'FontSize', 16);
    end
    axis off
    axis equal
    xlim([-3 3]*10^6);
    ylim([-2.5 2.5]*10^6);
end
cbar = colorbar('Position', [0.92, 0.1, 0.02, 0.8]); %[left bottom width height]
cbar.FontSize = 12;
cbar.Label.String = ['Pressure-Adjusted Temperature',char(176), 'C'];
cbar.Label.FontSize = 18;
cbar.Label.FontName = 'Sans Serif';
