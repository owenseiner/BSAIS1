%Generate agreement maps between methods
%Load subglacial routing 
%2 = frozen
%1 = melted
%0 = ocean/ice free land

x_8km = -3040:8:3040; %in km
y_8km = -3040:8:3040; %in km


%Import the active lakes from Livingstone et al and convert to table
lakesData2 = load('ActiveSubglacialLivingstone.mat');
%lakesData2 = readmatrix("Seven_Active_Lakes.xlsx");
lakesData2 = lakesData2.ActiveSubglacialLakesOnlyLivingstone;

%Import the stable lakes from Livingstone et al and convert to table
lakes_stable = load('StableSubglacialLivingstone.mat');
lakes_stable = lakes_stable.StableSubglacialLivingstone;

subglacial_binary_8km = load("binaryLakesRouting.mat").binaryLakesRouting;
BM3_mask = load("BM3Mask_AIS_8km.mat").downsampledBM3Mask;
x_BM3 = linspace(-3333, 3333, 761);
y_BM3 = linspace(-3333, 3333, 761);

[BM3_x, BM3_y] = meshgrid(x_BM3, y_BM3);
[x_target, y_target] = meshgrid(x_8km, y_8km);

BM3_mask_interp = interp2(BM3_x, BM3_y, BM3_mask, x_target, y_target, 'linear');

%subglacial_binary_8km_interp = interp2(x_BM3, y_BM3, subglacial_binary_8km, x_target, y_target, 'linear');
subglacial_binary_categorized = subglacial_binary_8km;
subglacial_binary_categorized(BM3_mask_interp == 0) = 0;
subglacial_binary_categorized(subglacial_binary_8km == 1) = 1;
subglacial_binary_categorized(subglacial_binary_8km == 0 & BM3_mask_interp == 2) = 2;
subglacial_binary_categorized(BM3_mask_interp == 1 | BM3_mask_interp == 3) = 3;
subglacial_binary_cateogrized(BM3_mask_interp == 0 | BM3_mask_interp == 1) = 0;

subglacial_binary_figure = zeros([761,761]);
subglacial_binary_figure(subglacial_binary_8km == 1) = 1;
%subglacial_binary_figure(BM3_mask_interp == 0 | BM3_mask_interp );


customColormap = [

    1.0 1.0 1.0;                    % 0: Ocean and Ice Free Land (White)
    190/255 0 0;                    % 1: Stream & lakes (Red)
    153/255 204/255 255/255;        % 2: Not stream & lakes ice (Light blue)
    0.5 0.5 0.5                     % 3: Ice shelves
];

%Create tables to store x and y (in m) of lakes on a 13333x13333 grid
x_excel_active = zeros(172, 1);
y_excel_active = zeros(172,1);
x_excel_stable = zeros(503, 1);
y_excel_stable = zeros(503,1);

%Convert lakes tables to arrays with only lat and lon
lakesTable2 = table2array(lakesData2(1:172, 3:4));
%lakesTable2 = lakesData2(1:8, 3:4);
lakesTableStable = table2array(lakes_stable(1:503, 3:4));

%Convert from lat and lon to x and y
%Use -1 to indicate southern hemisphere
%y is latitude (1st column), x is longitude (2nd column)
[x_excel_active, y_excel_active] = ll2xy(lakesTable2(:,1), lakesTable2(:, 2), -1);
[x_excel_stable, y_excel_stable] = ll2xy(lakesTableStable(:,1), lakesTableStable(:, 2), -1);

figure
imagesc(x_8km * 1e3, y_8km * 1e3, subglacial_binary_categorized)
colormap(customColormap)
hold on
%Plot the active and stable lakes
h1 = scatter(x_excel_active, y_excel_active * -1, 50, 'filled', 'MarkerEdgeColor', [1 1 1], 'MarkerFaceColor', [134/255 104/255 71/255], 'LineWidth', 1); % Decreased size
h2 = scatter(x_excel_stable, y_excel_stable * -1, 50, 'filled', 'MarkerEdgeColor', [1 1 1], 'MarkerFaceColor', [1 204/255 0], 'LineWidth', 1); % Decreased size
% Dummy plots for the map colors
h3 = plot(nan, nan, 's', 'MarkerEdgeColor', 'none', 'MarkerFaceColor', [190/255 0 0], 'MarkerSize', 15); % Blue square for legend
h4 = plot(nan, nan, 's', 'MarkerEdgeColor', 'none', 'MarkerFaceColor', [153/255 204/255 255/255], 'MarkerEdgeColor', [153/255 204/255 255/255], 'MarkerSize', 15); % Light blue
axis equal
axis off
x0 = -2.7E6;
y0 = 2.1E6;
lengthscale = 1E6;
widthscale = 50000;
fontcolor = 'black';
%Draw the scale bar
patch([x0, x0 + lengthscale, x0 + lengthscale, x0], ...
      [y0, y0, y0 + widthscale, y0 + widthscale], fontcolor, ...
      'EdgeColor', fontcolor);
%Add the scale bar label
text(x0 + lengthscale/2, 1.95E6, sprintf('%.1f km', lengthscale/1E3), ...
     'HorizontalAlignment', 'center', 'Color', fontcolor, 'FontSize', 16);
axis off
set(gcf,'color','w');
%Add the legend and position it closer to the plot
legend([h1, h2, h3, h4], {'Active Lakes', 'Stable Lakes', 'Presumed Thawed', 'Presumed Frozen'}, 'FontSize', 14, 'FontName', 'Sans Serif');
%Adjust the position of the legend
legend('boxoff'); % Remove the box around the legend
lgd = legend;
%lgd.Position = [0.75, 0.6, 0.1, 0.1]; % Custom position [left, bottom, width, height] % Custom position [left, bottom, width, height]
lgd.Position = [0.2, 0.25, lgd.Position(3), lgd.Position(4)];

figure
imagesc(x_8km * 1e3, y_8km * 1e3, subglacial_binary_categorized)
colormap(customColormap)
hold on
%Plot the active and stable lakes
h1 = scatter(x_excel_active, y_excel_active * -1, 50, 'filled', 'MarkerEdgeColor', [1 1 1], 'MarkerFaceColor', [134/255 104/255 71/255], 'LineWidth', 1); % Decreased size
h2 = scatter(x_excel_stable, y_excel_stable * -1, 50, 'filled', 'MarkerEdgeColor', [1 1 1], 'MarkerFaceColor', [1 204/255 0], 'LineWidth', 1); % Decreased size
% Dummy plots for the map colors
h3 = plot(nan, nan, 's', 'MarkerEdgeColor', 'none', 'MarkerFaceColor', [190/255 0 0], 'MarkerSize', 15); % Blue square for legend
h4 = plot(nan, nan, 's', 'MarkerEdgeColor', 'none', 'MarkerFaceColor', [153/255 204/255 255/255], 'MarkerEdgeColor', [153/255 204/255 255/255], 'MarkerSize', 15); % Light blue
axis equal
axis off
x0 = -2.7E6;
y0 = 2.1E6;
lengthscale = 1E6;
widthscale = 50000;
fontcolor = 'black';
%Draw the scale bar
patch([x0, x0 + lengthscale, x0 + lengthscale, x0], ...
      [y0, y0, y0 + widthscale, y0 + widthscale], fontcolor, ...
      'EdgeColor', fontcolor);
%Add the scale bar label
text(x0 + lengthscale/2, 1.95E6, sprintf('%.1f km', lengthscale/1E3), ...
     'HorizontalAlignment', 'center', 'Color', fontcolor, 'FontSize', 16);
axis off
set(gcf,'color','w');
%Add the legend and position it closer to the plot
legend([h1, h2, h3, h4], {'Active Lakes', 'Stable Lakes', 'Presumed Thawed', 'Presumed Frozen'}, 'FontSize', 14, 'FontName', 'Sans Serif');
%Adjust the position of the legend
legend('boxoff'); % Remove the box around the legend
lgd = legend;
lgd.ItemTokenSize = [60, 30]; % Custom position [left, bottom, width, height] % Custom position [left, bottom, width, height]
%title("Binary Basal Hydrology", 'FontSize', 16, 'FontWeight', 'bold', 'FontName', 'Sans Serif')


%title("Binary Basal Hydrology", 'FontSize', 16, 'FontWeight', 'bold', 'FontName', 'Sans Serif')


%All maps to be summed for agreement should have legend
%NaN = ocean/ice free land
%0 = frozen
%1 = melted

%Load lake_vostok_mask
%lake_vostok = load("lake_vostok_mask.mat").lake_vostok_mask;

%Original
%1 = frozen
%0 = melted
%NaN = ocean/ice free land
ismip6_binary_8km = load("ismip6_8km.mat").ismip6_8km;
agreement = ismip6_binary_8km + subglacial_binary_8km;
%Load Ant Smoothing binary agreement map
%BM2_mask_ais = load("BM2Mask_AIS_8km.mat").downsampledBM2Mask;

%For this upload, the legend is
%1 = melted
%0 = frozen
%NaN = ocean/ice free land
downsampledBM3Mask = load("BM3_mask_interp.mat").BM3_mask_interp;
ant_smoothing_binary_8km = load("binaryAntSmoothing_ColdIce.mat").binaryRatio_8km;
ant_smoothing_calc = ant_smoothing_binary_8km;
ant_smoothing_calc(downsampledBM3Mask == 0) = 2;

ant_smoothing_binary_8km(isnan(ismip6_binary_8km)) = NaN;
agreement3 = ant_smoothing_binary_8km + agreement;
agreement2 = ant_smoothing_binary_8km + ismip6_binary_8km;

indx_vostok = load("indx_vostok.mat").indx_vostok_ismip;
indy_vostok = load("indy_vostok.mat").indy_vostok_ismip;
indx_shelves = load("indx_shelves.mat").indx_shelves_ismip;
indy_shelves = load("indy_shelves.mat").indy_shelves_ismip;

BM3_mask_ais = zeros(761,761);
%1 = vostok
%2 = ice shelves
for i = 1:length(indx_vostok)
    BM3_mask_ais(indy_vostok(i), indx_vostok(i)) = 1;
end

for i = 1:length(indx_shelves)
    BM3_mask_ais(indy_shelves(i), indx_shelves(i)) = 2;
end

BM3_mask_ais = flipud(BM3_mask_ais);

save("BM3_mask_ais.mat", "BM3_mask_ais");

%Make lake vostok == 3 melted
agreement3(BM3_mask_ais == 1) = 3;
agreement2(BM3_mask_ais == 1) = 2;
%Add ice shelves as 4
agreement3(BM3_mask_ais == 2) = 4;
agreement2(BM3_mask_ais == 2) = 3;
%Set NaN to 5
agreement2(isnan(agreement2)) = 4;


%All maps should have legend
%NaN = ocean/ice free land
%0 = frozen
%1 = melted

%1 = frozen
%2 = ocean/ice free land
%Change this legend to match with ismip6/subglacial
%2 = frozen
%1 = melted
%0 = ocean/ice free land
%ant_smoothing_binary_8km(lake_vostok) = 1;

%Agreement 2 has values
%NaN = ocean/ice free land
%0 = frozen
%1 = 2/3 frozen
%2 = 2/3 melted
%3 = melted
%4 = ice shelf
%5 = ocean/ice free land

%Borehole verification for model & agreement
lat = [-66.033333; -80.016666; -75.1;
    -70.5; -66.76; -81.658; -90; -79.46765;
    -75.1];
lon = [-64.066666; -119.516666; 39.70;
    -65; 112.80; -148.81; 0; -112.08562; 123.4];
%1 = melted, 0 = frozen
basal_therm = [0; 1; 0; 0; 0; 0; 1; 1; 1];

x_coords = zeros(9, 1);
y_coords = zeros(9, 1);

%Assuming ll2xy function takes latitude and longitude and returns x and y
%Replace 'hemisphere' with appropriate hemisphere indicator if necessary
for i = 1:9
    [x_coords(i), y_coords(i)] = ll2xy(lat(i), lon(i), -1);
end


customColorMap_3 = [
    0 0.2 1;        % Blue for 0 (3/3 frozen)
    0.4 0.8 1;      % Sky Blue for 1 (2/3 frozen)
    1 0.8 0.8;      % Pink for 2 (2/3 melted)
    1 0 0;          % Red for 3 (melted)
    0.5 0.5 0.5;    % Gray for 4
];

customColorMap_2 = [
    0 0.2 1;        % Blue for 0 (2/2 frozen)
      1 204/255 0;        % Yellow for 1 (uncertain)
      1 0 0;        % Red for 2 (2/2 melted)
   0.5 0.5 0.5;     % Gray for 3
      1 1 1;        % White for 4
];

customColor_BinaryModel = [
    0 0.2 1;       %Blue for 0
    1 0 0;         %Red for 1
    0.5 0.5 0.5;   %Grey for 2
    1 1 1;         %White for 3
];

x_8km = -3040:8:3040; %in km
y_8km = -3040:8:3040; %in km

ismip6_binary_8km(isnan(ismip6_binary_8km)) = 3;
ismip6_binary_8km(BM3_mask_ais == 2) = 2;

figure
imagesc(x_8km * 1e3, y_8km * 1e3, ismip6_binary_8km)
colormap(customColor_BinaryModel)
hold on
for i = 1:length(basal_therm)
    if basal_therm(i) == 1
        scatter(x_coords(i), y_coords(i) * -1, 'filled', ...
                'MarkerEdgeColor', [1 1 1], 'MarkerFaceColor', [1 0 0], 'LineWidth', 2);
    elseif basal_therm(i) == 0
        scatter(x_coords(i), y_coords(i) * -1, 'filled', ...
                'MarkerEdgeColor', [1 1 1], 'MarkerFaceColor', [0 0 1], 'LineWidth', 2);
    end
end
axis equal
x0 = -2.7E6;
y0 = 2.1E6;
lengthscale = 1E6;
widthscale = 50000;
fontcolor = 'black';
hold on
%scatter(x_coords, y_coords * -1, 'filled', 'MarkerEdgeColor', [0 0 0], 'MarkerFaceColor', [0 0 0]);
%Draw the scale bar
p = patch([x0, x0 + lengthscale, x0 + lengthscale, x0], ...
          [y0, y0, y0 + widthscale, y0 + widthscale], fontcolor, ...
          'EdgeColor', fontcolor);
%Add the scale bar label
text(x0 + lengthscale/2, 1.95E6, sprintf('%.1f km', lengthscale/1E3), ...
     'HorizontalAlignment', 'center', 'Color', fontcolor);
axis off
set(gcf,'color','w');

%Plot agreement map 2
figure
imagesc(x_8km * 1e3, y_8km * 1e3, agreement2)
hold on
%Plot the boreholes
h1 = scatter(x_coords(basal_therm == 1), y_coords(basal_therm == 1) * -1, 100, 'filled', 'MarkerEdgeColor', [0 0 0], 'MarkerFaceColor', [1 0 0], 'LineWidth', 2); % Red points
h2 = scatter(x_coords(basal_therm == 0), y_coords(basal_therm == 0) * -1, 100, 'filled', 'MarkerEdgeColor', [0 0 0], 'MarkerFaceColor', [30/255, 144/255, 255/255], 'LineWidth', 2); % Dodger blue points
% Create dummy scatter plots for legend items representing colormap categories
h3 = plot(nan, nan, 's', 'MarkerEdgeColor', 'none', 'MarkerFaceColor', customColorMap_2(1,:), 'MarkerSize', 15); 
h4 = plot(nan, nan, 's', 'MarkerEdgeColor', 'none', 'MarkerFaceColor', customColorMap_2(2,:), 'MarkerSize', 15);
h5 = plot(nan, nan, 's', 'MarkerEdgeColor', 'none', 'MarkerFaceColor', customColorMap_2(3,:), 'MarkerSize', 15); 
h6 = plot(nan, nan, 's', 'MarkerEdgeColor', 'none', 'MarkerFaceColor', customColorMap_2(4,:), 'MarkerSize', 15); % Gray square for legend
colormap(customColorMap_2)
axis equal
x0 = -2.7E6;
y0 = 2.1E6;
lengthscale = 1E6;
widthscale = 50000;
fontcolor = 'black';
%Draw the scale bar
patch([x0, x0 + lengthscale, x0 + lengthscale, x0], ...
      [y0, y0, y0 + widthscale, y0 + widthscale], fontcolor, ...
      'EdgeColor', fontcolor);
%Add the scale bar label
text(x0 + lengthscale/2, 1.95E6, sprintf('%.1f km', lengthscale/1E3), ...
     'HorizontalAlignment', 'center', 'Color', fontcolor, 'FontSize', 14);
axis off
set(gcf, 'color', 'w');
%Add the legend
lgd = legend([h1, h2, h3, h4, h5, h6], {'Thawed Boreholes', 'Frozen Boreholes', 'Presumed Frozen', 'Uncertain', 'Presumed Thawed', 'Ice Shelves'}, 'Location', 'northeastoutside', 'FontSize', 14, 'FontName', 'Sans Serif');
%lgd.Position = [0.65, 0.7, 0.2, 0.2];
lgd.Position = [0.15, 0.2, lgd.Position(3), lgd.Position(4)];
lgd.Box = 'off';
%lgd.ItemTokenSize = [200, 18]; % Adjust these values as needed

figure
h = imagesc(x_8km .* 1e3, y_8km .* 1e3, agreement3);
set(h, 'AlphaData', ~isnan(agreement3));  % Make NaNs transparent
hold on
h1 = scatter(x_coords(basal_therm == 1), y_coords(basal_therm == 1) * -1, 100, 'filled', 'MarkerEdgeColor', [0 0 0], 'MarkerFaceColor', [1 0 0], 'LineWidth', 2); % Red points
h2 = scatter(x_coords(basal_therm == 0), y_coords(basal_therm == 0) * -1, 100, 'filled', 'MarkerEdgeColor', [0 0 0], 'MarkerFaceColor', [30/255, 144/255, 255/255], 'LineWidth', 2); % Dodger blue points
h3 = plot(nan, nan, 's', 'MarkerEdgeColor', 'none', 'MarkerFaceColor', customColorMap_3(1,:), 'MarkerSize', 15); % Blue square for legend
h4 = plot(nan, nan, 's', 'MarkerEdgeColor', 'none', 'MarkerFaceColor', customColorMap_3(2,:), 'MarkerSize', 15); % Light blue square for legend
h5 = plot(nan, nan, 's', 'MarkerEdgeColor', 'none', 'MarkerFaceColor', customColorMap_3(3,:), 'MarkerSize', 15); % Light red square for legend
h6 = plot(nan, nan, 's', 'MarkerEdgeColor', 'none', 'MarkerFaceColor', customColorMap_3(4,:), 'MarkerSize', 15); % Red square for legend
h7 = plot(nan, nan, 's', 'MarkerEdgeColor', 'none', 'MarkerFaceColor', customColorMap_3(5,:), 'MarkerSize', 15); % Gray square for legend
colormap(customColorMap_3)
axis equal
x0 = -2.7E6;
y0 = 2.1E6;
lengthscale = 1E6;
widthscale = 50000;
fontcolor = 'black';
%scatter(x_coords, y_coords * -1, 'filled', 'MarkerEdgeColor', [0 0 0], 'MarkerFaceColor', [0 0 0]);
%Draw the scale bar
p = patch([x0, x0 + lengthscale, x0 + lengthscale, x0], ...
          [y0, y0, y0 + widthscale, y0 + widthscale], fontcolor, ...
          'EdgeColor', fontcolor);
%Add the scale bar label
text(x0 + lengthscale/2, 1.95E6, sprintf('%.1f km', lengthscale/1E3), ...
     'HorizontalAlignment', 'center', 'Color', fontcolor, 'FontSize', 14);
axis off
lgd = legend([h1, h2, h3, h4, h5, h6, h7], {'Thawed Boreholes', 'Frozen Boreholes', 'Presumed Frozen', '2/3 Frozen', '2/3 Thawed', 'Presumed Thawed', 'Ice Shelves'}, 'Location', 'northeastoutside', 'FontSize', 14, 'FontName', 'Sans Serif');
%lgd.Position = [0.65, 0.7, 0.2, 0.2];
lgd.Position = [0.15, 0.22, lgd.Position(3), lgd.Position(4)];
lgd.Box = 'off';
set(gcf,'color','w');


totalNumPixels = sum((agreement3(:) == 0 | agreement3(:) == 1 | agreement3(:) == 2 | agreement3(:) == 3));
frozen_num = sum(agreement3(:) == 0);
uncertain_num = sum((agreement3(:) == 1) | (agreement3(:) == 2));
melted_num = sum(agreement3(:) == 3);

melted_basalslip = sum(ant_smoothing_calc(:) == 1);
frozen_basalslip = sum(ant_smoothing_calc(:) == 0);
totalice_basalslip = melted_basalslip + frozen_basalslip;
melted_basalslip_percent = melted_basalslip * 100 /totalice_basalslip;
frozen_basalslip_percent = frozen_basalslip * 100/totalice_basalslip;

melted_subglacial = sum(subglacial_binary_categorized(:) == 1);
frozen_subglacial = sum(subglacial_binary_categorized(:) == 2);
total_subglacial = melted_subglacial + frozen_subglacial;
melted_subglacial_percent = melted_subglacial * 100/total_subglacial;
frozen_subglacial_percent = frozen_subglacial * 100/total_subglacial;

totalNumPixels_2 = sum(agreement2(:) == 0 | agreement2(:) == 1 | agreement2(:) == 2);
melted_num_2 = sum(agreement2(:) == 2);
melted_percent = melted_num_2*100/totalNumPixels_2;
melted_percent_3 = melted_num * 100/totalNumPixels;
frozen_percent = frozen_num *100/ totalNumPixels;
uncertain_percent = 100 - frozen_percent - melted_percent;

disp(frozen_percent)
disp(uncertain_percent)
disp(melted_percent)

disp(melted_basalslip_percent)
disp(frozen_basalslip_percent)

disp(melted_subglacial_percent)
disp(frozen_subglacial_percent)
