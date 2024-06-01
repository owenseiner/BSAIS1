%This script constructs a matrix of -1 1 values to determine whether the
%each model predicts that the AIS is frozen or melted based on
%MacGregor's 2016 paper on Greenland
clearvars
clc
clf
k = 8.7e-4;
clear binaryImage
clear compileData

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

%Borehole verification for model & agreement
lat = [-66.033333; -80.016666; -75.1;
    -70.5; -66.76; -81.658; -90; -79.46765;
    -75.1];
lon = [-64.066666; -119.516666; 39.70;
    -65; 112.80; -148.81; 0; -112.08562; 123.4];
%1 = melted, 0 = frozen
basal_therm = [0; 1; 0; 0; 0; 0; 1; 1; 1];

x_coords_verif = zeros(height(lon), 1);
y_coords_verif = zeros(height(lat), 1);


%Assuming ll2xy function takes latitude and longitude and returns x and y
%Replace 'hemisphere' with appropriate hemisphere indicator if necessary
for i = 1:height(lat)
    [x_coords_verif(i), y_coords_verif(i)] = ll2xy(lat(i), lon(i), -1);
end


x_8km = -3040:8:3040; %in km
y_8km = -3040:8:3040; %in km

x_BM3 = linspace(-3333, 3333, 761);
y_BM3 = linspace(-3333, 3333, 761);

[BM3_x, BM3_y] = meshgrid(x_BM3, y_BM3);
[x_target, y_target] = meshgrid(x_8km, y_8km);

compileData = zeros(761);

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
        
        % Populate the binary matrix
        aboveBelowMelt(aboveBelowMelt > -1) = 1; %1 is melted
        aboveBelowMelt(aboveBelowMelt <= -1) = 0; %-1 is frozen
        compileData = compileData + aboveBelowMelt;
        
        %{
        figure(i)
            imagesc(aboveBelowMelt);
            colorbar;
            caxis([-3,3]);
            title(strrep(extractAfter(litempbotgr{i},'litempbotgr'),'_',' '));
            xlabel('meters');
            ylabel('meters');
        %}
            
     elseif (contains(litempbotgr(i),'ctrl'))
        temp = temps(:,:,1);
        temp(temp == 0) = NaN;
        thickness = thicknesses(:,:,1);
        celsius = temp - 273.15;
        pressureMeltingPoint = -(k * thickness);
        aboveBelowMelt = celsius - pressureMeltingPoint;
        % Populate the binary matrix
        aboveBelowMelt(aboveBelowMelt > -1) = 1; %1 is melted
        aboveBelowMelt(aboveBelowMelt < -1) = 0; %0 is frozen
        compileData = compileData + aboveBelowMelt;
            %{
            figure(i)
            imagesc(aboveBelowMelt);
            colorbar;
            caxis([-3,3]);
            title(strrep(extractAfter(litempbotgr{i},'litempbotgr'),'_',' '));
            xlabel('meters');
            ylabel('meters');
            %}
     end
end

x_coords = ncread("litempbotgr_AIS_DOE_MALI_ctrl_proj_std.nc","x");
y_coords = ncread("litempbotgr_AIS_DOE_MALI_ctrl_proj_std.nc","y");


figure(1)
    map = [ 0 0 180/255;...             %darkest blue -9  (9/9 frozen)
            59/255 123/255 223/255;...  %darker blue -7   (7/9 frozen)
            102/255 178/255 245/255;... %blue -5          (5/9 frozen)
            153/255 204/255 255/255;... %lighter blue -3  (3/9 frozen)
            204/255 229/255 255/255;... %lightest blue -1 (1/9 frozen)
            255/255 202/255 202/255;... %lightest red 1   (
            255/255 153/255 153/255;... %lighter red 3
            245/255 102/255 102/255;... %red 5
            230/255 70/255 70/255;...   %darker red 7
            187/255 33/255 33/255];     %darkest red 9
   datamin = 0;
   datamax = 9;
   


   % Plot bare rock and ice shelves around the figure
   filename = "BedMachineAntarctica_2020-07-15_v02.nc";
   [mask, thickness] = readGeometry(filename);
   mask(~isnan(compileData)) = NaN;

   mask_interp = interp2(BM3_x, BM3_y, mask, x_target, y_target, 'linear');
    
    
   colorshelf = [1 1 1]/1.5; %grey color in rgb 
   colorrock = [1 1 1]/2; %dark grey color in rgb
   white = [1 1 1]; %white in rgb
    
   %Create RGB matrix for the dataset
   image_rgb = ind2rgb(uint16((compileData - datamin)*(length(map)/(datamax-datamin))),map);
        
   %Find NaN and change to white
   [posi posj] = find(isnan(compileData));
   pos = sub2ind(size(image_rgb),repmat(posi,1,3),repmat(posj,1,3),repmat(1:3,size(posi,1),1));
   image_rgb(pos) = ones(size(pos,1),1)*white;  
    
   %Find position of ice shelves and change to grey
   %[posi posj] = find(mask == 3);
   %pos = sub2ind(size(image_rgb),repmat(posi,1,3),repmat(posj,1,3),repmat(1:3,size(posi,1),1));
   %image_rgb(pos) = ones(size(pos,1),1)*colorshelf; 

    ice_shelf_mask = (mask_interp == 3);
    for c = 1:3
        channel = image_rgb(:, :, c);
        channel(ice_shelf_mask) = colorshelf(c);
        image_rgb(:, :, c) = channel;
    end
    
    %Find position of rock and change to dark grey
    [posi posj] = find(mask_interp == 1);
    pos = sub2ind(size(image_rgb),repmat(posi,1,3),repmat(posj,1,3),repmat(1:3,size(posi,1),1));
    image_rgb(pos) = ones(size(pos,1),1)*colorrock; 
    
    %image_rgb(mask == 3) = 0;
    imagesc(x_8km .* 1e3, y_8km .* 1e3, image_rgb);
    caxis([0 10]);
    colorbar('Ticks',[0.5 1.5 2.5 3.5 4.5 5.5 6.5 7.5 8.5 9.5],... 
'TickLabels',["9/9 Frozen" "8/9 Frozen" "7/9 Frozen" "6/9 Frozen" "5/9 Frozen",...
                "5/9 Thawed" "6/9 Thawed" "7/9 Thawed" "8/9 Thawed" "9/9 Thawed"],...
                'TickLength',0, 'FontSize', 12.5);
    colormap(map);
    colorbar.Position = [0.9 0.1 0.02 0.8]; % Adjust these values as needed
    hold on
    for ind = 1:length(basal_therm)
        if basal_therm(ind) == 1
            scatter(x_coords_verif(ind), y_coords_verif(ind) * -1, 150, ... % Increased size
                'filled', 'MarkerEdgeColor', [0 0 0], ...
                'MarkerFaceColor', [1 0 0], 'LineWidth', 2); % Increased outline width
        elseif basal_therm(ind) == 0
                 scatter(x_coords_verif(ind), y_coords_verif(ind) * -1, 150, ... % Increased size
                'filled', 'MarkerEdgeColor', [0 0 0], ...
                'MarkerFaceColor', [30/255, 144/255, 255/255], 'LineWidth', 2); % Dodger Blue color and increased outline width
        end
    end
    %title('Agreement Between Nine ISMIP6 AIS Models');
    %xlabel('meters');
    %ylabel('meters');
        x0 = -2.7E6;
    y0 = 2.1E6;
    lengthscale = 1E6;
    widthscale = 50000;
    fontcolor = 'black';
    %Draw the scale bar
    p = patch([x0, x0 + lengthscale, x0 + lengthscale, x0], ...
          [y0, y0, y0 + widthscale, y0 + widthscale], fontcolor, ...
          'EdgeColor', fontcolor);
    %Add the scale bar label
    text(x0 + lengthscale/2, 1.90E6, sprintf('%.1f km', lengthscale/1E3), ...
        'HorizontalAlignment', 'center', 'Color', fontcolor, 'FontSize', 14);
    axis equal
    axis off
    set(gcf,'color','w')
    colorbar.Position = [0.87 0.2 0.02 0.6]; % Adjust these values as needed


    error
    
    
    %% Create a binary image of likely melted / frozen

binaryMap = [0 0 190/255; %blue
              190/255 0 0; %red
              1/1.5 1/1.5 1/1.5;]; %gray

allAgreedFrozen = sum(compileData(:) == 0); % Count occurrences of -9
allAgreedMelted = sum(compileData(:) == 9);  % Count occurrences of 9
totalNumPixelsIce = sum(~isnan(compileData(:))); % Count non-NaN elements
totalAgreement = ((allAgreedFrozen + allAgreedMelted)/totalNumPixelsIce) * 100
agreedMostlyFrozen = sum(compileData(:) == 2 | compileData(:) == 1);
agreedMostlyMelted = sum(compileData(:) == 7 | compileData(:) == 8);
agreedMostly = ((agreedMostlyFrozen + agreedMostlyMelted)/totalNumPixelsIce)*100
lowAgreedFrozen = sum(compileData(:) == 5| compileData(:) == 6); % Count occurrences of -9
lowAgreedMelted = sum(compileData(:) == 4|compileData(:) == 3);  % Count occurrences of 9
lowAgreed = ((lowAgreedFrozen + lowAgreedMelted)/totalNumPixelsIce)*100

smallerSize = [761, 761];
compileData(compileData >= 0) = 1; %melted
compileData(compileData < 0) = 0; %frozen
binaryImage = downsampleGrid_gpt(compileData, smallerSize);
%binaryImage legend
%NaN = ocean
%0 = frozen
%1 = melted

ismip6_8km = binaryImage;


subglacial = load("binaryRouting.mat").binaryLakesMethod;
%x_range = linspace(-3.333E6, 3.333E6, 761);
%y_range = linspace(-3.333E6, 3.333E6, 761);

%Create meshgrid for the original ismip6 data
%[X_ismip6, Y_ismip6] = meshgrid(x_coords, y_coords);

%Create meshgrid for the target grid (which is the same as subglacial map)
%[X_target, Y_target] = meshgrid(x_range, y_range);

%Interpolate ismip6 data onto the target grid
%ismip6_8km_interp = interp2(X_ismip6, Y_ismip6, ismip6_8km, X_target, Y_target, 'linear');
%lake_vostok = load("lake_vostok_mask.mat").lake_vostok_mask;
%ismip6_8km_interp(lake_vostok == 1) = 0;

pathname = fileparts('/Users/owenseiner/Desktop/Academics/Past_Terms/Antarctica_Research/Antarctica');
matfile = fullfile(pathname, 'ismip6_8km.mat');
save(matfile, "ismip6_8km")

ismip6_binary_8km = load("ismip6_8km.mat").ismip6_8km;
%ismip6_binary_8km(isnan(ismip6_binary_8km)) = NaN;

BM3_mask_ais = load("BM3_mask_ais.mat").BM3_mask_ais;
ismip6_binary_8km(BM3_mask_ais == 2) = 3;
    
    x_8km = -3040:8:3040; %in km
    y_8km = -3040:8:3040; %in km
    
    figure
    h = imagesc(x_8km * 1e3,y_8km * 1e3,ismip6_binary_8km);
    set(h, 'AlphaData', ~isnan(ismip6_binary_8km));  % Make NaNs transparent
    hold on
    h1 = scatter(x_coords_verif(basal_therm == 1), y_coords_verif(basal_therm == 1) * -1, 150, 'filled', 'MarkerEdgeColor', [0 0 0], 'MarkerFaceColor', [1 0 0], 'LineWidth', 2); % Red points
    h2 = scatter(x_coords_verif(basal_therm == 0), y_coords_verif(basal_therm == 0) * -1, 150, 'filled', 'MarkerEdgeColor', [0 0 0], 'MarkerFaceColor', [30/255, 144/255, 255/255], 'LineWidth', 2); % Dodger blue points
    h3 = plot(nan, nan, 's', 'MarkerEdgeColor', 'none', 'MarkerFaceColor', binaryMap(1,:), 'MarkerSize', 15); % Blue square for legend
    h4 = plot(nan, nan, 's', 'MarkerEdgeColor', 'none', 'MarkerFaceColor', binaryMap(2,:), 'MarkerSize', 15); % Red square for legend
    h5 = plot(nan, nan, 's', 'MarkerEdgeColor', 'none', 'MarkerFaceColor', binaryMap(3,:), 'MarkerSize', 15); % Gray square for legend
    colormap(binaryMap)
    x0 = -2.7E6;
    y0 = 2.1E6;
    lengthscale = 1E6;
    widthscale = 50000;
    fontcolor = 'black';
    %Draw the scale bar
    p = patch([x0, x0 + lengthscale, x0 + lengthscale, x0], ...
          [y0, y0, y0 + widthscale, y0 + widthscale], fontcolor, ...
          'EdgeColor', fontcolor);
    %Add the scale bar label
    text(x0 + lengthscale/2, 1.90E6, sprintf('%.1f km', lengthscale/1E3), ...
        'HorizontalAlignment', 'center', 'Color', fontcolor, 'FontSize', 14);
    axis equal
    lgd = legend([h1, h2, h3, h4, h5], {'Thawed Boreholes', 'Frozen Boreholes', 'Likely Frozen', 'Likely Thawed', 'Ice Shelves'}, 'Location', 'northeastoutside', 'FontSize', 14, 'FontName', 'Sans Serif');
    lgd.Position = [0.65, 0.7, 0.2, 0.2];
    lgd.Box = 'off';
    set(gcf,'color','w');
    axis off
    set(gcf,'color','w');
    



   
    