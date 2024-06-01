%Import the active lakes from Livingstone et al and convert to table
lakesData2 = load('ActiveSubglacialLivingstone.mat');
%lakesData2 = readmatrix("Seven_Active_Lakes.xlsx");
lakesData2 = lakesData2.ActiveSubglacialLakesOnlyLivingstone;

%Import the stable lakes from Livingstone et al and convert to table
lakes_stable = load('StableSubglacialLivingstone.mat');
lakes_stable = lakes_stable.StableSubglacialLivingstone;

%BedMachine v2 Antarctica, needed for water routing using laptop
BM3                         = struct;
BM3.x                       = 1e-3 .* double(ncread('/Users/owenseiner/Desktop/Academics/Past_Terms/Antarctica_Research/Antarctica/BedMachineAntarctica-v3.nc', 'x')); % projected x, km
BM3.y                       = 1e-3 .* double(flipud(ncread('/Users/owenseiner/Desktop/Academics/Past_Terms/Antarctica_Research/Antarctica/BedMachineAntarctica-v3.nc', 'y'))); % projected y, km
BM3.elev_bed                = rot90(ncread('/Users/owenseiner/Desktop/Academics/Past_Terms/Antarctica_Research/Antarctica/BedMachineAntarctica-v3.nc', 'bed')); % bed elevation, m
BM3.elev_bed_uncert         = rot90(ncread('/Users/owenseiner/Desktop/Academics/Past_Terms/Antarctica_Research/Antarctica/BedMachineAntarctica-v3.nc', 'errbed')); % bed elevation uncertainty, m
BM3.elev_surf               = double(rot90(ncread('/Users/owenseiner/Desktop/Academics/Past_Terms/Antarctica_Research/Antarctica/BedMachineAntarctica-v3.nc', 'surface'))); % surface elevation, m
BM3.mask_ais               = double(rot90(ncread('/Users/owenseiner/Desktop/Academics/Past_Terms/Antarctica_Research/Antarctica/BedMachineAntarctica-v3.nc', 'mask'))); % ice mask
BM3.thick                   = double(rot90(ncread('/Users/owenseiner/Desktop/Academics/Past_Terms/Antarctica_Research/Antarctica/BedMachineAntarctica-v3.nc', 'thickness'))); % ice thickness, m
BM3.thick(~BM3.thick)       = NaN; % NaN out zero ice thickness

%Determine the size of the lakes tables
[row_count_active, col_count_active] = size(lakesData2);
[row_count_stable, col_count_stable] = size(lakes_stable);

%Create tables to store x and y (in m) of lakes on a 13333x13333 grid
x_excel_active = zeros(row_count_active, 1);
y_excel_active = zeros(row_count_active,1);
x_excel_stable = zeros(row_count_stable, 1);
y_excel_stable = zeros(row_count_stable,1);

%Convert lakes tables to arrays with only lat and lon
lakesTable2 = table2array(lakesData2(1:172, 3:4));
%lakesTable2 = lakesData2(1:8, 3:4);
lakesTableStable = table2array(lakes_stable(1:503, 3:4));

%Convert from lat and lon to x and y
%Use -1 to indicate southern hemisphere
%y is latitude (1st column), x is longitude (2nd column)
[x_excel_active, y_excel_active] = ll2xy(lakesTable2(:,1), lakesTable2(:, 2), -1);
[x_excel_stable, y_excel_stable] = ll2xy(lakesTableStable(:,1), lakesTableStable(:, 2), -1);


%Begin topo toolbox analysis
%Create a DEM of the bathymetry using GRIDobj
%Head is used instead of elev_bed to better adjust for pressure from ice
%addpath([issmdir() '/externalpackages/topotoolbox/install']);
head = BM3.elev_bed + BM3.thick.*(917/1000);
[X Y]=meshgrid(BM3.x.*1000, BM3.y.*1000);
DEM = GRIDobj(X, Y, head);
lakesGrid = zeros(13333, 13333);

%Need to convert from x and y to indices on a 13333x13333 grid
%Determine the indices from x and y coordinates using an anonymous function
BM3_coord_to_index = @(coord) (2*10^-3*coord + 6667);

%Create arrays to store the indices
index_y_active = zeros(172,1);
index_x_active = zeros(172, 1);
%index_y_active = zeros(8,1);
%index_x_active = zeros(8, 1);
index_y_stable = zeros(503, 1);
index_x_stable = zeros(503,1);

%y references latitude --> row on lakesGrid
%x references longitude --> column on lakesGrid
index_x_active = round(BM3_coord_to_index(x_excel_active));
index_y_active = round(BM3_coord_to_index(y_excel_active));
for i=1:length(index_x_active)
	lakesGrid(index_y_active(i), index_x_active(i)) = 1;
end

index_x_stable = round(BM3_coord_to_index(x_excel_stable));
index_y_stable = round(BM3_coord_to_index(y_excel_stable));
  
%Create a DEM of the bathymetry using GRIDobj
%head is used instead of elev_bed to better adjust for pressure from ice
head = BM3.elev_bed + BM3.thick.*(917/1000);
%Dimensionless gradient of hydraulic head (Darcy slope)
[dHdx, dHdy] = gradient(head);

gradient_magnitude = sqrt(dHdx.^2 + dHdy.^2);

p75 = prctile(gradient_magnitude(:), 75);
%quiver(X, cos(flow_direction), sin(flow_direction), 'r');
%Visual water flow direction with regards to x axis
flow_direction = atan2(dHdy, dHdx);
head = BM3.elev_bed + BM3.thick.*(917/1000);
DEM = GRIDobj(X, Y, head);
%lakesGrid = readmatrix("lakesGrid.csv");


%Edit the lakesGrid to have 1s where subglacial 
%disp('Filling sinks...')
%DEMf = fillsinks(DEM);
%disp('Generating flow obj...')
%FD = FLOWobj(DEMf);
%disp('Generating basal melt...')
%basal = lakesGrid * 300;
%basalmelt = GRIDobj(X,Y, basal);
%disp('Calculating accumulation...')
%Accumw = flowacc(FD,basalmelt);
%disp('Generating stream object...')
%Aasign minimum upstream area to value of 1 pixel
%S = STREAMobj(FD, 'minarea', 1);
Accumw = load('accum100.mat').ACC;
stream = load('streamobj.mat').S;


%Convert back to 8km ISMIP6 grid
lakes_8km = zeros(761,761);
x_8km = -3040:8:3040; %in km
y_8km = -3040:8:3040; %in km

%ISMIP6 index from bedmachine index (valid for both x and y)
ISMIP6_index= @(BM3_index) round(1/16*(BM3_index -1) -285/8);

%All ice to frozen (0) then shelves and lakes changes
[indy_sheet_bm3 indx_sheet_bm3] = find(BM3.mask_ais==2);
indx_sheet_ismip = ISMIP6_index(indx_sheet_bm3);
indy_sheet_ismip = ISMIP6_index(indy_sheet_bm3);

for i = 1:length(indx_sheet_ismip)
	lakes_8km(indy_sheet_ismip(i),indx_sheet_ismip(i))=1;
end

Accumw.Z=flipud(Accumw.Z); %need to flip to stay consisten
[indy_bm3 indx_bm3] = find(Accumw.Z>0);

indx_ismip = ISMIP6_index(indx_bm3);
indy_ismip = ISMIP6_index(indy_bm3);

for i = 1:length(indx_ismip)
	lakes_8km(indy_ismip(i),indx_ismip(i))=2;
end

%Active lakes
indx_active_ismip = ISMIP6_index(index_x_active);
indy_active_ismip = ISMIP6_index(index_y_active);

for i = 1:length(indx_active_ismip)
	lakes_8km(indy_active_ismip(i),indx_active_ismip(i))=3;
end

%Stable lakes + Vostok
indx_stable_ismip = ISMIP6_index(index_x_stable);
indy_stable_ismip = ISMIP6_index(index_y_stable);

for i = 1:length(indx_stable_ismip)
	lakes_8km(indy_stable_ismip(i),indx_stable_ismip(i))=4;
end

[indy_vostok_bm3 indx_vostok_bm3] = find(BM3.mask_ais==4);

indx_vostok_ismip = ISMIP6_index(indx_vostok_bm3);
indy_vostok_ismip = ISMIP6_index(indy_vostok_bm3);

save("indx_vostok.mat", "indx_vostok_ismip")
save("indy_vostok.mat", "indy_vostok_ismip")

for i = 1:length(indx_vostok_ismip)
	lakes_8km(indy_vostok_ismip(i),indx_vostok_ismip(i))=4;
end

%Ice shelves
[indy_shelves_bm3 indx_shelves_bm3] = find(BM3.mask_ais==3);
indx_shelves_ismip = ISMIP6_index(indx_shelves_bm3);
indy_shelves_ismip = ISMIP6_index(indy_shelves_bm3);

save("indx_shelves.mat", "indx_shelves_ismip")
save("indy_shelves.mat", "indy_shelves_ismip")

for i = 1:length(indx_shelves_ismip)
	lakes_8km(indy_shelves_ismip(i),indx_shelves_ismip(i))=5;
end

%NaN ocean/no ice
%1 frozen
%2 routing
%3 active lakes
%4 stable lakes + vostok
%5 ice shelves


lakesrouting_8km = flipud(load("lakesrouting8km.mat").lakes_8km);
binaryLakesRouting = zeros(size(lakesrouting_8km));
binaryLakesRouting(lakesrouting_8km == 2 | lakesrouting_8km == 3 | lakesrouting_8km == 4) = 1;
binaryLakesRouting(lakesrouting_8km == 1) = 0;
binaryLakesRouting(lakesrouting_8km == 0) = NaN;

save("binaryLakesRouting.mat", "binaryLakesRouting");



figure(3); imagesc(x_8km, y_8km, lakes_8km); set(gca,'Ydir','Normal'); caxis([0 5]);

error
%[x,y] = STREAMobj2XY(stream);


%x = load('streamx.mat').S_x;
%y = load('streamy.mat').S_y;
%ixc_stream = stream.ixc;
%ix_stream = stream.ix; 
%save("ixc_stream.mat","ixc_stream");
%save("ix_stream.mat", "ix_stream");

%Create a matrix of accumulation from accum_matrix
%accum_matrix = Accumw.Z;
%stream = accum_matrix > 0;


%Define the size of the larger and smaller grids
smallerSize = [761, 761];

%stream_8km = downsampleGrid(stream, smallerSize);
%[y_coordinates, x_coordinates] = find(stream_8km);



% Define the grid size
gridSize = 13333;

%Scale and shift the coordinates to fit within the grid
x_scaled = round((x + 3.333E6) / 6E6 * (gridSize - 1)) + 1;
y_scaled = round((y + 3.333E6) / 6E6 * (gridSize - 1)) + 1;

%Ensure that indices are within the valid range
x_scaled = max(1, min(x_scaled, gridSize));
y_scaled = max(1, min(y_scaled, gridSize));


%Initialize a logical grid with zeros
stream_mask = false(gridSize);

% Convert subscripts to linear indices
indices = sub2ind(size(stream_mask), y_scaled, x_scaled);

% Set specific indices to 1
stream_mask(indices) = true;

stream_mask = flipud(stream_mask);


%Assuming S.ixc and S.ix represent connectivity information
%Convert linear indices to row and column indices
[row, col] = ind2sub([gridSize, gridSize], [stream.ixc; stream.ix]);

% Set the values in the sparse matrix where the stream network exists to true
stream_mask(sub2ind([gridSize, gridSize], row, col)) = true;

%downsampledStreamobj = downsampledGrid_sum(stream_mask, smallerSize);



%stream_reduced = load('streamreduced.mat');
%stream_X = load('streamx.mat');
%stream_Y = load('streamy.mat');


%FD = load('flowobj.mat').FD;
%D = load('flowdistance.mat').D;

%addpath('/Users/owenseiner/Desktop/Academics/Past_Terms/Antarctica_Research/Antarctica/topotoolbox/@STREAMobj')


%Create a matrix of accumulation from accum_matrix
accum_matrix = Accumw.Z;
%stream = accum_matrix > 0;


%Define the size of the larger and smaller grids
smallerSize = [761, 761];

downsampledBM3Mask = flipud(downsampleClassificationGrid(BM3.mask_ais, smallerSize));
save("BM3Mask_AIS_8km.mat", "downsampledBM3Mask");
%Downsample the gradient of magnitude (unitless)
gradient_8km = flipud(downsampleGrid(gradient_magnitude, smallerSize));
save("gradient_8km.mat", "gradient_8km")


%Call the downsampling function
downsampledStream = downsampleLakesGrid(stream_mask, smallerSize);


%Initialize stream_8km with zeros
stream_8km = zeros(size(downsampledStream));

%0 = ocean and ice free land
%1 = subglacial water route
%2 = floating ice and lake vostok
%3 = active lake
%4 = stable lake
%5 = presumed frozen ice (test with hydraulic head transfer momentum threshold)
%6 = beyond hydraulic head gradient threshold

%Set everywhere in stream_8km where downsampledGrid > 0 to be 1
stream_8km(downsampledStream > 0) = 1;

%Set everywhere where downsampledBM3Mask is 3 (floating ice) or 4 (Lake Vostok) to be 1 in stream_8km
lake_vostok_mask = (downsampledBM3Mask == 4);
%stream_8km(lake_vostok_mask) = 1;

save("lake_vostok_mask.mat", "lake_vostok_mask")

%save("lake_vostok_mask.mat", "lake_vostok_mask");
stream_8km(downsampledBM3Mask == 3) = 2;

%Set everywhere where BM3mask is 2 == 5 & not stream (presumed frozen bed)
frozen_ice_mask = (downsampledBM3Mask == 2) & (stream_8km ~= 1);

%Update the values in stream_8km based on the mask
stream_8km(frozen_ice_mask) = 5;

%Set where downsampledBM3Mask is 0 (ocean) or 1 (ice free land) to be 0 in stream_8km
stream_8km(downsampledBM3Mask == 0 | downsampledBM3Mask == 1 ) = 0;

hydraulic_head_mask = ((stream_8km == 5)  & (gradient_8km > p75));
stream_8km(hydraulic_head_mask) = 6;

%Separating lakesGridActive and stable works and produces correct num lakes
%The problem right now is the projection is slightly off
lakesGridActive = zeros(13333, 13333);
lakesGridStable = zeros(13333, 13333);
for n = 1:height(index_x_active)
    lakesGridActive(index_x_active(n,1), index_y_active(n,1)) = 1;
end

activeLake_8km = zeros(smallerSize);
activeLake_8km = downsampleLakesGrid(lakesGridActive, [761, 761]);

for n = 1:height(index_x_stable)
    lakesGridStable(index_x_stable(n,1), index_y_stable(n,1)) = 1;
end

stableLake_8km = zeros(smallerSize);
stableLake_8km = downsampleLakesGrid(lakesGridStable, [761,761]);

%Reverse the row on the subglacial lakes
stableLake_8km = flipud(stableLake_8km);
activeLake_8km = flipud(activeLake_8km);


%Set where lakesGrid_8km is 1 to be 1 in stream_8km
stream_8km(activeLake_8km == 1) = 3;
stream_8km(stableLake_8km == 1) = 4;

%Generate x and y for plotting
x_range = linspace(-3.333E6, 3.333E6, 761);
y_range = linspace(-3.333E6, 3.333E6, 761);

subglacial_binary_8km = load("binaryLakesRouting.mat").binaryLakesRouting;

%Set colormap based on your legend categories
customColormap = [
    1.0 1.0 1.0;   % 0: Ocean and Ice Free Land (White)
    0.2 0.8 0.2;   % 1: Subglacial Water Route (Green)
    1.0 1.0 1.0;   % 2: Floating Ice and Lake Vostok (White)
    1.0 0.5 0.0;   % 3: Active Lake (Orange)
    0.5 0.5 1.0;   % 4: Stable Lake (Light Blue)
    0.8 0.8 0.8;   % 5: Presumed Frozen Base (Gray)
    0.7 0.7 0.2    % 6: Hydraulic Head Threshold (yellowish-green)
    0.5 0.5 0.5    % 7: Ice shelves
];



%legendLabels = {'Ocean/Ice Free Land', 'Subglacial Water Route', 'Floating Ice/Lake Vostok', ...
 %   'Active Lake', 'Stable Lake', 'Presumed Frozen Base'};

stream_8km(downsampledBM3Mask == 1 | downsampledBM3Mask == 3) = 7;

figure
h = imagesc(x_8km * 1e3, y_8km * 1e3, stream_8km);
colormap(customColormap);


hold on

% Find the indices where the values are 3 or 4
indices3 = (stream_8km == 3);
indices4 = (stream_8km == 4);

km8_index_to_coord = @(index) index * (2 * 3.333E6) / (762 - 1) - 3.333E6;


%Create a scatter plot to visually emphasize the values
[row3, col3] = find(indices3);
scatter(km8_index_to_coord(col3), km8_index_to_coord(row3), 'MarkerEdgeColor', [1.0 0.5 0.0], 'MarkerFaceColor', [1.0 0.5 0.0], 'SizeData', 10);

[row4, col4] = find(indices4);
scatter(km8_index_to_coord(col4), km8_index_to_coord(row4), 'MarkerEdgeColor', [0.5 0.5 1.0], 'MarkerFaceColor', [0.5 0.5 1.0], 'SizeData', 10);
axis equal
xlabel("X (m)")
ylabel("Y (m)")
%title("Categorized Subglacial Lakes Routing & Hydraulic Gradient Threshold")
x0 = -2.7E6;
y0 = 2.1E6;
lengthscale = 1E6;
widthscale = 50;
fontcolor = 'black';

%Draw the scale bar
p = patch([x0, x0 + lengthscale, x0 + lengthscale, x0], ...
          [y0, y0, y0 + widthscale, y0 + widthscale], fontcolor, ...
          'EdgeColor', fontcolor);
%Add the scale bar label
text(x0 + lengthscale/2, 1.95E6, sprintf('%.1f km', lengthscale/1E3), ...
     'HorizontalAlignment', 'center', 'Color', fontcolor);
set(gcf,'color','w')
axis off
%0 = ocean and ice free land
%1 = subglacial water route
%2 = floating ice and lake vostok
%3 = active lake
%4 = stable lake
%5 = presumed frozen base 
%6 = hydraulic head gradient threshold

%Create a binary melted/frozen mask
%Create a matrix of zeros same size as stream_8km
%0 = ocean/ice free land/floating ice and lake vostok
%1 = melted base (lake, water route, gradient threshold)
%2 = frozen base 
%3 = ice free land
binaryLakesMethod = zeros(761,761);
binaryLakesMethod((stream_8km == 0) | (stream_8km == 2)) = 0;
indices = ((stream_8km == 1) | (stream_8km == 3) | (stream_8km == 4) | (stream_8km == 6)) & (downsampledBM3Mask == 2);
binaryLakesMethod(indices) = 1;
binaryLakesMethod(stream_8km == 5) = 2;
binaryLakesMethod(downsampledBM3Mask == 1 | downsampledBM3Mask == 3) = 3;
%binaryLakesMethod(downsampledBM3Mask == 4) = 1;

%Define the colormap
customColormap = [1 1 1; 1 0 0; 0 0 1; 0.5 0.5 0.5]; % White, Red, Blue, Grey

figure
imagesc(x_8km * 1e3, y_8km * 1e3, binaryLakesMethod)
colormap(customColormap)
title("Binary Frozen/Melted Map From Routing Method")
axis equal
axis off
%Add the scale bar label
% Add scale bar
x0 = -2.7E6;
y0 = 2.1E6;
lengthscale = 1E6;
widthscale = 50;
fontcolor = 'black';

%Draw the scale bar
p = patch([x0, x0 + lengthscale, x0 + lengthscale, x0], ...
          [y0, y0, y0 + widthscale, y0 + widthscale], fontcolor, ...
          'EdgeColor', fontcolor);
%Add the scale bar label
text(x0 + lengthscale/2, 1.95E6, sprintf('%.1f km', lengthscale/1E3), ...
     'HorizontalAlignment', 'center', 'Color', fontcolor);
set(gcf,'color','w')

ice_shelf_mask = zeros(size(BM3.mask_ais));
ice_shelf_mask(flipud(BM3.mask_ais == 3 | BM3.mask_ais == 1)) = 1; 
ice_shelf_8km = downsampleBinaryGrid_gpt(ice_shelf_mask, smallerSize);
save("binaryRouting.mat", "binaryLakesMethod");
save("ice_shelf_8km.mat", "ice_shelf_8km");







