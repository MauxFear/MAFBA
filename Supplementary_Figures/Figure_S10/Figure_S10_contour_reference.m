% Directory path
directory_path = 'D:/Code/MATLAB/Projects/MAFBA/triData/';
fileVersion = 'uProtZ_150823.mat';
% Load the data from the MAT file
fileName = sprintf('%s%s',directory_path, fileVersion);
load(fileName);
phi_max = 0.484;
triArray_Z(:,end+1) = triArray_Z(:,2) / phi_max;
% Create table to manipulate the data
z_df = array2table(triArray_Z, 'VariableNames', {'gamma', 'phi_Z', 't_ac', 'Z_divided_by_max'}  );
z_df = sortrows(z_df, {'gamma','Z_divided_by_max'});
disp(z_df(1:10,:))
 
fileVersion = 'uProtY_150823.mat';
% Load the data from the MAT file
fileName = sprintf('%s%s',directory_path, fileVersion);
load(fileName);
phi_max = 0.484;
triArray_Y(:,end+1) = triArray_Y(:,2) / phi_max; % Norm with phi max
triArray_Y(:,end+1) = triArray_Y(:,2) ./ (triArray_Y(:,1) .* (phi_max)); % Norm with phi M max
% Create table to manipulate the data
y_df = array2table(triArray_Y, 'VariableNames', {'gamma', 'phi_Y', 't_ac', 'Y_divided_by_max', 'Y_divided_by_Mmax'}  );
y_df = sortrows(y_df, {'gamma','Y_divided_by_max'});

disp(y_df(1:10,:))% Initialize an empty table to store the data


%%

fig = figure('Units', 'centimeters', 'Position', [1, 1, 13, 12.5]);
hold on;

% markers = {'o', '^', 's', 'd', '>', 'v', '<'}; %Y %%%% CHANGE THIS %%%%%
markers = {'-o', '-^', '-s', '-d', '->', '-v', '-<'}; %Y %%%% CHANGE THIS %%%%%
colors = {
        sscanf('C00217','%2x%2x%2x',[1 3])/255,... % dark red
        sscanf('6B8E23','%2x%2x%2x',[1 3])/255,... % Olive Green
        sscanf('008080','%2x%2x%2x',[1 3])/255,... % Teal
        sscanf('40E0D0','%2x%2x%2x',[1 3])/255,... % Turquoise
        sscanf('87CEEB','%2x%2x%2x',[1 3])/255,... % Sky Blue
        sscanf('6495ED','%2x%2x%2x',[1 3])/255,... % Cornflower Blue
        sscanf('1F56A8','%2x%2x%2x',[1 3])/255,... % blue
        [0.07450980392156863, 0.06666666666666667, 0.0196078431372549],... %black
        sscanf('FFFFFF','%2x%2x%2x',[1 3])/255,... %white
             };

% title(gamma_string)
% xlabel('\phi Z /\phi_{max} or \phi Y /\phi^{m}_{max}', 'FontSize', 14)
xlabel('\phi Z /\phi_{max} or \phi Y /\phi_{max}', 'FontSize', 14)

ylabel('Acetate Threshold (\lambda_{ac}) (h^{-1})', 'FontSize', 14)
set(gca, 'FontSize', 12)


%Plot Membrane proteins
ygammas = unique(y_df.gamma);
gamma_select =[1,0.35,0.3,0.25,0.23,0.2] ;
gamma_select= gamma_select(end:-1:1);

for i=1:numel(gamma_select)
    gammaVal = gamma_select(i);
    sub_y_df = y_df(y_df.gamma== gammaVal,:); 
    sub_y_df = sortrows(sub_y_df, {'gamma','Y_divided_by_max'});
    indexZero = find(sub_y_df.t_ac ==0,1);
    sub_y_df = sub_y_df(1:indexZero, :);
        
    plot(sub_y_df.Y_divided_by_max, sub_y_df.t_ac, '-', 'Color', colors{1+i}, 'LineWidth', 4, 'HandleVisibility', 'off');
    % Downsample the data to only keep 20 points with markers
    downsampleFactor = 3;  % Keep every 5th point
    numPoints = numel(sub_y_df.t_ac);
    indicesToShow = 1 : downsampleFactor : numPoints;
       
    if gammaVal == 1
        indicesToShow = sort ([indicesToShow, 54,55]);
  
    end
    
    gammaLabel = sprintf('Membrane (\\phi Y) \\gamma %s', num2str(gammaVal)); 
    plot(sub_y_df.Y_divided_by_max(indicesToShow), sub_y_df.t_ac(indicesToShow),...
    markers{1+i}, 'Color', colors{1+i}, 'LineWidth', 1.5,'MarkerEdgeColor', colors{1+i},...
    'MarkerFaceColor', colors{end}, 'MarkerSize', 8, 'DisplayName', gammaLabel);
end

% Plot Cytosolic Protein
gamma_select =[1,0.25,0.23,0.2] ;
for i=1:numel(gamma_select)
    gammaVal = gamma_select(i);
    sub_z_df = z_df(z_df.gamma== gammaVal,:); 
    plot(sub_z_df.Z_divided_by_max, sub_z_df.t_ac, '-', 'Color', colors{1}, 'LineWidth', 4, 'HandleVisibility', 'off');
    % Downsample the data to only keep 20 points with markers
    downsampleFactor = 4;  % Keep every 5th point
    numPoints = numel(sub_z_df.t_ac);
    indicesToShow = 1 : downsampleFactor : numPoints;
    indicesToShow = sort ([indicesToShow, 52,51]);
    h1=plot(sub_z_df.Z_divided_by_max(indicesToShow), sub_z_df.t_ac(indicesToShow),...
        markers{1}, 'Color', colors{1}, 'LineWidth', 1.5,'MarkerEdgeColor', colors{1},...
        'MarkerFaceColor', colors{end}, 'MarkerSize', 8, 'DisplayName', 'Cytosolic (\phi Z) ');
end

xlim([0, 1]);
y = ylim;
% Set the legend
lgd = legend('Location', 'southoutside','NumColumns', 2);
set(lgd, 'Box', 'off', 'FontSize', 12);
% Show the grid
grid on;

% Save the figure as SVG
filename = sprintf('uProtvsAceT_AlterV0_lbl.svg');
saveas(fig, filename, 'svg');


%%

fig = figure('Units', 'centimeters', 'Position', [1, 1, 13, 12.5]);
hold on;
markers = {'o', '^', 's', 'd', '>', 'v', '<'}; %Y %%%% CHANGE THIS %%%%%
colors = {
        sscanf('FFA500','%2x%2x%2x',[1 3])/255,... % Mango Orange 
        sscanf('FF6F61','%2x%2x%2x',[1 3])/255,... % Coral 
        sscanf('FF3E3E','%2x%2x%2x',[1 3])/255,... % Tomato Red
        sscanf('DC143C','%2x%2x%2x',[1 3])/255,... % Crimson
        sscanf('C00217','%2x%2x%2x',[1 3])/255,... % dark red
        [0.07450980392156863, 0.06666666666666667, 0.0196078431372549],... %black
        sscanf('FFFFFF','%2x%2x%2x',[1 3])/255,... %white
             };
         
xlabel('\phi Z /\phi_{max}', 'FontSize', 14)
ylabel('Acetate Threshold (\lambda_{ac}) (h^{-1})', 'FontSize', 14)
set(gca, 'FontSize', 12)

% Plot Cytosolic Protein
zgammas = unique(z_df.gamma);
gamma_select =[1,0.3,0.25,0.23,0.2] ;
gamma_select= gamma_select(end:-1:1);

for i=1:numel(gamma_select)
    gammaVal = gamma_select(i);
    sub_z_df = z_df(z_df.gamma== gammaVal,:); 
    sub_z_df = sortrows(sub_z_df, {'gamma','Z_divided_by_max'});
    indexZero = find(sub_z_df.t_ac ==0,1);
    sub_z_df = sub_z_df(1:indexZero, :);
    colorAdjusted = [colors{i}, 1];    
    plot(sub_z_df.Z_divided_by_max, sub_z_df.t_ac, '-', 'Color', colorAdjusted, 'LineWidth', 4, 'HandleVisibility', 'off');
    % Downsample the data to only keep 20 points with markers
    downsampleFactor = 4;  % Keep every 5th point
    numPoints = numel(sub_z_df.t_ac);
    indicesToShow = 1 : downsampleFactor : numPoints;
       
    if i <= 2
        indicesToShow = sort ([indicesToShow, 54,53]);
    else
        indicesToShow = sort ([indicesToShow, 52,51]);
    end
    
    gammaLabel = sprintf('Cytosolic (\\phi Z) \\gamma %s', num2str(gammaVal)); 
    plot(sub_z_df.Z_divided_by_max(indicesToShow), sub_z_df.t_ac(indicesToShow),...
    markers{i}, 'Color', colorAdjusted, 'LineWidth', 1.5,'MarkerEdgeColor', colors{i},...
    'MarkerFaceColor', colors{end}, 'MarkerSize', 8, 'DisplayName', gammaLabel);
end

xlim([0, 1]);
% Set the legend
lgd = legend('Location', 'southoutside','NumColumns', 2);
set(lgd, 'Box', 'off', 'FontSize', 12);
% Show the grid
grid on;

% Save the figure as SVG
filename = sprintf('uProtZvsAceT_AlterV0.svg');
saveas(fig, filename, 'svg');
%%

fig = figure('Units', 'centimeters', 'Position', [1, 1, 13, 12.5]);
hold on;

markers = {'o', '^', 's', 'd', '>', 'v', '<'}; %Y %%%% CHANGE THIS %%%%%
colors = {
        sscanf('6B8E23','%2x%2x%2x',[1 3])/255,... % Olive Green
        sscanf('008080','%2x%2x%2x',[1 3])/255,... % Teal
        sscanf('40E0D0','%2x%2x%2x',[1 3])/255,... % Turquoise
        sscanf('87CEEB','%2x%2x%2x',[1 3])/255,... % Sky Blue
        sscanf('6495ED','%2x%2x%2x',[1 3])/255,... % Cornflower Blue
        sscanf('1F56A8','%2x%2x%2x',[1 3])/255,... % blue
        [0.07450980392156863, 0.06666666666666667, 0.0196078431372549],... %black
        sscanf('FFFFFF','%2x%2x%2x',[1 3])/255,... %white
             };
         
xlabel('\phi Y /\phi_{max}', 'FontSize', 14)
ylabel('Acetate Threshold (\lambda_{ac}) (h^{-1})', 'FontSize', 14)
set(gca, 'FontSize', 12)

%Plot Membrane proteins
ygammas = unique(y_df.gamma);
gamma_select =[1,0.35,0.3,0.25,0.23,0.2] ;
gamma_select= gamma_select(end:-1:1);

for i=1:numel(gamma_select)
    gammaVal = gamma_select(i);
    sub_y_df = y_df(y_df.gamma== gammaVal,:); 
    sub_y_df = sortrows(sub_y_df, {'gamma','Y_divided_by_max'});
    indexZero = find(sub_y_df.t_ac ==0,1);
    sub_y_df = sub_y_df(1:indexZero, :);
        
    plot(sub_y_df.Y_divided_by_max, sub_y_df.t_ac, '-', 'Color', colors{i}, 'LineWidth', 4, 'HandleVisibility', 'off');
           
    if i == numel(gamma_select)
        downsampleFactor = 4;  % Keep every 5th point
        numPoints = numel(sub_y_df.t_ac);
        indicesToShow = 1 : downsampleFactor : numPoints;
        indicesToShow = sort ([indicesToShow, 54,55]);
    elseif i == numel(gamma_select)-1
        downsampleFactor = 4;  % Keep every 5th point
        numPoints = numel(sub_y_df.t_ac);
        indicesToShow = 1 : downsampleFactor : numPoints;
    else
        downsampleFactor = 3;  % Keep every 5th point
        numPoints = numel(sub_y_df.t_ac);
        indicesToShow = 1 : downsampleFactor : numPoints;
    end
    
    gammaLabel = sprintf('Membrane (\\phi Y) \\gamma %s', num2str(gammaVal)); 
    plot(sub_y_df.Y_divided_by_max(indicesToShow), sub_y_df.t_ac(indicesToShow),...
    markers{i+1}, 'Color', colors{i}, 'LineWidth', 1.5,'MarkerEdgeColor', colors{i},...
    'MarkerFaceColor', colors{end}, 'MarkerSize', 8, 'DisplayName', gammaLabel);
end

xlim([0, 1]);
y = ylim;
% Set the legend
lgd = legend('Location', 'southoutside','NumColumns', 2);
set(lgd, 'Box', 'off', 'FontSize', 12);
% Show the grid
grid on;

% Save the figure as SVG
filename = sprintf('uProtYvsAceT_AlterV0.svg');
saveas(fig, filename, 'svg');


%%
% ygammas = unique(y_df.gamma);
% y_normMax = unique(y_df.Y_divided_by_max);
% sub_y_df = y_df(y_df.gamma== 0.25,:);
% disp(sub_y_df);

%%
%Plotting Y Contour Plot
fig = figure();

% Define the colors for the colormap
white_color = [1 1 1];   % White color
blue_color =  sscanf('6495ED','%2x%2x%2x',[1 3])/255;  % Blue color
green_color = sscanf('6B8E23','%2x%2x%2x',[1 3])/255;  % Green color
aqua_color = sscanf('40E0D0','%2x%2x%2x',[1 3])/255;     % Aqua color
% Number of color points in the colormap
num_colors = 20;
% Create the custom colormap with a smooth gradient
Ycustom_color_map = [linspace(white_color(1), green_color(1), num_colors);
                    linspace(white_color(2), green_color(2), num_colors);
                    linspace(white_color(3), green_color(3), num_colors)].';
% Number of color points in the colormap
num_colors = 20;
% Append the red part of the colormap
Ycustom_color_map = [Ycustom_color_map; [linspace(green_color(1), aqua_color(1), num_colors);
                                       linspace(green_color(2), aqua_color(2), num_colors);
                                       linspace(green_color(3), aqua_color(3), num_colors)].'];
% Number of color points in the colormap
num_colors = 20;
% Append the red part of the colormap
Ycustom_color_map = [Ycustom_color_map; [linspace(aqua_color(1), blue_color(1), num_colors);
                                       linspace(aqua_color(2), blue_color(2), num_colors);
                                       linspace(aqua_color(3), blue_color(3), num_colors)].';
                                       sscanf('1F56A8','%2x%2x%2x',[1 3])/255];
% Set the custom colormap
colormap(Ycustom_color_map);

% Extract data from the table
x_values = y_df.Y_divided_by_max;
y_values = y_df.gamma;
color_values = y_df.t_ac;

% Define grid for contour plot
x_grid = linspace(min(x_values), max(x_values), 100);
y_grid = linspace(min(y_values), max(y_values), 100);
[X, Y] = meshgrid(x_grid, y_grid);

% Interpolate color values onto grid
Z_color = griddata(x_values, y_values, color_values, X, Y);
% Create a filled contour plot with the custom colormap
contourf(X, Y, Z_color, 20, 'LineStyle', ':', 'LineWidth', 1);
cb = colorbar;
ylabel(cb, 'Acetate Threshold (\lambda_{ac}) (h^{-1})', 'FontSize', 14);
set(gca, 'FontSize', 12)
% Set x-axis limits to show only the range from 0 to 1
xlim([0, 1]);

% Add labels and title
xlabel('\phi Y /\phi_{max}', 'FontSize', 14)
ylabel('\gamma ratio (\phi^{m}_{max}/\phi_{max})', 'FontSize', 14)

% Show the grid
% grid on;

% Save the figure as SVG
filename = sprintf('Contour_uProtYvsAceT_AlterV0.svg');
saveas(fig, filename, 'svg');


%%
 
%Plotting Z Contour Plot
fig = figure();

% Define the colors for the colormap
white_color = [1 1 1];   % White color
yellow_color = sscanf('FFD700','%2x%2x%2x',[1 3])/255;  % Yellow color
orange_color = sscanf('FF6F61','%2x%2x%2x',[1 3])/255;  % Yellow color
red_color = sscanf('DC143C','%2x%2x%2x',[1 3])/255;     % Red color

% Number of color points in the colormap
num_colors = 20;
% Create the custom colormap with a smooth gradient
Zcustom_color_map = [linspace(white_color(1), yellow_color(1), num_colors);
                    linspace(white_color(2), yellow_color(2), num_colors);
                    linspace(white_color(3), yellow_color(3), num_colors)].';
% Number of color points in the colormap
num_colors = 20;
% Append the red part of the colormap
Zcustom_color_map = [Zcustom_color_map; [linspace(yellow_color(1), orange_color(1), num_colors);
                                       linspace(yellow_color(2), orange_color(2), num_colors);
                                       linspace(yellow_color(3), orange_color(3), num_colors)].'];
% Number of color points in the colormap
num_colors = 20;
% Append the red part of the colormap
Zcustom_color_map = [Zcustom_color_map; [linspace(orange_color(1), red_color(1), num_colors);
                                       linspace(orange_color(2), red_color(2), num_colors);
                                       linspace(orange_color(3), red_color(3), num_colors)].';
                                       sscanf('C00217','%2x%2x%2x',[1 3])/255];
% Set the custom colormap
colormap(Zcustom_color_map);
% Extract data from the table
x_values = z_df.Z_divided_by_max;
y_values = z_df.gamma;
color_values = z_df.t_ac;

% Define grid for contour plot
x_grid = linspace(min(x_values), max(x_values), 100);
y_grid = linspace(min(y_values), max(y_values), 100);
[X, Y] = meshgrid(x_grid, y_grid);

% Interpolate color values onto grid
Z_color = griddata(x_values, y_values, color_values, X, Y);
% Create a filled contour plot with the custom colormap
contourf(X, Y, Z_color, 20, 'LineStyle', ':', 'LineWidth', 1);
colormap(Zcustom_color_map); % Apply the custom colormap
cb = colorbar;
ylabel(cb, 'Acetate Threshold (\lambda_{ac}) (h^{-1})', 'FontSize', 14);
% Set x-axis limits to show only the range from 0 to 1
xlim([0, 1]);
set(gca, 'FontSize', 12)
% Add labels and title
xlabel('\phi Z /\phi_{max}', 'FontSize', 14)
ylabel('\gamma ratio (\phi^{m}_{max}/\phi_{max})', 'FontSize', 14)

% Save the figure as SVG
filename = sprintf('Contour_uProtZvsAceT_AlterV0.svg');
saveas(fig, filename, 'svg');


%%
fig= figure();
% Create a surface plot
surf(X, Y, Z_color);
colormap(Zcustom_color_map); % Choose a colormap (e.g., jet)
colorbar; % Add a colorbar

% Add labels and title
% Add labels and title
xlabel('\phi Z /\phi_{max}', 'FontSize', 14)
ylabel('\gamma ratio (\phi^{m}_{max}/\phi_{max})', 'FontSize', 14)
zlabel('Acetate Threshold (\lambda_{ac})');

% Set axis limits if desired
xlim([min(x_values), max(x_values)]);
ylim([min(y_values), max(y_values)]);

% Save the figure as SVG
filename = sprintf('3D_uProtZvsAceT_AlterV0.svg');
saveas(fig, filename, 'svg');

%%

% PLOTING NORMALIZED DATA WITH THEIR OWN UPPER LIMITS



fig = figure('Units', 'centimeters', 'Position', [1, 1, 13, 12.5]);
hold on;
%     markers = {'o', '>', '^', 'd', 's', 'v', '<'}; %Z %%%% CHANGE THIS %%%%%
%     markers = {'-o', '->', '-^', '-d', '-s', '-v', '-<'}; %Z

markers = {'o', '^', 's', 'd', '>', 'v', '<'}; %Y %%%% CHANGE THIS %%%%%
colors = {
        sscanf('C00217','%2x%2x%2x',[1 3])/255,... % dark red
        sscanf('6B8E23','%2x%2x%2x',[1 3])/255,... % Olive Green
        sscanf('008080','%2x%2x%2x',[1 3])/255,... % Teal
        sscanf('40E0D0','%2x%2x%2x',[1 3])/255,... % Turquoise
        sscanf('87CEEB','%2x%2x%2x',[1 3])/255,... % Sky Blue
        sscanf('6495ED','%2x%2x%2x',[1 3])/255,... % Cornflower Blue
        sscanf('1F56A8','%2x%2x%2x',[1 3])/255,... % blue
        [0.07450980392156863, 0.06666666666666667, 0.0196078431372549],... %black
        sscanf('FFFFFF','%2x%2x%2x',[1 3])/255,... %white
             };        
         
% title(gamma_string)
% xlabel('\phi Z /\phi_{max} or \phi Y /\phi^{m}_{max}', 'FontSize', 14)
xlabel('\phi Z /\phi_{max} or \phi Y /\phi^{m}_{max}', 'FontSize', 14)

ylabel('Acetate Threshold (\lambda_{ac}) (h^{-1})', 'FontSize', 14)
set(gca, 'FontSize', 12)


%Plot Membrane proteins
ygammas = unique(y_df.gamma);
gamma_select =[1, 0.35,0.3,0.25,0.23,0.2] ;
gamma_select= gamma_select(end:-1:1);
for i=1:numel(gamma_select)
    gammaVal = gamma_select(i);
    sub_y_df = y_df(y_df.gamma== gammaVal,:); 
    sub_y_df = sortrows(sub_y_df, {'gamma','Y_divided_by_Mmax'});
    indexZero = find(sub_y_df.t_ac ==0,1);
    sub_y_df = sub_y_df(1:indexZero, :);
        
    plot(sub_y_df.Y_divided_by_Mmax, sub_y_df.t_ac, '-', 'Color', colors{1+i}, 'LineWidth', 4, 'HandleVisibility', 'off');
    % Downsample the data to only keep 20 points with markers
    downsampleFactor = 3;  % Keep every 5th point
    numPoints = numel(sub_y_df.t_ac);
    indicesToShow = 1 : downsampleFactor : numPoints;
       
    if gammaVal == 1
        indicesToShow = sort ([indicesToShow, 54,55]);
  
    end
    
    gammaLabel = sprintf('Membrane (\\phi Y) \\gamma %s', num2str(gammaVal)); 
    plot(sub_y_df.Y_divided_by_Mmax(indicesToShow), sub_y_df.t_ac(indicesToShow),...
    markers{1+i}, 'Color', colors{1+i}, 'LineWidth', 1.5,'MarkerEdgeColor', colors{1+i},...
    'MarkerFaceColor', colors{end}, 'MarkerSize', 8, 'DisplayName', gammaLabel);
end

% Plot Cytosolic Protein
sub_z_df = z_df(z_df.gamma== 1,:); 
plot(sub_z_df.Z_divided_by_max, sub_z_df.t_ac, '-', 'Color', colors{1}, 'LineWidth', 4, 'HandleVisibility', 'off');
% Downsample the data to only keep 20 points with markers
downsampleFactor = 4;  % Keep every 5th point
numPoints = numel(sub_z_df.t_ac);
indicesToShow = 1 : downsampleFactor : numPoints;
indicesToShow = sort ([indicesToShow, 52,51]);
h1=plot(sub_z_df.Z_divided_by_max(indicesToShow), sub_z_df.t_ac(indicesToShow),...
    markers{1}, 'Color', colors{1}, 'LineWidth', 1.5,'MarkerEdgeColor', colors{1},...
    'MarkerFaceColor', colors{end}, 'MarkerSize', 8, 'DisplayName', 'Cytosolic (\phi Z)');


xlim([0, 1]);
y = ylim;
% Set the legend
lgd = legend('Location', 'southoutside','NumColumns', 2);
set(lgd, 'Box', 'off', 'FontSize', 12);
% Show the grid
grid on;

% Save the figure as SVG
filename = sprintf('NormuProtvsAceT_AlterV0.svg');
saveas(fig, filename, 'svg');


%%

fig = figure('Units', 'centimeters', 'Position', [1, 1, 13, 12.5]);
hold on;

markers = {'o', '^', 's', 'd', '>', 'v', '<'}; %Y %%%% CHANGE THIS %%%%%
colors = {
        sscanf('6B8E23','%2x%2x%2x',[1 3])/255,... % Olive Green
        sscanf('008080','%2x%2x%2x',[1 3])/255,... % Teal
        sscanf('40E0D0','%2x%2x%2x',[1 3])/255,... % Turquoise
        sscanf('87CEEB','%2x%2x%2x',[1 3])/255,... % Sky Blue
        sscanf('6495ED','%2x%2x%2x',[1 3])/255,... % Cornflower Blue
        sscanf('1F56A8','%2x%2x%2x',[1 3])/255,... % blue
        [0.07450980392156863, 0.06666666666666667, 0.0196078431372549],... %black
        sscanf('FFFFFF','%2x%2x%2x',[1 3])/255,... %white
             };
         
xlabel('\phi Y /\phi^{m}_{max}', 'FontSize', 14)
ylabel('Acetate Threshold (\lambda_{ac}) (h^{-1})', 'FontSize', 14)
set(gca, 'FontSize', 12)

%Plot Membrane proteins
ygammas = unique(y_df.gamma);
gamma_select =[1,0.35,0.3,0.25,0.23,0.2] ;
gamma_select= gamma_select(end:-1:1);

for i=1:numel(gamma_select)
    gammaVal = gamma_select(i);
    sub_y_df = y_df(y_df.gamma== gammaVal,:); 
    sub_y_df = sortrows(sub_y_df, {'gamma','Y_divided_by_Mmax'});
    indexZero = find(sub_y_df.t_ac ==0,1);
    sub_y_df = sub_y_df(1:indexZero, :);
        
    plot(sub_y_df.Y_divided_by_Mmax, sub_y_df.t_ac, '-', 'Color', colors{i}, 'LineWidth', 4, 'HandleVisibility', 'off');
           
    if i == numel(gamma_select)
        downsampleFactor = 4;  % Keep every 5th point
        numPoints = numel(sub_y_df.t_ac);
        indicesToShow = 1 : downsampleFactor : numPoints;
        indicesToShow = sort ([indicesToShow, 54,55]);
    elseif i == numel(gamma_select)-1
        downsampleFactor = 4;  % Keep every 5th point
        numPoints = numel(sub_y_df.t_ac);
        indicesToShow = 1 : downsampleFactor : numPoints;
    else
        downsampleFactor = 3;  % Keep every 5th point
        numPoints = numel(sub_y_df.t_ac);
        indicesToShow = 1 : downsampleFactor : numPoints;
    end
    
    gammaLabel = sprintf('Membrane (\\phi Y) \\gamma %s', num2str(gammaVal)); 
    plot(sub_y_df.Y_divided_by_Mmax(indicesToShow), sub_y_df.t_ac(indicesToShow),...
    markers{i+1}, 'Color', colors{i}, 'LineWidth', 1.5,'MarkerEdgeColor', colors{i},...
    'MarkerFaceColor', colors{end}, 'MarkerSize', 8, 'DisplayName', gammaLabel);
end

xlim([0, 1]);
y = ylim;
% Set the legend
lgd = legend('Location', 'southoutside','NumColumns', 2);
set(lgd, 'Box', 'off', 'FontSize', 12);
% Show the grid
grid on;

% Save the figure as SVG
filename = sprintf('NormuProtYvsAceT_AlterV0.svg');
saveas(fig, filename, 'svg');


%%
% ygammas = unique(y_df.gamma);
% y_normMax = unique(y_df.Y_divided_by_max);
% sub_y_df = y_df(y_df.gamma== 0.25,:);
% disp(sub_y_df);

%%
%Plotting Y Contour Plot
fig = figure();

% Define the colors for the colormap
white_color = [1 1 1];   % White color
blue_color =  sscanf('6495ED','%2x%2x%2x',[1 3])/255;  % Blue color
green_color = sscanf('6B8E23','%2x%2x%2x',[1 3])/255;  % Green color
aqua_color = sscanf('40E0D0','%2x%2x%2x',[1 3])/255;     % Aqua color
% Number of color points in the colormap
num_colors = 30;
% Create the custom colormap with a smooth gradient
Ycustom_color_map = [linspace(white_color(1), green_color(1), num_colors);
                    linspace(white_color(2), green_color(2), num_colors);
                    linspace(white_color(3), green_color(3), num_colors)].';
% Number of color points in the colormap
num_colors = 30;
% Append the red part of the colormap
Ycustom_color_map = [Ycustom_color_map; [linspace(green_color(1), aqua_color(1), num_colors);
                                       linspace(green_color(2), aqua_color(2), num_colors);
                                       linspace(green_color(3), aqua_color(3), num_colors)].'];
% Number of color points in the colormap
num_colors = 30;
% Append the red part of the colormap
Ycustom_color_map = [Ycustom_color_map; [linspace(aqua_color(1), blue_color(1), num_colors);
                                       linspace(aqua_color(2), blue_color(2), num_colors);
                                       linspace(aqua_color(3), blue_color(3), num_colors)].';
                                       sscanf('1F56A8','%2x%2x%2x',[1 3])/255];
% Set the custom colormap
colormap(Ycustom_color_map);

% Extract data from the table
y_values = y_df.Y_divided_by_Mmax;
x_values = y_df.gamma;
color_values = y_df.t_ac;

% Define grid for contour plot
x_grid = linspace(min(x_values), max(x_values), 200);
y_grid = linspace(min(y_values), max(y_values), 200);
[X, Y] = meshgrid(x_grid, y_grid);

% Interpolate color values onto grid
Z_color = griddata(x_values, y_values, color_values, X, Y);
% Create a filled contour plot with the custom colormap
contourf(X, Y, Z_color, 100, 'LineStyle', ':', 'LineWidth', 1);
cb = colorbar;
ylabel(cb, 'Acetate Threshold (\lambda_{ac}) (h^{-1})', 'FontSize', 14);
set(gca, 'FontSize', 12)
% Set x-axis limits to show only the range from 0 to 1
ylim([0, 1]);

% Add labels and title
ylabel('\phi Y /\phi^{m}_{max}', 'FontSize', 14)
xlabel('\gamma ratio (\phi^{m}_{max}/\phi_{max})', 'FontSize', 14)

% Show the grid
% grid on;

% Save the figure as SVG
filename = sprintf('NormContour_uProtYvsAceT_AlterV0.svg');
saveas(fig, filename, 'svg');


%%

% Combine x_values and y_values into a single coordinate array
coordinates = [x_values, y_values];

% Find unique coordinates
unique_coords = unique(coordinates, 'rows');



%%
fig = figure('Units', 'centimeters', 'Position', [1, 1, 13, 12.5]);
hold on;

% Define the colors for the colormap
white_color = [1 1 1];   % White color
blue_color =  sscanf('6495ED','%2x%2x%2x',[1 3])/255;  % Blue color
green_color = sscanf('6B8E23','%2x%2x%2x',[1 3])/255;  % Green color
aqua_color = sscanf('40E0D0','%2x%2x%2x',[1 3])/255;     % Aqua color
% Number of color points in the colormap
num_colors = 10;
% Create the custom colormap with a smooth gradient
Ycustom_color_map = [linspace(white_color(1), green_color(1), num_colors);
                    linspace(white_color(2), green_color(2), num_colors);
                    linspace(white_color(3), green_color(3), num_colors)].';
% Number of color points in the colormap
num_colors = 25;
% Append the red part of the colormap
Ycustom_color_map = [Ycustom_color_map; [linspace(green_color(1), aqua_color(1), num_colors);
                                       linspace(green_color(2), aqua_color(2), num_colors);
                                       linspace(green_color(3), aqua_color(3), num_colors)].'];
% Number of color points in the colormap
num_colors = 42;
% Append the red part of the colormap
Ycustom_color_map = [Ycustom_color_map; [linspace(aqua_color(1), blue_color(1), num_colors);
                                       linspace(aqua_color(2), blue_color(2), num_colors);
                                       linspace(aqua_color(3), blue_color(3), num_colors)].';
                                       sscanf('1F56A8','%2x%2x%2x',[1 3])/255];
% Set the custom colormap
colormap(Ycustom_color_map);

xlabel('\phi Z /\phi_{max} or \phi Y /\phi^{m}_{max}', 'FontSize', 14)
ylabel('Acetate Threshold (\lambda_{ac}) (h^{-1})', 'FontSize', 14)

% Initialize an empty table to store concatenated sub_y_df tables
concatenatedTable = table();

% Depuring the data
ygammas = unique(y_df.gamma);
for i=1:numel(ygammas)
    gammaVal = ygammas(i);
    sub_y_df = y_df(y_df.gamma== gammaVal,:); 
    sub_y_df = sortrows(sub_y_df, {'gamma','Y_divided_by_Mmax'});
    indexZero = find(sub_y_df.t_ac ==0,1);
    sub_y_df = sub_y_df(1:indexZero, :);
    colorMapSize = size(Ycustom_color_map);
    indexColor = round(gammaVal*(colorMapSize(1)-1))+1;
  
    plot(sub_y_df.Y_divided_by_Mmax, sub_y_df.t_ac, '-', 'LineWidth', 4, 'HandleVisibility', 'off', 'Color', Ycustom_color_map(indexColor,:));
    
    % Vertically concatenate the current sub_y_df to the concatenatedTable
    concatenatedTable = vertcat(concatenatedTable, sub_y_df);

end

xlim([0, 1]);
y = ylim;
% Set the legend
lgd = legend('Location', 'southoutside','NumColumns', 2);
set(lgd, 'Box', 'off', 'FontSize', 12);
% Show the grid
grid on;
cb = colorbar;
ylabel(cb, '\gamma ratio (\phi^{m}_{max}/\phi_{max})', 'FontSize', 14);
set(gca, 'FontSize', 12)

% Save the figure as SVG
filename = sprintf('NormContuProtYvsAceT_AlterV0.svg');
saveas(fig, filename, 'svg');


%%
%Plotting Y Contour Plot
fig = figure();

% Define the colors for the colormap
white_color = [1 1 1];   % White color
blue_color =  sscanf('6495ED','%2x%2x%2x',[1 3])/255;  % Blue color
green_color = sscanf('6B8E23','%2x%2x%2x',[1 3])/255;  % Green color
aqua_color = sscanf('40E0D0','%2x%2x%2x',[1 3])/255;     % Aqua color
% Number of color points in the colormap
num_colors = 30;
% Create the custom colormap with a smooth gradient
Ycustom_color_map = [linspace(white_color(1), green_color(1), num_colors);
                    linspace(white_color(2), green_color(2), num_colors);
                    linspace(white_color(3), green_color(3), num_colors)].';
% Number of color points in the colormap
num_colors = 30;
% Append the red part of the colormap
Ycustom_color_map = [Ycustom_color_map; [linspace(green_color(1), aqua_color(1), num_colors);
                                       linspace(green_color(2), aqua_color(2), num_colors);
                                       linspace(green_color(3), aqua_color(3), num_colors)].'];
% Number of color points in the colormap
num_colors = 30;
% Append the red part of the colormap
Ycustom_color_map = [Ycustom_color_map; [linspace(aqua_color(1), blue_color(1), num_colors);
                                       linspace(aqua_color(2), blue_color(2), num_colors);
                                       linspace(aqua_color(3), blue_color(3), num_colors)].';
                                       sscanf('1F56A8','%2x%2x%2x',[1 3])/255];
% Set the custom colormap
colormap(Ycustom_color_map);

% Extract data from the table
x_values = concatenatedTable.Y_divided_by_Mmax;
y_values = concatenatedTable.gamma;
color_values = concatenatedTable.t_ac;

% Define grid for contour plot
x_grid = linspace(min(x_values), max(x_values), 100);
y_grid = linspace(min(y_values), max(y_values), 100);
[X, Y] = meshgrid(x_grid, y_grid);

% Interpolate color values onto grid
Z_color = griddata(x_values, y_values, color_values, X, Y);
% Create a filled contour plot with the custom colormap
contourf(X, Y, Z_color, 20, 'LineStyle', ':', 'LineWidth', 1);
cb = colorbar;
ylabel(cb, 'Acetate Threshold (\lambda_{ac}) (h^{-1})', 'FontSize', 14);
set(gca, 'FontSize', 12)
% Set x-axis limits to show only the range from 0 to 1
xlim([0, 1]);

% Add labels and title
xlabel('\phi Y /\phi^{m}_{max}', 'FontSize', 14)
ylabel('\gamma ratio (\phi^{m}_{max}/\phi_{max})', 'FontSize', 14)

% Show the grid
% grid on;

% Save the figure as SVG
filename = sprintf('NormContour_uProtYvsAceT_AlterV0.svg');
saveas(fig, filename, 'svg');