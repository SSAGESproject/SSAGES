% Import original data
if exist('fes.dat','file')
    Original_Data = importdata('fes.dat');
else
    disp('File fes.dat does not exist, ignoring original results');
    Original_Data = [];
end

% Import SSAGES metadynamics data in the form:
% xcenter, ycenter, sigma, sigma, height
if exist('hills.out','file')
    SSAGES_Data = importdata('hills.out');
else
    clear Original_Data;
    error('File hills.out does not exist, exiting');
end

if isempty(SSAGES_Data)
    clear Original_Data;
    error('No data in hills.out! Exiting'); 
end

% Set local variables xcenter, ycenter, sigma, and height
xc = SSAGES_Data(:,1);
yc = SSAGES_Data(:,2);

sigma = SSAGES_Data(1,3);
height = SSAGES_Data(1,5);

% Sum hills using gaussian function
SSAGES_Gaussian = @(x, y) sum(height*exp(-((x-xc).^2 + (y-yc).^2)./(2*sigma^2)));

% Create a Mesh to generate a plot on
% Bounds from lammps input is box from -1.5 to 2.0
[xg,yg] = meshgrid(-1.5:.05:2, -1.5:.05:2);

% Run the function on the mesh to get the free energy landscape 
for i=1:size(xg,1)
    for j=1:size(xg,2)
        SSAGES_Z(i,j) = SSAGES_Gaussian(xg(i,j),yg(i,j));
    end
end

% Generate the inverse free energy landscape
figure(1);
surf(xg,yg,SSAGES_Z);
title('Inverse free energy SSAGES');
xlabel('x position')
ylabel('y position') 
zlabel('K_bT') 

% Generate the free energy landscape
figure(2);
surf(xg,yg,-SSAGES_Z);
title('Free energy SSAGES');
xlabel('x position')
ylabel('y position') 
zlabel('K_bT') 

% Generate the free energy landscape from original data
if ~isempty(Original_Data)
    for i = 1:100
        for j = 1:100
            Original_Z(i,j) = Original_Data(100*(i-1)+j,3);
        end
    end
    Original_X = unique(Original_Data(:,1));
    Original_Y = unique(Original_Data(:,2));
    figure(3);
    surf(Original_X,Original_Y,-Original_Z);
    title('Free energy Reference');
    xlabel('x position')
    ylabel('y position') 
    zlabel('K_bT') 
end