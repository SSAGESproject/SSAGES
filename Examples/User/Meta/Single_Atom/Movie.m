set(0,'defaulttextInterpreter','latex');
set(groot, 'defaultAxesTickLabelInterpreter','latex');
set(groot, 'defaultLegendInterpreter','latex');

clear vars;
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

% Create a Mesh to generate a plot on
% Bounds from lammps input is box from -3.14159 to 3.14159
[xg,yg] = meshgrid(-3.14159:0.1:3.14159, -3.14159:0.1:3.14159);
frames = struct('cdata',[],'colormap',[]);
ii=length(SSAGES_Data)-1;
while ii < length(SSAGES_Data)
    % Set local variables xcenter, ycenter, sigma, and height
    xc = SSAGES_Data(1:ii,1);
    yc = SSAGES_Data(1:ii,2);

    sigma = SSAGES_Data(1,3);
    height = SSAGES_Data(1,5);

    % Sum hills using gaussian function
    % SSAGES_Gaussian = @(x, y) sum(height*exp(-((x-xc).^2 + (y-yc).^2)./(2*sigma^2)));

    % Run the function on the mesh to get the free energy landscape 
    for i=1:size(xg,1)
        for j=1:size(yg,2)
            SSAGES_Z(i,j) = SG(xg(i,j),yg(i,j),xc,yc,height,sigma);
        end
    end

    % Generate the free energy landscape
    figure(2);
    surf(xg,yg,-SSAGES_Z);
    view(52.1287,21.6839);
    p = Plot(gcf, true);
    %p.ZLim = [-60, 0];
    p.BoxDim = [8,8];
    p.ZLabel='F ($K_bT$)';
    p.XLabel='X Position';
    p.YLabel='Y Position';
    
    set(gca,'FontSize',20);
    drawnow limitrate;
    frames(end+1) = getframe;
    ii=int16(ii*1.01+1);
    return;
end
frames(1) = [];