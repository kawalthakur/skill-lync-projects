% Usteady State Heat conduction by Explicit method 
clear 
clc
close all

% Defining Number of Intervals  
Nx = 20;
Ny = 20;

%Defining the number of time steps 
nt = 2000;

%Defining Diffusivity Constant 
alpha = 1.4;


% Number of grid points 
nx = Nx+1;
ny = Ny+1;

% Rectangular Dimensions 
Lx = 1;
Ly = 1;

% Grid spacing 
dx = Lx/Nx;
dy = Ly/Ny;

% Defining dt by the CFL criteria for explicit method 
% CFL = (2*alpha*dt)/(dx^2)
% CFL < 0.5 for stability 
CFL = 0.5;
dt = (CFL*dx^2)/(2*alpha);

%dt = 0.0002;
D = alpha*dt/(dx^2);

x = (0:Nx)*dx;
y = (0:Ny)*dy;
[X,Y] = meshgrid(x,y);

% Error criteria
eps = 1.e-4;

% Initial value for the temperature
temp =300*zeros(nx,ny,nt);

% Defining the Boundary Conditions Temperature in (K)
temp(:,ny,:) = 600;   % Top
temp(:,1,:)  = 900;   % Bottom
temp(1,:,:)  = 400;   % Left
temp(nx,:,:) = 800;   % Right
temp(1,1,:)   = (900+400)/2;
temp(1,ny,:)  = (600+400)/2;
temp(nx,1,:)  = (900+800)/2;
temp(nx,ny,:) = (600+800)/2;

temp_old = temp(:,:,1);

% Solving by Explicit method 
for t = 1:nt
    for i = 2:nx-1
        for j = 2:ny-1
            % Solving by Explicit Method 
            temp(i,j,t+1)=(1-4*D)*temp_old(i,j)+...
                D*(temp_old(i-1,j)+temp_old(i+1,j)+temp_old(i,j-1)+temp_old(i,j+1)); 
        end
    end
    temp_old = temp(:,:,t);
end

% Plotting the results by a countour plot 
i = 0;
time = 0;
for t = 1:10:nt-1
    i = i + 1;
    % Plotting the countour plot 
    [c,h]=contourf(Y,X,temp(:,:,t),'ShowText','on');
    c = colorbar;
    colormap(jet);
    c.FontSize = 11;
    set(get(c,'label'),'string','Contour:Temperature(K)'...
        ,'FontSize',11);
    xlabel('X(m)','Fontsize',14);
    ylabel('Y(m)','Fontsize',14);
    title(sprintf("Unsteady State Heat Conduction(Explicit) \n Time - %f sec ",time));
    %h.LineColor = 'none';
    %clabel(c,h);
    pause(0.003);
    M(i) = getframe(gcf);
    time = dt+dt*(t-1);
    if(max(max(abs(temp(:,:,t)-temp(:,:,t+1))))<eps)
        % This means that system has reached a Steady State 
        break;
    end
    
end  

%Creating a Video File with decreased Frame Rate for the animation 
videofile = VideoWriter('Unsteady_state_heat_conduction_explicit.avi','Uncompressed AVI');
videofile.FrameRate = 9;
open(videofile);
writeVideo(videofile,M);
close(videofile);