% Unsteady State heat conduction Equation 
clear 
clc
close all

% Defining Number of Intervals  
Nx = 100;
Ny = 100;

%Defining the number of time steps 
nt = 500;

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

% No Stability criteria for an Implicit method 
% Implicit Methods are inherentily Stable
dt = 1.e-3;
D = alpha*dt/(dx^2);

x = (0:Nx)*dx;
y = (0:Ny)*dy;
[X,Y] = meshgrid(x,y);

% Convergence Criteria / Tolerance 
eps = 1.e-4;

% Method Number 
% 1 - Gauss Jacobi Method 
% 2 - Gauss Siedel Method
% 3 - SOR Method using Gauss Siedel
method = 3;
if(method==1)
    type = 'Gauss Jacobi';
end
if(method==2)
    type = 'Gauss Siedel';
end
if(method==3)
    type = 'SOR method';
end

%Iterative Method 
temp = zeros(nx,ny,nt);

% Defining the Boundary Conditions Temperature in (K)
temp(:,ny,:) = 600;   % Top
temp(:,1,:)  = 900;   % Bottom
temp(1,:,:)  = 400;   % Left
temp(nx,:,:) = 800;   % Right

% Solving by Implicit Scheme
for t=2:nt
    error = 2;
    temp(:,:,t)=temp(:,:,t-1);
    temp_old = temp(:,:,t);
    ncount = 0;
    while(error > eps)
        ncount = ncount + 1;
        for i = 2:nx-1
            for j = 2:ny-1
                if(method==1)
                    %Gauss Jacobi Method
                    temp(i,j,t)= (1/(1+4*D))*(temp(i,j,t-1)+(D)*(temp_old(i-1,j)...
                       +temp_old(i+1,j)+temp_old(i,j-1)+temp_old(i,j+1))); 
                end
                if(method==2)
                    %Gauss Siedel Method
                    temp(i,j,t)= (1/(1+4*D))*(temp(i,j,t-1)+(D)*(temp(i-1,j,t)...
                       +temp_old(i+1,j)+temp(i,j-1,t)+temp_old(i,j+1)));               
                end
                if(method==3)
                    % SOR Method
                    omega  = 1.85;
                    DT = temp(i-1,j,t)+temp_old(i+1,j)+temp(i,j-1,t)+temp_old(i,j+1);
                    temp(i,j,t)=(1-omega)*temp_old(i,j)+(omega/(1+4*D))*(temp(i,j,t-1)+D*(DT));
                    %temp(i,j,t) = (1/(1-omega*(1+4*D)))*(temp_old(i,j)+omega*(temp(i,j,t-1)+(D)*(temp(i-1,j,t)...
                       %+temp_old(i+1,j)+temp(i,j-1,t)+temp_old(i,j+1))));
                end
                
            end
        end
        error = max(max(abs(temp(:,:,t)-temp_old)));
        if any(any(isinf(temp(:,:,t))) | any(isnan(temp(:,:,t))))
            fprintf("Iterations Diverged.!");
            return;
        end
        temp_old = temp(:,:,t);
        fprintf('%g - %e\n',ncount,error);
    end
    fprintf('\n%g\n',ncount);
    Iterations(t-1) = ncount;
end

% Plotting the results by a countour plot 
time = 0;
i = 0;
for t = 1:4:nt-1
     i = i + 1;
    % Plotting the countour plot 
    [C,h]=contourf(Y,X,temp(:,:,t),'ShowText','on');
    c = colorbar;
    colormap(jet);
    c.FontSize = 11;
    clabel(C,h,'LabelSpacing',150);
    set(get(c,'label'),'string','Contour:Temperature(K)'...
        ,'FontSize',11);
    xlabel('X(m)','Fontsize',14);
    ylabel('Y(m)','Fontsize',14);
    title(sprintf("Unsteady State Heat Conduction- %s \n Time = %f sec ",type,time));
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
fprintf("Average Number of Iterations per time steps = %g\n",round(mean(Iterations)));

%Creating a Video File with decreased Frame Rate for the animation 
videofile = VideoWriter('Unsteady_state_heat_conduction_Implicit.avi','Uncompressed AVI');
videofile.FrameRate = 9;
open(videofile);
writeVideo(videofile,M);
close(videofile);   