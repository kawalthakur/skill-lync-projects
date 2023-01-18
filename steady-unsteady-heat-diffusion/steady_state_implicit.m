% Solving Steady State Heat Conduction Equation
close all
clear 
clc

% Defining Number of Intervals  
Nx = 100;
Ny = 100;

% Number of grid points 
nx = Nx+1;
ny = Ny+1;

% Rectangular Dimensions 
Lx = 1;
Ly = 1;

% Grid spacing 
dx = Lx/Nx;
dy = Ly/Ny;

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


% Temperature value for the grid points 
temp = zeros(nx,ny);


% Defining the Boundary Conditions Temperature in (K)
temp(:,ny) = 600;   % Top
temp(:,1)  = 900;   % Bottom
temp(1,:)  = 400;   % Left
temp(nx,:) = 800;   % Right
temp(1,1)   = (900+400)/2;
temp(1,ny)  = (600+400)/2;
temp(nx,1)  = (900+800)/2;
temp(nx,ny) = (600+800)/2;

% Solving with Iterative Method
temp_old = temp;    % Initial guess for the values of temperature 
error = 2*eps;      % Error value at the beginning 
ncount = 0;         % Iteration counter 

while(error>eps)
    ncount = ncount + 1;
    for i = 2:nx-1
        for j = 2:ny-1
            if(method==1)
                % Solving by Gauss Jacobi Method
                temp(i,j) = (1/4)*(temp_old(i-1,j)+temp_old(i+1,j)...
                    +temp_old(i,j-1)+temp_old(i,j+1));
            end
          
            if(method==2)
                % Solving By Gauss Siedel Method
                temp(i,j) = (1/4)*(temp(i-1,j)+temp_old(i+1,j)...
                    +temp(i,j-1)+temp_old(i,j+1));
            end
            
            if(method==3)
                % SOR Method with Gauss Siedel method
                omega = 1.95;
                temp(i,j) = (1-omega)*temp_old(i,j)+ (omega/4)*(...
                    temp_old(i+1,j)+temp(i-1,j)+temp_old(i,j+1)+...
                    temp(i,j-1));
            end
            
        end
    end
    
    error = max(max(abs(temp(:,:)-temp_old(:,:))));
    
    if any(any(isinf(temp(:,:))) | any(isnan(temp(:,:))))
        % To check whether the Iterations are Diverging 
        % If Diverging then stop the program
        fprintf("Iteration Diverged.!");
        return;
    end
    
    % Updating the pervious value of temperature
    temp_old = temp;
    
    % Displaying the Iteration count with the Error 
    fprintf("\n%g - %e",ncount,error);
    
end
fprintf("\n Number of Iterations - %g",ncount);
% Plotting the countour plot 
[c,h]=contourf(Y,X,temp(:,:),'ShowText','on');
c = colorbar;
colormap(jet);
c.FontSize = 11;
set(get(c,'label'),'string','Contour:Temperature(K)'...
    ,'FontSize',11);
xlabel('X(m)','Fontsize',14);
ylabel('Y(m)','Fontsize',14);
%h.LineColor = 'none';
%clabel(c,h);
title(sprintf('Steady State Heat Conduction over a Rectangular Domain \n No. of Iterations = %d, %s'...
,ncount,type));