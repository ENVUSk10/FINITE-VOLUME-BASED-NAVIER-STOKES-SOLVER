% STREAM VORTICITY FORMULATION FOR LID DRIVEN CAVITY
clc
clear
close all
%% DEFINING THE MESH
domain_size = 1;
grid_points = 51;
h = domain_size/(grid_points-1);
U = 1;
nu = 1e-3;
Re = (U*domain_size)/nu;
dt = 0.001;
Beta = 1.95;
CFL = (U*dt)/h;
F = (nu*dt)/h^2;


if CFL^2 <=2*F && CFL <= 1 && F<=0.5
    disp("Stable")
else 
    disp("Unstable. Pls change the values")
end

x = linspace(0,domain_size,grid_points);
y = linspace(0,domain_size,grid_points);

error_req = 1e-7;
error_mag = 1;
iterations = 0;
%% INITIALIZING THE VARIABLES
psi = zeros(grid_points,grid_points);
w = zeros(grid_points,grid_points);
u = zeros(grid_points,grid_points);
u(1,:) = 1;
v = zeros(grid_points,grid_points);

oldpsi = psi;
oldw = w;
%% NUMERICAL SOLUTION
while error_mag > error_req
    for i = 2:grid_points-1
        for j = 2:grid_points-1
            psi(i,j) = 0.25*Beta*(psi(i+1,j) + psi(i-1,j) +psi(i,j+1) + ...
                       psi(i,j-1) + h^2*w(i,j)) + (1 - Beta)*(psi(i,j));
        end
    end
    
    for i = 2: grid_points-1
        for j = 2:grid_points-1
            w(1,j) = -2.0*psi(2,j)/h^2- ((2*U)/h);
            w(i,1) = -2.0*psi(i,2)/h^2;
            w(i,grid_points) = -2.0*psi(i, grid_points-1)/h^2 ;
            w(grid_points,j) = -2.0*psi(grid_points-1,j)/h^2;
        end
    end
    
    error_mag = 0;
    for i =2:grid_points-1
        for j= 2:grid_points-1
            w(i,j) = w(i,j) - (dt/(4*h^2))*((psi(i, j+1) - psi(i,j-1))*(w(i+1,j) - w(i-1,j))) + (dt/(4*h^2))...
                     *((psi(i+1,j) - psi(i-1,j))*(w(i,j+1) - w(i,j-1))) + (dt/(Re*h^2))*(w(i+1,j) + w(i-1,j)... 
                     + w(i,j+1) + w(i,j-1) -4*w(i,j));
            
             error_mag = error_mag + abs(w(i,j) - oldw(i,j)); % error calculation
        end
    end
    
    if(rem(iterations, 2000)) == 0
       figure(1);
       hold on
       plot(iterations, error_mag, '-ko')
       xlabel('Iterations')
       ylabel('Residual Error')
       set(gca,'yscale','log')
    end
    
    oldpsi = psi;
    oldw = w;
    iterations = iterations + 1;
end
 %% POST-PROCESSING
for i =2:grid_points-1
    for j= 2:grid_points-1
        u(i,j) = -(psi(i+1,j) - psi(i-1,j))/(2*h);
        v(i,j) = (psi(i,j-1) - psi(i,j-1))/(2*h);
     end
end
plot(u(:,((grid_points+1)/2)),1-y, 'LineWidth', 1)
contour(x,1-y,psi)   
        
        
        
        
        
        
        
        
        
        
        
        
        
        










