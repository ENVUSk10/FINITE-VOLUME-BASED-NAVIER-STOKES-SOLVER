%FINITE VOLUME BASED NAVIER STOKES SOLVER FOR SIMULATING A FLOW THROUGH LID
%DRIVEN CAVITY USING THE ARTIFICIAL COMPRESSIBILITY METHOD
clc
close all
clear all
%% CREATING THE MESH
domain_size=1;
grid_points=41;
h=domain_size/(grid_points-1);

dt=0.001;
Re=100;
delta=4.5;
%% INITIALIZING THE PROBLEMS
%collocated variables
u_final(grid_points,grid_points)=0;
v_final(grid_points,grid_points)=0;
p_final(grid_points,grid_points)=1;
u_final(1,:)=1;
%staggered variables
u(grid_points+1,grid_points)=0;
v(grid_points,grid_points+1)=0;
p(grid_points+1,grid_points+1)=1;
u(1,:)=2;
% new staggered
u_new(grid_points+1,grid_points)=0;
v_new(grid_points,grid_points+1)=0;
p_new(grid_points+1,grid_points+1)=1;
u_new(1,:)=2;
%error analysis
error=1;
error_req=1e-6;
iteration=0;
error_record=0;
%% ARTIFICIAL COMPRESSSION METHOD
figure(1);
while error > error_req
    for i=2:grid_points
        for j=2:grid_points-1
            pressure = -(p(i,j+1) - p(i,j))/h;
            diffusion = (1/Re)*((u(i+1,j) - 2*u(i,j) + u(i-1,j))/(h*h) + (u(i,j+1) - 2*u(i,j) + u(i,j-1))/(h*h));
            advection_x = ((0.5*(u(i,j)+u(i,j+1)))^2 - (0.5*(u(i,j)+u(i,j-1)))^2)/h;
            advection_y = ((0.25*(u(i,j)+u(i-1,j))*(v(i-1,j)+v(i-1,j+1))) - (0.25*(u(i,j)+u(i+1,j))*(v(i,j)+v(i,j+1))))/h;
            u_new(i,j) = u(i,j) + dt*(diffusion - advection_x - advection_y + pressure);
        end
    end
    % U--> boundary conditions
     u_new(1,:) = 2 - u_new(2,:);
     u_new(grid_points + 1,:) = -u_new(grid_points,:);
     u_new(2:grid_points,1) = 0;
     u_new(2:grid_points,grid_points) = 0;
    
     for i=2:grid_points-1
        for j=2:grid_points
            pressure = -(p(i,j) - p(i+1,j))/h;
            diffusion = (1/Re)*((v(i+1,j) - 2*v(i,j) + v(i-1,j))/(h*h) + (v(i,j+1) - 2*v(i,j) + v(i,j-1))/(h*h));
            advection_y = ((0.5*(v(i,j)+v(i-1,j)))^2 - (0.5*(v(i,j)+v(i+1,j)))^2)/h;
            advection_x = ((0.25*(u(i,j)+u(i+1,j))*(v(i,j)+v(i,j+1))) - (0.25*(u(i,j-1)+u(i+1,j-1))*(v(i,j)+v(i,j-1))))/h;
            v_new(i,j) = v(i,j) + dt*(diffusion - advection_x - advection_y + pressure);
        end
     end
     % V--> boundary conditions
     v_new(:,1) = -v_new(:,2);
     v_new(:,grid_points + 1) = -v_new(:,grid_points);
     v_new(1,2:grid_points) = 0;
     v_new(grid_points,2:grid_points) = 0;
    for i=2:grid_points
        for j=2:grid_points
          p_new(i,j) = p(i,j) - delta*dt*(u(i,j) - u(i,j-1) + v(i-1,j) - v(i,j))/h;
        end
    end
    % P--> boundary conditions
    p_new(1,:) = p_new(2,:);
    p_new(grid_points + 1,:) = p_new(grid_points,:);
    p_new(:,1) = p_new(:,2);
    p_new(:,grid_points + 1) = p_new(:,grid_points);
    
    
    error=0;
    iteration=iteration+1;
    
    for i=2:grid_points-1
        for j=2:grid_points-1
            error = error + abs((u_new(i,j) - u_new(i,j-1) + v_new(i-1,j) - v_new(i,j))/h);
        end
    end
    % Convergence visualization
    if rem(iteration,5000)==0
        figure(1);
        semilogy(iteration,error,'-ko');
        hold on;
    end
    
    
    u=u_new;
    v=v_new;
    p=p_new;
end    
% substituting back to the collocated grid  
for i=1:grid_points
    for j=2:grid_points
        u_final(i,j)=(u(i,j)+u(i+1,j))/2;
        v_final(i,j)=(v(i,j)+v(i,j+1))/2;
        p_final(i,j)=(p(i,j)+p(i+1,j+1)+p(i+1,j)+p(i,j+1))/4;

    end
end
%% VISUALIZATION
x=((1:grid_points)-1).*h;
y=1-((1:grid_points)-1).*h;

[X,Y]=meshgrid(x,y)
contourf(X,Y,u_final,21,'linestyle','none')
colormap(jet)
colorbar

quiver(X,Y,u_final,v_final,4,'k')





