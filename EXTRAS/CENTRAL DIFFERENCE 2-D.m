clc
close all
clear all
%% INITIALIZING 

domain_size=1;
grid_points=51;
h=domain_size/(grid_points-1);

x=0:h:domain_size;

T(grid_points,grid_points)=0;
T_new(grid_points,grid_points)=0;

T(1,1:grid_points)=1;
T(1:grid_points,1)=1;


T_new(1,1:grid_points)=1;
T_new(1:grid_points,1)=1;

rho=1
u=1
v=1
gamma=0.002
p=rho*u*h/gamma;
dx=h;
dy=h;

error_mag=1;
error_req=1e-7;
iteration=0;

%% CENTRAL DIFFERENCING

a_e=gamma-rho*u*h/2;
a_n=gamma-rho*v*h/2;
a_w=gamma+rho*u*h/2;
a_s=gamma+rho*v*h/2;
a_p=a_s+a_n+a_w+a_e;

while error_mag>error_req
    for i=2:grid_points-1
        for j=2:grid_points-1
            T_new(i,j)=(a_n*T(i,j-1)+a_s*T(i,j+1)+a_w*T(i-1,j)+a_e*T(i+1,j))/a_p;
        end
    end
    iteration=iteration+1;
    error_mag=0;
    for i=2:grid_points-1
        for j=2:grid_points-1
            error_mag=error_mag+abs(T(i,j)-T_new(i,j));
        end
        
    end

    T=T_new;
end

%%
x=((1:grid_points)-1).*h
y=1-((1:grid_points)-1).*h
[X,Y]=meshgrid(x,y);
contourf(X,Y,T,50);
colorbar

plot(y,T(:,(grid_points+1)/2),'--o')


% A 1 2 3 4 5 B








