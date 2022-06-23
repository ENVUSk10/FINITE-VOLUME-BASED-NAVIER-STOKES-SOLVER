clc
close all
clear all
%% CREATING THE MESH
domain_size=1;
grid_points=31;
h=domain_size/(grid_points-1);
dt=0.0001;

T(grid_points,grid_points)=0;
T_new(grid_points,grid_points)=0;

T_transient=0;

T(:,1)=100;
T(grid_points,:)=200;
T(:,grid_points)=300;
T(1,:)=400;

T_new(:,1)=100;
T_new(grid_points,:)=200;
T_new(:,grid_points)=300;
T_new(1,:)=400;

error_mag=1;
error_req=1e-6;
iteration=0;

error_record=0;
%% CALCULATIONS
while error_mag>error_req
    for i=2:grid_points-1
        for j=2:grid_points-1
            T_new(i,j)=T_new(i,j)+(T_new(i+1,j)+T_new(i-1,j)+T_new(i,j+1)+T_new(i,j-1)-4*T_new(i,j))*(dt/(h*h));
            T_transient(iteration+1,1:grid_points,1:grid_points)=T_new;
        end
    end
    iteration=iteration+1;
    error_mag=0; 
    for i=2:grid_points-1
        for j=2:grid_points-1
            error_mag=error_mag+abs(T(i,j)-T_new(i,j));
            error_record(iteration)=error_mag;
        end
    end
    
    if rem(iteration,1000)==0
        iteration
        error_record(iteration)
    end
T=T_new;
end
%% VISUALIZATION
x=((1:grid_points)-1).*h;
y=1-((1:grid_points)-1).*h;
[X,Y]=meshgrid(x,y);
contourf(X,Y,T,12)
colorbar
%% TIMEBASED VISUALIZATION
time=(1:iteration).*dt
plot(time,error_record)

% FOR FINDING A PLOT AT PATICULAR TIME-INTERVAL
time_iteration=8185;
x=((1:grid_points)-1).*h;
y=1-((1:grid_points)-1).*h;
T_time_iteration=T_transient(time_iteration,:,:);
T_time_iteration=reshape(T_time_iteration,[grid_points,grid_points]);
contourf(X,Y,T_time_iteration,12)
colorbar;

% ANIMATION OF EVENTUAL STEADY STATE FORMATION
array=1:100:iteration
for i=1:length(array)
    time_iteration=array(i);
    T_time_iteration=T_transient(time_iteration,:,:);
    T_time_iteration=reshape(T_time_iteration,[grid_points,grid_points]);
    contourf(X,Y,T_time_iteration,12)
    colorbar;
    pause(1);
end



