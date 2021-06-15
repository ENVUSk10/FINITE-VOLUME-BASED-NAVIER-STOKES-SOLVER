clc
close all
clear all
%% INITIALISING 
domain_size=1;
node_points=42;
h=domain_size/(node_points-1);

x=0:h:domain_size;

phi(1)=0;
phi(node_points)=1;

phi_new(1)=0;
phi_new(node_points)=1;

rho=1;
u=1;
gamma=0.02;
p=rho*u*domain_size/gamma;

error_mag=1;
error_req=1e-6;
iteration=0;
%% EXACT SOLUTION
phi_exact=phi(1)+(phi(node_points)-phi(1))/(exp(p)-1)*(exp(p*x/domain_size)-1)
figure(1);
plot(x,phi_exact,'ro');
a_p=rho*u/2-rho*u/2+gamma/h+gamma/h;
a_e=gamma/h-rho*u/2;
a_w=gamma/h+rho*u/2;
%% CENTRAL DIFFERENCE
while error_mag>error_req;
    
    for i=2:node_points-1
        phi_new(i)=(a_e*phi(i+1)+a_w*phi(i-1))/a_p;
    end
    iteration=iteration+1;
    error_mag=0;
    for i=2:node_points-1;
        error_mag=error_mag+abs(phi(i)-phi_new(i));
    end   
 phi=phi_new;
end
figure(1);
hold on
plot()
%% UPWIND
error=1;
req=1e-6
a_p=rho*u+gamma/h+gamma/h;
a_e=gamma/h;
a_w=gamma/h+rho*u;

while error>req;
    
    for i=2:node_points-1
        phi_new(i)=(a_e*phi(i+1)+a_w*phi(i-1))/a_p;
    end
    
    error=0;
    for i=2:node_points-1;
        error_mag=error_mag+abs(phi(i)-phi_new(i));
    end   
 phi=phi_new;
end
figure(1);
hold on
plot(x,phi,'b','linewidth',2)






