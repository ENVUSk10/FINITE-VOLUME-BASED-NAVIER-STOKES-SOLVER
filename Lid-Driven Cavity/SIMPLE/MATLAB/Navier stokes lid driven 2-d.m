clear all
close all
clc
made changes in Master from newname
%% DEFINING THE MESH
grid_points =51;
dom_length = 1;
h = dom_length/(grid_points-1);
Re = 200; 
nu = 1/Re;
alpha = 0.8;
%% Initializing the variables

u_final(grid_points,grid_points)=0;
v_final(grid_points,grid_points)=0;
p_final(grid_points,grid_points)=1;
u_final(1,:)=2;


u(grid_points+1,grid_points)=0;
u_star(grid_points+1,grid_points)=0;
de(grid_points+1,grid_points)=0;

v(grid_points,grid_points+1)=0;
v_star(grid_points,grid_points+1)=0;
dn(grid_points,grid_points+1)=0;

p(grid_points+1,grid_points+1)=1;
pressure_ccorrection(grid_points+1,grid_points+1)=0;
b(grid_points+1,grid_points+1)=0;

u(1,:)=2;

u_new(grid_points+1,grid_points)=0;
v_new(grid_points,grid_points+1)=0;
p_new(grid_points+1,grid_points+1)=1;
u_new(1,:)=2;
%% Solving the governing equations
error_mag=1;
iterations=0;
error_req=1e-7;

while error_mag > error_req
   
    for i = 2:grid_points
        for j = 2:grid_points - 1
            uE = 0.5*(u(i,j) + u(i,j+1));
            uW = 0.5*(u(i,j) + u(i,j-1));
            vN = 0.5*(v(i-1,j) + v(i-1,j+1));
            vS = 0.5*(v(i,j) + v(i,j+1));
            
            aE = -0.5*uE*h + nu;
            aW = 0.5*uW*h + nu;
            aN = -0.5*vN*h + nu;
            aS = 0.5*vS*h + nu;
            ae = 0.5*uE*h - 0.5*uW*h + 0.5*vN*h - 0.5*vS*h + 4*nu;
            de(i,j) = -h/ae;
            
            u_star(i,j) = (aE*u(i,j+1) + aW*u(i,j-1) + aN*u(i-1,j) + aS*u(i+1,j))/ae + de(i,j)*(p(i,j+1) - p(i,j));
        end
    end
    
    
    u_star(1,:)=2-u_star(2,:);
    u_star(grid_points + 1,:)=-u_star(grid_points,:);
    u_star(2:grid_points,1)=0;
    u_star(2:grid_points,grid_points)=0;
    
    
    for i=2:grid_points - 1
        for j=2:grid_points
            uE=0.5*(u(i,j)+u(i+1,j));
            uW=0.5*(u(i,j-1)+u(i+1,j-1));
            vN=0.5*(v(i-1,j)+v(i,j));
            vS=0.5*(v(i,j)+v(i+1,j));
            
            aE=-0.5*uE*h + nu;
            aW=0.5*uW*h + nu;
            aN=-0.5*vN*h + nu;
            aS=0.5*vS*h + nu;
            an=0.5*uE*h-0.5*uW*h+0.5*vN*h-0.5*vS*h+4*nu;
            dn(i,j) = -h/an;
            
            v_star(i,j) = (aE*v(i,j+1)+aW*v(i,j-1)+aN*v(i-1,j)+aS*v(i+1,j))/an +dn(i,j)*(p(i,j)-p(i+1,j));
        end
    end
    

    v_star(:,1)= -v_star(:,2);
    v_star(:,grid_points + 1)=-v_star(:,grid_points);
    v_star(1,2:grid_points)=0;
    v_star(grid_points,2:grid_points)=0;
    
  
    pressure_ccorrection(1:grid_points+1,1:grid_points+1)=0;
    
 
    for i = 2:grid_points
        for j = 2:grid_points
            aE=-de(i,j)*h;
            aW=-de(i,j-1)*h;
            aN=-dn(i-1,j)*h;
            aS=-dn(i,j)*h;
            a_P=aE+aW+aN+aS;
            b(i,j) = -(u_star(i,j) - u_star(i,j-1))*h + (v_star(i,j) - v_star(i-1,j))*h;
            
            pressure_ccorrection(i,j)=(aE*pressure_ccorrection(i,j+1)+aW*pressure_ccorrection(i,j-1)+aN*pressure_ccorrection(i-1,j)+aS*pressure_ccorrection(i+1,j)+b(i,j))/a_P;
        end
    end
    
    
    for i=2:grid_points
        for j=2:grid_points
            p_new(i,j)=p(i,j)+alpha*pressure_ccorrection(i,j);
        end
    end
    
   
    p_new(1,:)=p_new(2,:);
    p_new(grid_points + 1,:)=p_new(grid_points,:);
    p_new(:,1)=p_new(:,2);
    p_new(:,grid_points + 1)=p_new(:,grid_points);
    
    
    for i=2:grid_points
        for j=2:grid_points-1
            u_new(i,j)=u_star(i,j)+alpha*de(i,j)*(pressure_ccorrection(i,j+1)-pressure_ccorrection(i,j));
        end
    end

    
  
    u_new(1,:)=2-u_new(2,:);
    u_new(grid_points + 1,:)=-u_new(grid_points,:);
    u_new(2:grid_points,1)=0;
    u_new(2:grid_points,grid_points)=0;
    
    for i=2:grid_points - 1
        for j=2:grid_points
            v_new(i,j)=v_star(i,j)+alpha*dn(i,j)*(pressure_ccorrection(i,j) - pressure_ccorrection(i+1,j));
        end
    end
    
    
    v_new(:,1)=-v_new(:,2);
    v_new(:,grid_points + 1)=-v_new(:,grid_points);
    v_new(1,2:grid_points)=0;
    v_new(grid_points,2:grid_points)=0;
            
    

    error_mag = 0;
    for i = 2:grid_points
        for j = 2:grid_points
            error_mag = error_mag + abs(b(i,j));
        end
    end
    u = u_new;
    v = v_new;
    p = p_new;
    iterations = iterations + 1;
    
end

for i = 1:grid_points
    for j = 1:grid_points
        u_final(i,j) = 0.5*(u(i,j) + u(i+1,j));
        v_final(i,j) = 0.5*(v(i,j) + v(i,j+1));
        p_final(i,j) = 0.25*(p(i,j) + p(i,j+1) + p(i+1,j) + p(i+1,j+1));
    end
end
%% VISUALIZATION
figure(3);
hold on;
y_dom = 1-((1:grid_points)-1).*h;
plot(u_final(:,(grid_points+1)/2),1-y_dom, 'LineWidth', 1)
legend('Re=200,102x102','Re=100,52x52')
print(gcf,'centrleline.png','-dpng','-r300dpi')

x_dom = ((1:grid_points)-1).*h;
y_dom = 1-((1:grid_points)-1).*h;
[X,Y] = meshgrid(x_dom,y_dom);
contourf(X,Y,v_final, 21, 'LineStyle', 'none')
colorbar
colormap('jet')
xlabel('x')
ylabel('y')

figure(3);
hold on;
quiver(X, Y, u_final, v_final, 5, 'k')
print(gcf,'re100','-dpng','-r300')