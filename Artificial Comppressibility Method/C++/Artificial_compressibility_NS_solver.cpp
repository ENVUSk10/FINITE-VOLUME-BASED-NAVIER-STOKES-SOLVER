/*
2-D incompressible flow Navier-Stokes solver using the Artificial Compressibility method.
Note : This code won't run for higher Reynolds numbers, for that you must to change dt and delta
for Re = 400 try to run with use dt = 0.0001 delta = 100.0
*/
#include<bits/stdc++.h>
using namespace std;
int main()
{
// Defining the Mesh
double domain_size = 1, grids = 51;    
double dx = domain_size/(grids-1);
double dy = dx;
double Beta = dx/dy;
double Re = 100.0;
double dt = 0.001;
double delta = 4.5;
double pressure, advection_x, advection_y, diffusion;
int grid_points = 51;

// Error definition
double error_mag = 1,error_req = 1e-8;
int iterations = 1;

// Declaration of the variables (Collocated and Staggered)
double u_final[grid_points][grid_points], v_final[grid_points][grid_points], p_final[grid_points][grid_points];
double u[grid_points+1][grid_points], v[grid_points][grid_points+1], p[grid_points+1][grid_points+1];
double u_new[grid_points+1][grid_points], v_new[grid_points][grid_points+1], p_new[grid_points+1][grid_points+1];
double x[grid_points], y[grid_points];

// Initializing Collocated variables
for (int i = 0; i<grid_points; i++)
{
	for(int j = 0; j<grid_points; j++)
	{
		
		if (i == 0)					  // Top wall
		{
			u_final[i][j] = 1.0;
			v_final[i][j] = 0.0;
			p_final[i][j] = 0.0;	 // You can Initialize the pressure to any value but should be equal throughout
		
		}
		
		else if (i == grid_points-1) // Bottom Wall
		{
			u_final[i][j] = 0.0;
			v_final[i][j] = 0.0;
			p_final[i][j] = 0.0;	
		
		}
		
		else if (j == 0)			// Left Wall
		{
			u_final[i][j] = 0.0;
			v_final[i][j] = 0.0;
			p_final[i][j] = 0.0;
			
		}
		
		else if (j == grid_points-1) // Right Wall
		{
			u_final[i][j] = 0.0;
			v_final[i][j] = 0.0;
			p_final[i][j] = 0.0;
	
		}
		
		else						// Interior Points
		{
			u_final[i][j] = 0.0;
			v_final[i][j] = 0.0;
			p_final[i][j] = 0.0;
		}
	
	}
	
}

// Staggered u-velocity
for (int i = 0; i<grid_points+1; i++)
{
	for(int j = 0; j<grid_points; j++)
	{   
		
		
		if (i == 0)
		{
			u[i][j] = 2.0;
			u_new[i][j] = 2.0;
		
		}
		
		else if (i == grid_points)
		{
			u[i][j] = 0.0;
			u_new[i][j] = 0.0;
		
		}
		else if (j == 0)
		{
			u[i][j] = 0.0;
			u_new[i][j] = 0.0;
			
		}
		
		else if (j == grid_points-1)
		{
			u[i][j] = 0.0;
			u_new[i][j] = 0.0;
	
		}
		else
		{
			u[i][j] = 0.0;
			u_new[i][j] = 0.0;;
		}
	
	}

}

// Staggered v-velocity
for (int i = 0; i<grid_points; i++)
{
	for(int j = 0; j<grid_points+1; j++)
	{   
		
		
		v[i][j] = 0.0;
		v_new[i][j] = 0.0;
	}

}

// Staggered pressure
for (int i = 0; i<grid_points+1; i++)
{
	for(int j = 0; j<grid_points+1; j++)
	{   
		
		
		p[i][j] = 0.0;
		p_new[i][j] = 0.0;
	}

}

// Gauss iteration starts
do 
{
	for (int i = 1; i<grid_points; i++)
    {	
		for (int j = 1; j<grid_points-1; j++)
        {    
			pressure = -(p[i][j+1] - p[i][j])/dx;
            diffusion = (1/Re)*((u[i+1][j] - 2*u[i][j] + u[i-1][j])/(dx*dx) + (u[i][j+1] - 2*u[i][j] + u[i][j-1])/(dx*dx));
            advection_x = (pow((0.5*(u[i][j]+u[i][j+1])),2) - pow((0.5*(u[i][j]+u[i][j-1])),2))/dx;
            advection_y = ((0.25*(u[i][j]+u[i-1][j])*(v[i-1][j]+v[i-1][j+1])) - (0.25*(u[i][j]+u[i+1][j])*(v[i][j]+v[i][j+1])))/dy;
            u_new[i][j] = u[i][j] + dt*(diffusion - advection_x - advection_y + pressure);
        }
    }

	// Boundary conditions for x-momentum
	for (int i = 0; i<grid_points+1; i++)
	{
		for(int j = 0; j<grid_points; j++)
		{   
			if(i == 0)								                        // Top Wall
			{
				u_new[i][j] = 2 - u_new[i+1][j];
			}
			
			if(i == grid_points)					                        // Bottom Wall
			{
				u_new[i][j] = -u_new[i-1][j];;
			}
			
			if((j == 0 || j == grid_points-1) && (i>0 && i<grid_points))	// Side Walls
			{
				u_new[i][j] = 0;;
			}
		}
	}
	
	// y - momentum equation 
	for (int i = 1; i<grid_points-1; i++)
    {	
		for (int j = 1; j<grid_points; j++)
        {    
			pressure = -(p[i][j] - p[i+1][j])/dy;
            diffusion = (1/Re)*((v[i+1][j] - 2*v[i][j] + v[i-1][j])/(dy*dy) + (v[i][j+1] - 2*v[i][j] + v[i][j-1])/(dx*dx));
            advection_y = (pow((0.5*(v[i][j]+v[i-1][j])),2) - pow((0.5*(v[i][j]+v[i+1][j])),2))/dy;
            advection_x = ((0.25*(u[i][j]+u[i+1][j])*(v[i][j]+v[i][j+1])) - (0.25*(u[i][j-1]+u[i+1][j-1])*(v[i][j]+v[i][j-1])))/dx;
            v_new[i][j] = v[i][j] + dt*(diffusion - advection_x - advection_y + pressure);
        }
    }
	
	// Boundary conditions
	for (int i = 0; i<grid_points; i++)
	{
		for(int j = 0; j<grid_points+1; j++)
		{   
			if(i == 0 || i == grid_points - 1)				   // Top and Bottom walls
			{
				v_new[i][j] = 0;
			}
			
			if(j == 0 && i>0 && i<grid_points-1 )			  // Left Wall
			{		
				v_new[i][j] = -v_new[i][j+1];;
			}
			
			if(j == grid_points && i>0 && i<grid_points-1 )   // Right Wall
			{
				v_new[i][j] = -v_new[i][j-1];;
			}
			
			
		}
	}
	
	// Continuity Equation
	for (int i = 1; i<grid_points; i++)
	{
		for(int j = 1; j<grid_points; j++)
		{   
			p_new[i][j] = p[i][j] - delta*dt*(u_new[i][j] - u_new[i][j-1] + v_new[i-1][j] - v_new[i][j])/dx;
			
		}
	}
	

	for (int i = 0; i<grid_points; i++)
	{
		for(int j = 0; j<grid_points+1; j++)
		{   
			if(i == 0)
			{
				p_new[i][j] = p_new[i+1][j];
			}
			
			if(i == grid_points)
			{
				p_new[i][j] = p_new[i-1][j];;
			}
			
			if(j == 0 && i>0 && i<grid_points )
			{
				p_new[i][j] = p_new[i][j+1];;
			}
			
			if(j == grid_points && i>0 && i<grid_points )
			{
				p_new[i][j] = p_new[i][j-1];;
			}
			
			
		}
	}
	
	error_mag = 0;
    for (int i = 1; i<grid_points; i++)
	{
		for(int j = 1; j<grid_points; j++)
		{   
            error_mag = error_mag + abs((u_new[i][j] - u_new[i][j-1] + v_new[i-1][j] - v_new[i][j])/dx);
        }
    }
    
	if(iterations%1000 == 0 || iterations == 1)
	{	
		cout << "Residual error after " << iterations <<" iterations is" <<" " << error_mag << endl;
	}
	iterations+=1;
	for (int i = 0; i<grid_points; i++)
	{
		for(int j = 0; j<grid_points; j++)
		{
			u[i][j] = u_new[i][j];
			v[i][j] = v_new[i][j];
			p[i][j] = p_new[i][j];
	
		}

	}

}while(error_mag > error_req);

// Reverting back to the collocated grid
for (int i  = 0; i<grid_points-1; i++)
{    
	for (int j = 0; j<grid_points; j++)
    {    
		u_final[i][j] = 0.5*(u[i][j] + u[i+1][j]);
        v_final[i][j] = 0.5*(v[i][j] + v[i][j+1]);
        p_final[i][j] = 0.25*(p[i][j] + p[i][j+1] + p[i+1][j] + p[i+1][j+1]);
    }
}

// Creating Meshgrids
for(int i = 0; i<grid_points;i++)
{
	y[i] = i*dy;
	x[i] = i*dx;
}

// Output of Midline U velocity
ofstream outfile;
outfile.open("MidVelocities.dat");
if(outfile.is_open())
{
    for(int i = 0; i<grid_points; i++)
	{
		outfile << u_final[i][(grid_points+1)/2] << " " << v_final[(grid_points+1)/2][i] <<  " "  << x[i] << " " << y[grid_points - (i+1) ];
		outfile << endl;
    }
    outfile.close();
}
else cout << "ERROR";

}