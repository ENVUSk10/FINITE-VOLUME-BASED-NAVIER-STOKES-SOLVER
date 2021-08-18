// C++ Code For Lid driven Cavity using Stream-Vorticity formulation
#include<bits/stdc++.h>
using namespace std;
int main()
{
double domain_size = 1.0, grids = 51.0,U = 1.0;
double h = domain_size/(grids-1);
double nu = 1e-2;
double Beta = 1.95;
double dt = 0.001;
double CFL = (U*dt)/h;
double Re = (U*domain_size)/nu;
double F = (nu*dt)/pow(h,2);
int grid_points = 51;

double error_mag = 1,error_req = 1e-8;
int iterations = 1;

double v[grid_points][grid_points], u[grid_points][grid_points], psi[grid_points][grid_points], w[grid_points][grid_points];
double w_old[grid_points][grid_points],y[grid_points],x[grid_points];

// CFL number check
if (pow(CFL,2) <=2*F && CFL <= 1 && F<=0.5)
    cout << "Stable Condition" << endl;
else 
   cout << "Unstable Condition" << endl;

// Meshpoints
for(int i = 0; i<grid_points;i++)
{
	y[i] = i*h;
	x[i] = i*h;
}


// Intializing the variables
for (int i = 0; i<grid_points; i++)
{
	for(int j = 0; j<grid_points; j++)
	{
		if (i == 0) //top wall
		{
			u[i][j] = 1.0; 
			v[i][j] = 0.0;
			psi[i][j] = 0.0;
			w[i][j] = 0.0;
			w_old[i][j] = 0.0;
		
		}
		
		else if (i == grid_points-1) // bottom wall
		{
			u[i][j] = 0.0;
			v[i][j] = 0.0;
			psi[i][j] = 0.0;
			w[i][j] = 0.0;
			w_old[i][j] = 0.0;
		
		}
		
		else if (j == 0) // Left wall
		{
			u[i][j] = 0.0;
			v[i][j] = 0.0;
			psi[i][j] = 0.0;
			w[i][j] = 0.0;
			w_old[i][j] = 0.0;
			
		}
		
		else if (j == grid_points-1) // right wall
		{
			u[i][j] = 0.0;
			v[i][j] = 0.0;
			psi[i][j] = 0.0;
			w[i][j] = 0.0;
			w_old[i][j] = 0.0;
	
		}
		else
		{
			u[i][j] = 0.0;
			v[i][j] = 0.0;
			psi[i][j] = 0.0;
			w[i][j] = 0.0;
			w_old[i][j] = 0.0;
		}
	
	}

}





// Numerical Solution
while(error_mag > error_req)
{
	// Stream Function
	for (int i = 1; i<grid_points-1; i++)
	{
		for(int j = 1; j<grid_points-1; j++)
		{
			psi[i][j] = 0.25*Beta*(psi[i+1][j] + psi[i-1][j] +psi[i][j+1] +  psi[i][j-1] + pow(h,2)*w[i][j]) + (1 - Beta)*(psi[i][j]);
	
		}

	}
	
	// Boundary Condition for Vorticity
	for (int i = 0; i<grid_points; i++)
	{
		for(int j = 0; j<grid_points; j++)
		{
			w[0][j] = -2.0*psi[1][j]/pow(h,2)- ((2*U)/h);					// Top Wall
            w[i][0] = -2.0*psi[i][1]/pow(h,2);								// Left Wall
            w[i][grid_points-1] = -2.0*psi[i][grid_points-2]/pow(h,2) ;		// Right Wall
            w[grid_points-1][j] = -2.0*psi[grid_points-2][j]/pow(h,2);		// Bottom Wall
	
		}

	}
	
	// Vorticity and Error Calculation
	error_mag = 0;	
	for (int i = 1; i<grid_points-1; i++)
	{
		for(int j = 1; j<grid_points-1; j++)
		{
			w[i][j] = w[i][j] - (dt/(4*pow(h,2)))*((psi[i][j+1] - psi[i][j-1])*(w[i+1][j] - w[i-1][j])) + (dt/(4*pow(h,2)))*((psi[i+1][j] 
					- psi[i-1][j])*(w[i][j+1] - w[i][j-1])) + (dt/(Re*pow(h,2)))*(w[i+1][j] + w[i-1][j] 
                     + w[i][j+1] + w[i][j-1] -4*w[i][j]);
            
             error_mag = error_mag + abs(w[i][j] - w_old[i][j]); 
	
		}

	}
	
	
	if(iterations%1000 == 0 || iterations == 1)
	{	
		cout << "Error after " << iterations << " " << error_mag << endl;
	}
	
	// Update the old Value of Vorticity
	for (int i = 0; i<grid_points; i++)
	{
		for(int j = 0; j<grid_points; j++)
		{
			w_old[i][j] = w[i][j];
	
		}

	}
	iterations += 1;
}

// Velocity Calculation
for (int i = 1; i<grid_points-1; i++)
	{
		for(int j = 1; j<grid_points-1; j++)
		{
			u[i][j] = (0.5/h)*-(psi[i+1][j] - psi[i-1][j]);
			v[i][j] = (0.5/h)*(psi[i][j+1]) - psi[i][j-1];
			
		}

	}

// Output midline U-velocities for Validation with the Benchmark
ofstream outfile;
outfile.open("Midvelocity.csv");
if(outfile.is_open()){
    for(int i = 0; i<grid_points; i++){
        
        outfile << y[grid_points - (i+1) ] << "," <<  u[i][(grid_points+1)/2] ;
        outfile << endl;
    }
    outfile.close();
}
else cout << "ERROR";

}
