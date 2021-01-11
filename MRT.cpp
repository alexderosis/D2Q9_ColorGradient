#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <cmath>
#include <cassert>
#include <vector>
#include <iostream>
#include <fstream>
#include <sstream>
#include <ctime>
using namespace std;
///----------------------------------------------------------------------------------------------------------------------------------
///----------------------------------------------------------------------------------------------------------------------------------
const bool plot_vtk = true;
const int n_phase = 2;
const int nx = 500, ny = 2*nx, np = 9;
const double NX = (double)nx-1, gravity = 0.05*0.05/NX, rho0_b = 1., At = 0.85, rho0_r = -rho0_b*(At+1)/(At-1), Reynolds = 300000., nu = sqrt(NX*gravity)*NX/Reynolds, T = sqrt(NX/gravity);
const double cs2 = 1./3., beta = 0.7, alpha_b = 4./9., alpha_r = 1.-(1.-alpha_b)*rho0_b/rho0_r;
vector<const int> cx = {0, 1, 0, -1, 0, 1, -1, -1, 1},
									cy = {0, 0, 1, 0, -1, 1, 1, -1, -1},
									opp= {0, 3, 4, 1, 2, 7, 8, 5, 6};
vector<const double> wf = {4./9., 1./9., 1./9., 1./9., 1./9., 1./36., 1./36., 1./36., 1./36.};
vector<const double> B = {-4/27., 2/27., 2/27., 2/27., 2/27., 5/108., 5/108., 5/108., 5/108.};
vector<const double> alphaK = {alpha_b, alpha_r}, ni = {nu, nu};
const int nsteps = (int)(3*T+1), n_out = (int)(T/100);
vector<double> f(nx*ny*np,0.), f_old(nx*ny*np,0.), rhoK(nx*ny*n_phase,0.), rhoK_old(nx*ny*n_phase,0.), rho(nx*ny,0.), u(nx*ny,0.), v(nx*ny,0.), rho_old(nx*ny,0.);
vector<double> gradx_rhoK(nx*ny*n_phase,0.), grady_rhoK(nx*ny*n_phase,0.);
vector<double> gradx_rho(nx*ny,0.), grady_rho(nx*ny,0.);
double U, V, W, ftemp, A, C, R;
double colorGradX, colorGradY, colorGradNorm, tempx, tempy, norm_c, tmp, tmp3, phi;
double alpha, totalMass, streamed;
vector<double> kTotalMass(nx*ny,0.);
const double surfaceTension = 1E-5;
double A_, R1, nu_eff, omega_eff, GX, GY, Ckl;
int check;
///-CMS
vector<double> temp_pop(np,0.);
double k0, k1, k2, k3, k4, k5, k6, k7, k8;
double r0, r1, r2, r3, r4, r5, r6, r7, r8;
double CX, CY, U2, V2, UV, U3, V3;
double value;
int id, idn, newx, newy, id1, id2;
double FX, FY;
///----------------------------------------------------------------------------------------------------------------------------------
///----------------------------------------------------------------------------------------------------------------------------------
void write_fluid_vtk(int time)
{
	stringstream output_filename;
	output_filename << "vtk_fluid_MRT/fluid_t" << time << ".vtk";
	ofstream output_file;
	output_file.open(output_filename.str().c_str());

	output_file << "# vtk DataFile Version 3.0\n";
	output_file << "fluid_state\n";
	output_file << "ASCII\n";
	output_file << "DATASET RECTILINEAR_GRID\n";
	output_file << "DIMENSIONS " << nx << " " << ny << " " << "1" << " " << "\n";
	output_file << "X_COORDINATES " << nx << " float\n";
	for(int x = 0; x < nx; ++x)
		output_file << x << " ";
	output_file << "\n";
	output_file << "Y_COORDINATES " << ny  << " float\n";
	for(int y = 0; y < ny ; ++y)
		output_file << y  << " ";
	output_file << "\n";
	output_file << "Z_COORDINATES " << 1  << " float\n";
	output_file << 0  << " ";
	output_file << "\n";
	output_file << "POINT_DATA " << (nx) * (ny) * 1  << "\n";

	output_file << "SCALARS density float 1\n";
	output_file << "LOOKUP_TABLE default\n";
	for(int Y = 0; Y < ny ; ++Y)
		for(int X = 0; X < nx; ++X)
    {
      id = X*ny+Y;
			output_file << rho[id] << "\n";
    }
	output_file << "SCALARS order_param float 1\n";
	output_file << "LOOKUP_TABLE default\n";
	for(int Y = 0; Y < ny ; ++Y)
		for(int X = 0; X < nx; ++X)
    {
      id = X*ny+Y;
      output_file << (rhoK[id*n_phase+1]-rhoK[id*n_phase+0])/(rhoK[id*n_phase+0]+rhoK[id*n_phase+1]) << "\n";
    }

	output_file << "VECTORS velocity_vector float\n";
	for(int Y = 0; Y < ny ; ++Y)
		for(int X = 0; X < nx; ++X)
    {
      id = X*ny+Y;
			output_file << u[id] << " " << v[id] << " " << "0" << "\n";
    }

	output_file.close();
}
///----------------------------------------------------------------------------------------------------------------------------------
void initial_state()
{
	double phi, X;
	int max = 1, min = -1;
	double h, an, bn;
	for(int x=0; x<nx; x++)
    for(int y=0; y<ny; y++)
		{
      id = x*ny+y;
			U = u[id] = 0.;
      V = v[id] = 0.;
      X = (double)x / ((double)nx-1);
			h = 0;
			for(int n=30; n<40; n++)
			{
				an = rand()%(max-min+1)+min;
				bn = rand()%(max-min+1)+min;
				h += an*cos(2.*M_PI*n*X)+bn*sin(2.*M_PI*n*X);
			}
			h = 0.5*ny+nx*0.002*h;
      if(y>h)
      {
	      rhoK[id*n_phase+0] = 0.;
        rhoK[id*n_phase+1] = rho0_r;
	    }
  	  else
  	  {
  	    rhoK[id*n_phase+0] = rho0_b;
  	    rhoK[id*n_phase+1] = 0.;
  	  }
			rho[id] = 0.;
			for(int k=0; k<n_phase; k++)
				rho[id] += rhoK[id*n_phase+k];
		  alpha = 0.;
			for(int k=0; k<n_phase; k++)
				alpha += rhoK[id*n_phase+k]*alphaK[k];
			alpha /= rho[id];
			C = -1.5*(U*U+V*V);
			for(int n=0; n<np; n++)
			{
        A = U*cx[n]+V*cy[n];
        if(n==0)
          phi = alpha;
        else if(n<5)
          phi = (1.-alpha)/5.;
        else
          phi = (1.-alpha)/20.;
				f[id*np+n] = rho[id]*(phi+wf[n]*(3.*A+4.5*A*A+C));
		  }
		}
}
///----------------------------------------------------------------------------------------------------------------------------------
void compute_gradient_rho(int x, int y)
{
  id = x*ny+y;
  for(int k=0; k<n_phase; k++)
  {
  	gradx_rhoK[id*n_phase+k] = 0.;
    grady_rhoK[id*n_phase+k] = 0.;
		if(k==0)
		{
			gradx_rho[id] = 0.;
			grady_rho[id] = 0.;
		}
		if(y>0 && y<ny-1)
	    for(int n=1; n<np; n++)
		  {
			  newx = x+cx[n];
			  newy = y+cy[n];
			  if(x==0 || x==nx-1)
				  newx = (newx+nx)%nx;
			  if(y==0 || y==ny-1)
				 newy = (newy+ny)%ny;
        idn = newx*ny+newy;
	  		value = rhoK[idn*n_phase+k]/rho_old[idn];
			  gradx_rhoK[id*n_phase+k] += 3.*wf[n]*cx[n]*value;
			  grady_rhoK[id*n_phase+k] += 3.*wf[n]*cy[n]*value;
				if(k==0)
				{
					gradx_rho[id] += 3.*wf[n]*cx[n]*rho[idn];
					grady_rho[id] += 3.*wf[n]*cy[n]*rho[idn];
				}
			}
		if(y==0) // SOUTH wall
		{
      id1 = ((x+1+nx)%nx)*ny+y;
      id2 = ((x-1+nx)%nx)*ny+y;
			gradx_rhoK[id*n_phase+k] = 0.5*(rhoK[id1*n_phase+k]/rho_old[id1] - rhoK[id2*n_phase+k]/rho_old[id2]);
			if(x==0)
				gradx_rhoK[id*n_phase+k] = -rhoK[id*n_phase+k]/rho_old[id] + rhoK[id1*n_phase+k]/rho_old[id1];
			if(x==nx-1)
				gradx_rhoK[id*n_phase+k] = rhoK[id*n_phase+k]/rho_old[id] - rhoK[id2*n_phase+k]/rho_old[id2];

      id1 = x*ny+(y+1);
			id2 = x*ny+(y+2);
			grady_rhoK[id*n_phase+k] = -3./2*rhoK[id*n_phase+k]/rho_old[id] + 2.*rhoK[id1*n_phase+k]/rho_old[id1] - 1./2.*rhoK[id2*n_phase+k]/rho_old[id2];
		}
		if(y==ny-1) // NORTH wall
		{
      id1 = ((x+1+nx)%nx)*ny+y;
      id2 = ((x-1+nx)%nx)*ny+y;
			gradx_rhoK[id*n_phase+k] = 0.5*(rhoK[id1*n_phase+k]/rho_old[id1] - rhoK[id2*n_phase+k]/rho_old[id2]);
			if(x==0)
				gradx_rhoK[id*n_phase+k] = -rhoK[id*n_phase+k]/rho_old[id] + rhoK[id1*n_phase+k]/rho_old[id1];
			if(x==nx-1)
				gradx_rhoK[id*n_phase+k] = rhoK[id*n_phase+k]/rho_old[id] - rhoK[id2*n_phase+k]/rho_old[id2];

      id1 = x*ny+(y-1);
			id2 = x*ny+(y-2);
			grady_rhoK[id*n_phase+k] = 3./2*rhoK[id*n_phase+k]/rho_old[id] - 2.*rhoK[id1*n_phase+k]/rho_old[id1] + 1./2.*rhoK[id2*n_phase+k]/rho_old[id2];
		}
	}
}
///----------------------------------------------------------------------------------------------------------------------------------
void perturbation(int x, int y, double omega)
{
  id = x*ny+y;
	A_ = 9*surfaceTension*omega*0.25;
	R1 = 1./rho[id];
  for(int k=0; k<n_phase; k++)
    for(int l=0; l<n_phase; l++)
    	if(k!=l)
    	{
				colorGradX = R1*(gradx_rhoK[id*n_phase+k]*rhoK[id*n_phase+l]-
                 		     gradx_rhoK[id*n_phase+l]*rhoK[id*n_phase+k]);
    		colorGradY = R1*(grady_rhoK[id*n_phase+k]*rhoK[id*n_phase+l]-
          		        	 grady_rhoK[id*n_phase+l]*rhoK[id*n_phase+k]);
    		colorGradNorm = max(sqrt(colorGradX*colorGradX+
  	          							     colorGradY*colorGradY),1e-12);
    		Ckl = min(1E6*rhoK[id*n_phase+l]*rhoK[id*n_phase+k]/rho0_r/rho0_b,1.0);
    		for(int n=0; n<np; n++)
    		{
    			tmp = colorGradX*cx[n]+colorGradY*cy[n];
      		tmp3 = 0.5*A_*colorGradNorm*Ckl*(wf[n]*tmp*tmp/colorGradNorm/colorGradNorm-B[n]);
      		f[id*np+n] += tmp3;
    		}
		  }
}
///----------------------------------------------------------------------------------------------------------------------------------
void recoloring_and_streaming()
{
	f_old = f;
  rhoK_old = rhoK;
  std::fill(f.begin(), f.end(), 0.);
	for(int k=0; k<n_phase; k++)
  {
		std::fill(kTotalMass.begin(), kTotalMass.end(), 0.);
    for(int x=0; x<nx; x++)
			for(int y=0; y<ny; y++)
			{
        id = x*ny+y;
				R = rho[id];
				R1 = 1./R;
  	    tempx = tempy = 0.;
  	    for(int l=0; l<n_phase; l++)
  	    	if(k!=l)
  	      {
  	      	colorGradX = gradx_rhoK[id*n_phase+k]*rhoK_old[id*n_phase+l]-
              	     		     gradx_rhoK[id*n_phase+l]*rhoK_old[id*n_phase+k];
    				colorGradY = grady_rhoK[id*n_phase+k]*rhoK_old[id*n_phase+l]-
          		          		 grady_rhoK[id*n_phase+l]*rhoK_old[id*n_phase+k];
    				colorGradNorm = max(sqrt(colorGradX*colorGradX+
  	          	  						         colorGradY*colorGradY),1e-12);
  	        tempx += beta*colorGradX/colorGradNorm*rhoK_old[id*n_phase+l];
  	        tempy += beta*colorGradY/colorGradNorm*rhoK_old[id*n_phase+l];
					}
				alpha = 0.;
				for(int kk=0; kk<n_phase; kk++)
					alpha += rhoK_old[id*n_phase+kk]*alphaK[kk];
				alpha *= R1;
				for(int n=0; n<np; n++)
				{
  	      norm_c = sqrt(cx[n]*cx[n]+cy[n]*cy[n]);
  	      if(n==0)
  	      	norm_c = 0.;
  	      else
  	      	norm_c = 1./norm_c;
  	      tmp = tempx*cx[n]*norm_c+tempy*cy[n]*norm_c;
          if(n==0)
  	        phi = alpha;
  	      else if(n<5)
  	        phi = (1.-alpha)/5.;
        	else
    	    	phi = (1.-alpha)/20.;
  	      tmp *= phi;
  	      newx = x+cx[n];
  	      newy = y+cy[n];
  	      if(x==0 || x==nx-1)
  	      	newx = (newx+nx)%nx;
  	      if(y==0 || y==ny-1)
  	        newy = (newy+ny)%ny;
  	      idn = newx*ny+newy;
        	streamed = (f_old[id*np+n]+tmp)*rhoK_old[id*n_phase+k]*R1;
					f[idn*np+n] += streamed;
					kTotalMass[idn] += streamed;
					rhoK[idn*n_phase+k] = kTotalMass[idn];
				}
			}
  }
}
///----------------------------------------------------------------------------------------------------------------------------------
int algorithm_MRT()
{
	rho_old = rho;
	check = 0;
	totalMass = 0.;
	/// 1. Compute color-blind macroscopic variables
	for(int x=0; x<nx; x++)
		for(int y=0; y<ny; y++)
		{
      id = x*ny+y;
			R = U = V = 0.;
			FX = FY = 0.;
	    for(int n=0; n<np; n++)
	    {
	    	ftemp = temp_pop[n] = f[id*np+n];
	      R += ftemp;
	      U += ftemp*cx[n];
	      V += ftemp*cy[n];
	    }
			R1 = 1./R;
		  FY = -(R - 0.5*(rho0_r+rho0_b))*gravity;
	    U *= R1;
	    V = (V + 0.5*FY)*R1;
	    rho[id] = R;
	    u[id] = U;
			v[id] = V;
			totalMass += R;
			compute_gradient_rho(x, y);
			GX = gradx_rho[id];
			GY = grady_rho[id];
			alpha = nu_eff = 0.;
			for(int k=0; k<n_phase; k++)
			{
				 alpha += rhoK[id*n_phase+k]*alphaK[k];
				 nu_eff += rhoK[id*n_phase+k]*ni[k];
			}
			alpha *= R1;
			nu_eff *= R1;
			omega_eff = 1./(3.*nu_eff+0.5);
			if(fabs(U)>100.)
        check = 1;
      ///compute moments
			U2 = U*U;
			V2 = V*V;
			UV = U*V;
			U3 = U2*U;
			V3 = V2*V;
			/*k4 = k5 = 0;
		  for(int n=0; n<np; n++)
      {
        CX = (double)cx[n];
        CY = (double)cy[n];
        ftemp = f[id*np+n];
        k4 += ftemp*(CX*CX-CY*CY);
				k5 += ftemp*CX*CY;
			}*/
			r4 = f[id*np+1]-f[id*np+2]+f[id*np+3]-f[id*np+4];
			r5 = f[id*np+5]-f[id*np+6]+f[id*np+7]-f[id*np+8];
      FY = -(R-0.5*(rho0_r+rho0_b))*gravity;
			r0 = R;
			r1 = R*U+0.5*FX/R;
			r2 = R*V+0.5*FY/R;
			r3 = (FX*U+FY*V)/R+R*(5.*U2+5.*V2-6.*alpha+6.)/5.;
			r4 = r4*(1.-omega_eff)+R*omega_eff*(U2-V2)-2.*(0.5*omega_eff-1.)*(FX*U-FY*V)/R;
			r5 = r5*(1.-omega_eff)+R*U*V*omega_eff-(0.5*omega_eff-1.)*(FY*U+FX*V)/R;
			r6 = 0.5*cs2*(3.*FY*U2+6.*FX*V*U+FY)/R+R*V*(3.*U2+1.)*cs2;
			r7 = 0.5*cs2*(3.*FX*V2+6.*FY*U*V+FX)/R+R*U*(3.*V2+1.)*cs2;
			r8 = R*(U2*V2+cs2*U2+cs2*V2-alpha/5.+1./5.)+cs2*(FY*U2*V+FX*U*V2+cs2*FX*U+cs2*FY*V)/R;

			//r0 += 3.*nu_eff*(U*GY+V*GX);
			//r3 += 4.*nu_eff*(GX+GY)*(U+V);
			//r4 += 2.*nu_eff*omega_eff*(U*GX-V*GY);
      //r5 += nu_eff*omega_eff*(U+V)*(GX+GY);
 		  //r8 += nu_eff*cs2*(GX*(4.*U+3.*V)+GY*(3.*U+4.*V));

			f[id*np+0] = r0-r3+r8;
      f[id*np+1] = 0.5*(r1-r7-r8)+0.25*(r3+r4);
      f[id*np+2] = 0.5*(r2-r6-r8)+0.25*(r3-r4);
      f[id*np+3] = 0.25*(r3+r4)+0.5*(-r1+r7-r8);
      f[id*np+4] = 0.25*(r3-r4)+0.5*(-r2+r6-r8);
      f[id*np+5] = 0.25*(r5+r6+r7+r8);
      f[id*np+6] = 0.25*(r6-r5-r7+r8);
      f[id*np+7] = 0.25*(r5-r6-r7+r8);
      f[id*np+8] = 0.25*(r7-r6-r5+r8);
			perturbation(x, y, omega_eff);
		}
	recoloring_and_streaming();
	return check;
}
///----------------------------------------------------------------------------------------------------------------------------------
void boundary()
{
	for(int i=0; i<nx; i++)
		for(int n=0; n<np; n++)
		{
      if(cy[n]>0)
			{
  			id = i*ny+0;
  			f[id*np+n] = f_old[id*np+opp[n]];
  		}
  		if(cy[n]<0)
  		{
  			id = i*ny+ny-1;
  			f[id*np+n] = f_old[id*np+opp[n]];
  		}
		}
}
///----------------------------------------------------------------------------------------------------------------------------------
int main(int argc, char *argv[])
{
	FILE *data_output;
	data_output = fopen("data_MRT.txt","wt");
	system("mkdir vtk_fluid_MRT");
	initial_state();
	int check_mach = 0;
	printf("%lf %lf\n", rho0_b , rho0_r);
	for(int i=0; i<nsteps; i++)
  {
    check_mach = algorithm_MRT();
    boundary();
    if(check_mach==1)
      goto labelA;
		if(plot_vtk==true && i%n_out==0)
			write_fluid_vtk(i);
    if(i%1==0)
      printf("Iteration %d of %d. TotMass=%e\n", i, nsteps, totalMass);
    fprintf(data_output,"%d    %e\n", i, totalMass);
  }
  labelA:
  fclose(data_output);
  return 0;
}
///----------------------------------------------------------------------------------------------------------------------------------
///----------------------------------------------------------------------------------------------------------------------------------
