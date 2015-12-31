#include <iostream>
#include <cmath>
#include <time.h>
#include <algorithm>
#include <vector>
#include <ctime>
#include <cstdio>
#include <cstdlib>
#include <cassert>
#include <string>
#include "miniapp.hh"
#include "auxiliary_function.hh"

using namespace std;

typedef std::vector<real> v_dbl;
#ifdef FLOAT
typedef float real;
const int factor = 16;
#else
typedef double real;
const int factor = 8;
#endif


//----------------------------------------------------------------------------//
Solver::Solver(int n_eg_in, int n_a_in, int cm_xy_in, int fm_xy_in,
		int cm_z_in, int fm_z_in, int upscatter_in, int iter_in,
		int xbs, int ybs, int zbs, int nTs_in)
: n_eg(n_eg_in)
, n_a(n_a_in)
, upscatter(upscatter_in)
, iter(iter_in)
, blocksize_z(zbs)
, nTs(nTs_in)
, mesh(cm_xy_in, fm_xy_in, cm_z_in, fm_z_in, xbs, ybs, zbs, nTs)
{

	cout << "input parameters \n" << "# of energy groups: " << n_eg << "\n"
			<< "# of angles in each octant: " << n_a << "\n"
			<< "# of coarse cells along x & y axes: " << cm_xy_in << " \n"
			<< "# of fine cells in one coarse cell along x & y axes: "<< fm_xy_in << "\n"
			<< "# of coarse cells along z axis: " << cm_z_in << " \n"
			<< "# of fine cells in one coarse cell along z axis: "<< fm_z_in << "\n"
			<< "# of upscatter energy groups: " << upscatter << "\n"
			<< "# of running iterations: " << iter << "\n"
			<< "x block size: " << xbs << endl
			<< "y block size: " << ybs << endl
			<< "z block size: " << zbs << endl
      << "number threads: " << nTs << endl;

	//mesh.print_summary();
	phi_size = mesh.ncell * n_eg;
	phi.resize(phi_size);
	Q.resize(phi_size);

	//allocate the total cross section point and assign all to 0.5
	SigT.resize(mesh.nmat * n_eg, 0.5);

	// allocate the scattering cross section to a low tridiagonal matrix and
	// assign the value 0.2
	SigS = new real** [mesh.nmat];
	for(int i = 0; i < mesh.nmat; i++)
	{
		SigS[i] = new real* [n_eg];
		//for energy groups without upscattering, only energy groups with energy >=
		//(group ID <=) them can scatter into them
		for (int j = 0; j < n_eg - upscatter; j++)
			SigS[i][j] = new real [j + 1];
		//for energy groups with upscattering
		for(int j = n_eg - upscatter; j < n_eg; j++)
			SigS[i][j] = new real [n_eg];
	}

	for(int i = 0; i < mesh.nmat; i++)
	{
		for(int j = 0; j < n_eg - upscatter; j++)
			for(int k = 0; k < j + 1; k++)
				SigS[i][j][k] = 0.2;
		for(int j = n_eg - upscatter; j < n_eg; j++)
			for(int k = 0; k < n_eg; k++)
				SigS[i][j][k] = 0.2;
	}

	//Initialize mu,xi,eta;
	get_quadrature();
	print_line();
	cout << "initialization finished \n";
	print_line();
}

Solver::~Solver()
{
	if(SigS)
	{
		for(int i = 0; i < mesh.nmat; i++)
		{
			for(int j = 0; j < n_eg; j++)
				delete[] SigS[i][j];
			delete[] SigS[i];
		}
		delete[] SigS;
	}
}

//----------------------------------------------------------------------------//
void Solver::get_quadrature()
{
	/*  This provides an easy but sort of realistic quadrature defined only
	 *  by the numer of angles per octant.
	 */
	mu.resize(n_a);
	eta.resize(n_a);
	xi.resize(n_a);
	wt.resize(n_a);
	for (int i = 0; i < n_a; ++i)
	{
		mu[i] = 1.0/(i+0.5/i/n_a);
		int j = n_a - i - 1;
		eta[i] = 1.0/(j+0.5/j/n_a);
		xi[i] = std::sqrt(1.0-mu[i]*mu[i]-eta[i]*eta[i]);
		wt[i] = 1.5707963267948966 / n_a;
	}
}

void Solver::Calculate(std::string order)
{
	//nTs = num_threads;
	cout << "# of threads is " << nTs << endl;
	// initiate values for old sweep
	N = sqrt(nTs);
	blocksize_xy = mesh.nx / N;
	n_b = mesh.nz / blocksize_z;
	n_p = 2 * N - 2 + n_b;
	remain = mesh.nz % blocksize_z;

	// make sure nTs is square num
	std:vector<int> square_num(10);
	for (int i = 0; i < square_num.size(); i++)
		square_num[i] = (i + 1) * (i + 1);
	if (order == "aes" or order == "ase" or order == "eas" or
			order == "esa" or order == "sae" or order == "sea")
	{
		assert(std::find(square_num.begin(), square_num.end(), nTs) != square_num.end());
		assert(mesh.nx % N == 0);
	}


	int start_TID[(2 * N - 1) * 2];// start_TID -- TID info for planes having starting threads;
	Set_start_TID(start_TID, N);

	time_used_total = omp_get_wtime();
	assert(order == "aes" or order == "ase" or order == "eas" or
			order == "esa" or order == "sae" or order == "sea" or
			order == "aes_mod" or order == "ase_mod" or order == "eas_mod" or
			order == "esa_mod" or order == "sae_mod" or order == "sea_mod");

	// initial guess of phi, 0.0
	std::fill(phi.begin(), phi.end(), 0.0);

	time_used_sweep = 0.0;
	//convergence parameters
	real eps_phi = 1.0e-5;
	real err_phi = 1.0;
	int it = 0;

	while (it < iter)
	{
		// store old phi
		v_dbl phi0(phi);

		// calculate Q
		for (int k = 0; k < mesh.nz; k++)
			for (int j = 0; j < mesh.ny; j++)
				for (int i = 0; i < mesh.nx; i++)
				{
					int cell = mesh.ijk_to_cell(i, j, k);
					int mat = mesh.mat_id[cell];
					// non-upscatter E range
					for(int e = 0; e < n_eg - upscatter; ++e)
					{
						real totsca = 0.0;
						for(int ee = 0; ee < e + 1; ee++)
							totsca += phi[cell * n_eg + ee] * SigS[mat][e][ee];
						// using pi/4, not S[i][j][k][e], since it is
						// uniform and isotropic source.
						Q[cell * n_eg + e] = 1.0 / (M_PI * 4.0) +
								totsca / (M_PI * 4.0);
					}
					// upscatter E range
					for(int e = n_eg - upscatter; e < n_eg; e++)
					{
						real totsca = 0.0;
						for(int ee = 0; ee < n_eg; ee++)
							totsca += phi[cell * n_eg + ee] * SigS[mat][e][ee];
						Q[cell * n_eg + e] = 1.0 / (M_PI * 4.0) +
								totsca / (M_PI * 4.0);
					}
				}


		double sweep_time = omp_get_wtime();

		if (order == "aes")
			sweep_aes(start_TID);
		else if (order == "ase")
			sweep_ase(start_TID);
		else if (order == "eas")
			sweep_eas(start_TID);
		else if (order == "esa")
			sweep_esa(start_TID);
		else if(order == "sae")
			sweep_sae(start_TID);
		else if(order == "sea")
			sweep_sea(start_TID);

		else if (order == "aes_mod")
			sweep_aes_mod();
		else if (order == "ase_mod")
			sweep_ase_mod();
		else if (order == "eas_mod")
			sweep_eas_mod();
		else if (order == "esa_mod")
			sweep_esa_mod();
		else if (order == "sae_mod")
			sweep_sae_mod();
		else if (order == "sea_mod")
			sweep_sea_mod();


		sweep_time = omp_get_wtime() - sweep_time;
		time_used_sweep += sweep_time;

		real max = fabs(phi[0] - phi0[0]) / phi0[0];
		for (int i = 1; i < phi_size; i++)
			if (fabs(phi[i] - phi0[i]) / phi0[i] > max)
				max = fabs(phi[i] - phi0[i]) / phi0[i];
		err_phi = max;
		cout << "iteration " << it << " error of phi is " << err_phi << endl;
		it = it + 1;
	}

	time_used_total = omp_get_wtime() - time_used_total;
	print_line();
	cout << "summary\n";
	print_line();

	if (it <= iter && err_phi < eps_phi)
	{
		cout << "\nconverged in " << it << " iterations" << "\n";
		cout << "time used in " << order << " sweep of " << it << " iterations is " << time_used_sweep << " s" << endl;
		cout << "elapsed time per sweep is " << time_used_sweep / it << " s" << endl;
		cout << "total time is " << time_used_total << " s" << endl;

		if (phi_size < 10){
			cout << "\nthe resulting phi is \n";
			for(int i = 0; i < phi_size; i++)
				cout << phi[i] << " ";
			cout << '\n';}
		else{
			cout << "\nthe first 10 resulting phi is \n";
			for(int i = 0; i < 10; i++)
				cout << phi[i] << " ";
			cout << endl;}
		cout << "\nerr_phi is " << err_phi << "\n\n";
	}

	else
	{
		cout << "\ndo not converge and time used in "<< order << " sweep in " << iter << " iterations is " << time_used_sweep << " s" << endl;
		cout << "elapsed time per sweep is " << time_used_sweep / iter << " s" << endl;
		cout << "total time used is " << time_used_total << " s" << endl;

		if (phi_size < 10){
			cout << "\nthe resulting phi is \n";
			for (int i = 0; i < phi_size; i++)
				cout << phi[i] << " ";
			cout << endl;}
		else{
			cout << "\nthe first 10 resulting phi is \n";
			for(int i = 0; i < 10; i++)
				cout << phi[i] << " ";
			cout << endl;}
		cout << "\nerr_phi is " << err_phi << "\n\n";
	}

	double flop = mesh.nx*mesh.ny*mesh.nz*n_a*8.*n_eg*iter*25;
	double gflops = flop/1.e9/time_used_sweep;
    printf("TIME    = %f \n", time_used_sweep);
    printf("FLOPS   = %f \n", flop);
	printf("GLOPS/s = %f \n", gflops);
	double peak = factor*3.1*nTs;
	printf("PEAK    = %f \n", peak);
    printf("EFF o/o = %f \n", gflops/peak*100.);
}

real Solver::get_sweeptime(){
	return time_used_sweep;
}

real Solver::get_totaltime(){
	return time_used_total;
}
