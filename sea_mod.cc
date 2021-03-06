#include "miniapp.hh"
#include "auxiliary_function.hh"
#include <numeric>

void Solver::sweep_sea_mod()
{
	std::fill(phi.begin(), phi.end(), 0.0);
	// use, e.g., ybs[0] because only ybs[-1] could be smaller
	// x-directed  Z         Y         z size       y size
	real bd_info_x[mesh.nbz][mesh.nby][mesh.zbs[0]][mesh.ybs[0]][n_eg][n_a];
	// y-directed  Z         X         z size       x size
	real bd_info_y[mesh.nbz][mesh.nbx][mesh.zbs[0]][mesh.xbs[0]][n_eg][n_a];
	// z-directed  Y         X         y size       x size
	real bd_info_z[mesh.nby][mesh.nbx][mesh.ybs[0]][mesh.xbs[0]][n_eg][n_a];

	real muDelta[n_a], etaDelta[n_a], xiDelta[n_a], sum[n_a];
	for(int a = 0; a < n_a; a++)
	{
		muDelta[a]  = 2.0 * mu[a] / mesh.dx;
		etaDelta[a] = 2.0 * eta[a] / mesh.dy;
		xiDelta[a]  = 2.0 * xi[a] / mesh.dz;
		sum[a]      = muDelta[a] + etaDelta[a] + xiDelta[a];
	}

	// material index
	const int *mat_id = &mesh.mat_id[0];

	// OCTANTS
	for (int o = 0; o < 8; o++)
	{
		// set all cell edge fluxes to zero.  this is consistent with vacuum.
		SetValue(&bd_info_x[0][0][0][0][0][0],
				mesh.nbz*mesh.nby*mesh.zbs[0]*mesh.ybs[0]*n_eg*n_a, 0.0);
		SetValue(&bd_info_y[0][0][0][0][0][0],
				mesh.nbz*mesh.nbx*mesh.zbs[0]*mesh.xbs[0]*n_eg*n_a, 0.0);
		SetValue(&bd_info_z[0][0][0][0][0][0],
				mesh.nby*mesh.nbx*mesh.ybs[0]*mesh.xbs[0]*n_eg*n_a, 0.0);
#pragma omp parallel num_threads(nTs)
		{
			int TID = omp_get_thread_num();
			// PLANES
			for (int p = 0; p < mesh.planes[o].size(); ++p)
			{
				// BLOCKS
				for (int b = 0; b < mesh.planes[o][p].size(); ++b)
				{
					// SKIP THIS BLOCK IF NOT MINE
					if (b % nTs != TID)
						continue;
					// block id's
					int i_b = mesh.planes[o][p][b][mesh.X];
					int j_b = mesh.planes[o][p][b][mesh.Y];
					int k_b = mesh.planes[o][p][b][mesh.Z];
					// Z
					for (int kk = 0; kk < mesh.zbs[k_b]; ++kk)
					{
						int k = mesh.get_block_k(o, k_b, kk);
						// Y
						for (int jj = 0; jj < mesh.ybs[j_b]; ++jj)
						{
							int j = mesh.get_block_j(o, j_b, jj);
							// X
							for (int ii = 0; ii < mesh.xbs[i_b]; ++ii)
							{
								int i = mesh.get_block_i(o, i_b, ii);
								int cell = mesh.ijk_to_cell(i, j, k);
								int m = mat_id[cell];
								// GROUPS
								for (int g = 0; g < n_eg; ++g)
								{
									real c[n_a];
									v_dbl phi_tmp(n_a);
#pragma omp simd
									// ANGLES
									for(int a = 0; a < n_a; a++)
									{
										//========================================================//
										// THE WORK!
										c[a] = ( muDelta[a] * bd_info_x[k_b][j_b][kk][jj][g][a] +
												    etaDelta[a] * bd_info_y[k_b][i_b][kk][ii][g][a] +
												     xiDelta[a] * bd_info_z[j_b][i_b][jj][ii][g][a] +
												     Q[cell * n_eg + g]) / (SigT[m * n_eg + g] + sum[a]);
										phi_tmp[a] = wt[a] * c[a];
										bd_info_x[k_b][j_b][kk][jj][g][a] = 2.0*c[a] - bd_info_x[k_b][j_b][kk][jj][g][a];
										bd_info_y[k_b][i_b][kk][ii][g][a] = 2.0*c[a] - bd_info_y[k_b][i_b][kk][ii][g][a];
										bd_info_z[j_b][i_b][jj][ii][g][a] = 2.0*c[a] - bd_info_z[j_b][i_b][jj][ii][g][a];
										//========================================================//
									} // ANGLES
									phi[cell * n_eg + g] += std::accumulate(phi_tmp.begin(), phi_tmp.end(), 0.0);
								} // GROUPS
							} // I
						} // J
					} // K
				} // BLOCKS
#pragma omp barrier
			} // PLANES
		} // PARALLEL
	} // OCTANTS
}
