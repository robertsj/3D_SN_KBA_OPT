/*
 * eas.cc
 *
 *  Created on: Dec 9, 2015
 *      Author: kevin
 */
#include "miniapp.hh"
#include "auxiliary_function.hh"

void Solver::sweep_eas(int start_TID[])
{
	std::fill(phi.begin(), phi.end(), 0.0);
	bool forward_x, forward_y, forward_z;

	//The arrangement of bd_info_x or _y is [blockID_X][blockID_Y][blockID_Z][z][x or y] to realize unit strid
	//one row is wasted along x & y direction
	real bd_info_x[(N + 1) * (N + 1) * n_b * blocksize_z * blocksize_xy], bd_info_y[(N + 1) * (N + 1) * (n_b) * blocksize_z * blocksize_xy];
	real ch, cv[mesh.ny], cz[mesh.nx * mesh.ny];//used in the remain sweep after block sweep

	real muDelta[n_a], etaDelta[n_a], xiDelta[n_a], sum[n_a];
	for(int a = 0; a < n_a; a++)
	{
		muDelta[a] = 2.0 * mu[a] / mesh.dy;
		etaDelta[a] = 2.0 * eta[a] / mesh.dx;
		xiDelta[a] = 2.0 * xi[a] / mesh.dz;
		sum[a] = muDelta[a] + etaDelta[a] + xiDelta[a];
	}

	//Energy loop
	for(int e = 0; e < n_eg; ++e)
	{
		//Octant loop
		for(int o = 0; o < 8; o++)
		{
			DetermineDir(o, forward_z, forward_x, forward_y);
			//Angle loop
			for(int a = 0; a < n_a; a++)
			{
				//for simplicity, set all to 0, in fact only the start bd_info need to set 0
				SetValue(bd_info_x, (N + 1) * (N + 1) * n_b * blocksize_z * blocksize_xy, 0.0);
				SetValue(bd_info_y, (N + 1) * (N + 1) * n_b * blocksize_z * blocksize_xy, 0.0);

#pragma omp parallel num_threads(nTs)
				{
					//thread-private variables
					int TID = omp_get_thread_num();
					int start_plane, end_plane;
					Find_start_plane(start_TID, N, TID, start_plane);
					end_plane = start_plane + n_b; //working plane for TID is [start_plane, end_plane);
					//block ID
					int blockID_x, blockID_y, blockID_z;
					Set_block_ID_xy(start_TID, N, start_plane, TID, blockID_x, blockID_y);
					real bd_info_z[blocksize_xy * blocksize_xy]; //[X][Y]
					SetValue(bd_info_z, blocksize_xy * blocksize_xy, 0.0);
					real c;

					for(int p = 0; p < n_p; p++)
					{
						if (p >= start_plane && p < end_plane)//select working threads
						{
							blockID_z = p - start_plane;
							//block sweep
							for(int z0 = blockID_z * blocksize_z; z0 < blockID_z * blocksize_z + blocksize_z; z0++) //always use global ordinates
							{
								const int z_global = forward_z ? z0 : mesh.nz - 1 - z0;
								const int z_local = z_global % blocksize_z;
								for(int x0 = blockID_x * blocksize_xy; x0 < blockID_x * blocksize_xy + blocksize_xy; x0++)
								{
									const int x_global = forward_x ? x0 : mesh.nx - 1 - x0;
									const int x_local = x_global % blocksize_xy;
									const int index_x = Get_index(N, n_b, blocksize_z, blocksize_xy, 1, blockID_x, blockID_y, blockID_z, z_local, x_local);
									for(int y0 = blockID_y * blocksize_xy; y0 < blockID_y * blocksize_xy + blocksize_xy; y0++)
									{
										const int y_global = forward_y ? y0 : mesh.ny - 1 - y0;
										const int y_local = y_global % blocksize_xy;
										const int m = mesh.mat_id[z_global * mesh.nx * mesh.ny + x_global * mesh.ny + y_global];
										const int index_y = Get_index(N, n_b, blocksize_z, blocksize_xy, 1, blockID_x, blockID_y, blockID_z, z_local, y_local);
										const int index_z = x_local * blocksize_xy + y_local;
										c = (muDelta[a] * bd_info_x[index_x] +	etaDelta[a] * bd_info_y[index_y] +	xiDelta[a] * bd_info_z[index_z] +
												Q[z_global * mesh.nx * mesh.ny * n_eg + x_global * mesh.ny * n_eg + y_global * n_eg + e]) / (SigT[m * n_eg + e] + sum[a]);
										phi[z_global * mesh.nx * mesh.ny * n_eg + x_global * mesh.ny * n_eg + y_global * n_eg + e] += wt[a] * c;
										bd_info_x[index_x] = 2.0 * c - bd_info_x[index_x];
										bd_info_y[index_y] = 2.0 * c - bd_info_y[index_y];
										bd_info_z[index_z] = 2.0 * c - bd_info_z[index_z];
									}
								}
							}
							//after block sweep, copy the bd_info_x and _y from the memory of TID into the corresponding memory of (TID + 1)
							copy(bd_info_x, Get_index(N, n_b, blocksize_z, blocksize_xy, 1, blockID_x, blockID_y + 1, blockID_z, 0, 0),
									bd_info_x, Get_index(N, n_b, blocksize_z, blocksize_xy, 1, blockID_x, blockID_y, blockID_z, 0, 0), blocksize_z * blocksize_xy);
							copy(bd_info_y, Get_index(N, n_b, blocksize_z, blocksize_xy, 1, blockID_x + 1, blockID_y, blockID_z, 0, 0),
									bd_info_y, Get_index(N, n_b, blocksize_z, blocksize_xy, 1, blockID_x, blockID_y, blockID_z, 0, 0), blocksize_z * blocksize_xy);
						}
#pragma omp barrier
					}
					if(remain != 0)
					{
						if (forward_x == true && forward_y == true)
							for(int x0 = 0; x0 < blocksize_xy; x0++)
								for(int y0 = 0; y0 < blocksize_xy; y0++)
									cz[(blockID_x * blocksize_xy + x0) * mesh.ny + blockID_y * blocksize_xy + y0] =
											bd_info_z[x0 * blocksize_xy + y0];
						else if(forward_x == true && forward_y == false)
							for(int x0 = 0; x0 < blocksize_xy; x0++)
								for(int y0 = 0; y0 < blocksize_xy; y0++)
									cz[(blockID_x * blocksize_xy + x0) * mesh.ny + (N - 1 - blockID_y) * blocksize_xy + y0] =
											bd_info_z[x0 * blocksize_xy + y0];
						else if(forward_x == false && forward_y == true)
							for(int x = 0; x < blocksize_xy; x++)
								for(int y = 0; y < blocksize_xy; y++)
									cz[((N - 1 - blockID_x) * blocksize_xy + x) * mesh.ny + blockID_y * blocksize_xy + y] =
											bd_info_z[x * blocksize_xy + y];
						else
							for(int x = 0; x < blocksize_xy; x++)
								for(int y = 0; y < blocksize_xy; y++)
									cz[((N - 1 - blockID_x) * blocksize_xy + x) * mesh.ny + (N - 1 - blockID_y) * blocksize_xy + y] =
											bd_info_z[x * blocksize_xy + y];
					}
				}
				//remain sweep
				if (remain != 0)
				{
					real c;
					for(int z0 = 0; z0 < remain; z0++)
					{
						const int z = forward_z ? n_b * blocksize_z + z0 : remain - 1 - z0;
						SetValue(cv, mesh.ny, 0.0);
						for(int x0 = 0; x0 < mesh.nx; x0++)
						{
							const int x = forward_x ? x0 : mesh.nx - 1 - x0;
							ch = 0.0;
							//Space loop:y
							for(int y0 = 0; y0 < mesh.ny; y0++)
							{
								const int y = forward_y ? y0 : mesh.ny - 1 - y0;
								const int m = mesh.mat_id[z * mesh.nx * mesh.ny + x * mesh.ny + y];
								c = (muDelta[a] * ch + etaDelta[a] * cv[y] + xiDelta[a] * cz[x * mesh.ny + y] +
										Q[z * mesh.nx * mesh.ny * n_eg + x * mesh.ny * n_eg + y * n_eg + e])
																																																								/ (SigT[m * n_eg + e] + sum[a]);
								phi[z * mesh.nx * mesh.ny * n_eg + x * mesh.ny * n_eg + y * n_eg + e] += 1.5707963267948966 / n_a * c;
								ch = 2.0 * c - ch;
								cv[y] = 2.0 * c - cv[y];
								cz[x * mesh.ny + y] = 2.0 * c - cz[x * mesh.ny + y];
							}
						}
					}
				}
			}
		}
	}
}



