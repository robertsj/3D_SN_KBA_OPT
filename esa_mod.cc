#include "miniapp.hh"
#include "auxiliary_function.hh"
#include <numeric>
#include <string.h>
//#define NOBOUNDX
//#define NOBOUNDY
#define NOBOUNDZ

#define ALIGN __attribute__((aligned(64)))
//#define ALIGN

void Solver::sweep_esa_mod()
{
  std::fill(phi.begin(), phi.end(), 0.0);
  real phi_t[phi.size()] ALIGN;
  for (int i = 0; i < phi.size(); ++i)
    phi_t[i] = 0.0;

#ifndef NOBOUNDX
  real bd_info_x[mesh.nbz][mesh.nby][mesh.zbs[0]][mesh.ybs[0]][n_a] ALIGN;
#endif
#ifndef NOBOUNDY
  real bd_info_y[mesh.nbz][mesh.nbx][mesh.zbs[0]][mesh.xbs[0]][n_a] ALIGN;
#endif
#ifndef NOBOUNDZ
  real bd_info_z[mesh.nby][mesh.nbx][mesh.ybs[0]][mesh.xbs[0]][n_a] ALIGN;
#endif


  real muDelta[n_a] ALIGN;
  real etaDelta[n_a] ALIGN;
  real xiDelta[n_a]  ALIGN;
  real sum[n_a] ALIGN;
  real wtt[n_a] ALIGN;
  for(int a = 0; a < n_a; a++)
  {
    muDelta[a]  = 2.0 * mu[a];
    etaDelta[a] = 2.0 * eta[a];
    xiDelta[a]  = 2.0 * xi[a];
    sum[a]      = muDelta[a] + etaDelta[a] + xiDelta[a];
    wtt[a]      = wt[a];
  }



  // material index
  const int *mat_id = &mesh.mat_id[0];

#pragma omp parallel num_threads(nTs)
  {
    int TID = omp_get_thread_num();

#ifdef NOBOUNDZ
    real bd_info_z[mesh.ybs[0]][mesh.xbs[0]][n_a];
#endif

    // GROUPS
    for (int g = 0; g < n_eg; ++g)
    {
      // OCTANTS
      for (int o = 0; o < 8; o++)
      {
#pragma omp master
        {
          // set all cell edge fluxes to zero.  this is consistent with vacuum.

#ifndef NOBOUNDX
          SetValue(&bd_info_x[0][0][0][0][0],
              mesh.nbz*mesh.nby*mesh.zbs[0]*mesh.ybs[0]*n_a, 0.0);
#endif
#ifndef NOBOUNDY
          SetValue(&bd_info_y[0][0][0][0][0],
              mesh.nbz*mesh.nbx*mesh.zbs[0]*mesh.xbs[0]*n_a, 0.0);
#endif
#ifndef NOBOUNDZ
          SetValue(&bd_info_z[0][0][0][0][0],
              mesh.nby*mesh.nbx*mesh.ybs[0]*mesh.xbs[0]*n_a, 0.0);
#else
          SetValue(&bd_info_z[0][0][0],
              mesh.ybs[0]*mesh.xbs[0]*n_a, 0.0);
#endif
        }

#ifdef NOBOUNDZ
        SetValue(&bd_info_z[0][0][0], mesh.ybs[0]*mesh.xbs[0]*n_a, 0.0);
#endif
        // PLANES
        for (int p = 0; p < mesh.planes[o].size(); ++p)
        {
          // BLOCKS
#pragma omp for
          for (int b = 0; b < mesh.planes[o][p].size(); ++b)
          {
            // block id's
            int i_b = mesh.planes[o][p][b][mesh.X];
            int j_b = mesh.planes[o][p][b][mesh.Y];
            int k_b = mesh.planes[o][p][b][mesh.Z];
            // Z
            for (int kk = 0; kk < mesh.zbs[k_b]; ++kk)
            {
              int k = mesh.get_block_k(o, k_b, kk);
              real one_over_dz = 1/mesh.dz[k];

              // Y
              for (int jj = 0; jj < mesh.ybs[j_b]; ++jj)
              {
                int j = mesh.get_block_j(o, j_b, jj);
                real one_over_dy = 1/mesh.dy[j];

                // X
                for (int ii = 0; ii < mesh.xbs[i_b]; ++ii)
                {
                  int i = mesh.get_block_i(o, i_b, ii);
                  int cell = mesh.ijk_to_cell(i, j, k);
                  int m = mat_id[cell];
                  real sig_t = SigT[m * n_eg + g];
                  real q = Q[cell * n_eg + g];
                  real c[n_a] ALIGN;
                  real phi_tmp[n_a] ALIGN;
#ifdef NOBOUNDX
                  real bdx = mesh.dx[i];
#endif
#ifdef NOBOUNDY
                  real bdy = mesh.dy[j];
#endif
#ifdef NOBOUNDZ
                  real bdz = mesh.dz[k];
#endif
                  real one_over_dx = 1/mesh.dx[i];

                  // ANGLES
#pragma vector aligned
#pragma omp simd
                  for (int a = 0; a < n_a; a++)
                  {

                    //========================================================//
                    // THE WORK!
                    real eps_x = one_over_dx * muDelta[a];     // 1
                    real eps_y = one_over_dy * etaDelta[a];    // 1
                    real eps_z = one_over_dz * xiDelta[a];     // 1
                    real sum_eps = eps_x + eps_y + eps_z;      // 2
                    real denom = 1.0 / (sig_t + sum_eps);      // 1 + div  ... 6+div
                    c[a]  = q;
                    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~//
#ifdef NOBOUNDX
                    c[a] += eps_x * bdx;
#else
                    c[a] += eps_x * bd_info_x[k_b][j_b][kk][jj][a]; // 2
#endif
#ifdef NOBOUNDY
                    c[a] += eps_y * bdy;
#else
                    c[a] += eps_y * bd_info_y[k_b][i_b][kk][ii][a]; // 2
#endif
#ifdef NOBOUNDZ
                    //c[a] += eps_z * bdz;
                    c[a] += eps_z * bd_info_z[jj][ii][a]; // 2
#else
                    c[a] += eps_z * bd_info_z[j_b][i_b][jj][ii][a]; // 2     ... 12+div
#endif
                    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~//

                    c[a] *= denom; // 1
                    phi_t[cell * n_eg + g] += wtt[a] * c[a]; // 2 ... 18+div

                    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~//
#ifdef NOBOUNDX
                    bdx = 2.0*c[a] - bdx; // 2
#else
                    bd_info_x[k_b][j_b][kk][jj][a] = 2.0*c[a] - bd_info_x[k_b][j_b][kk][jj][a]; // 2
#endif
#ifdef NOBOUNDY
                    bdy = 2.0*c[a] - bdy; // 1
#else
                    bd_info_y[k_b][i_b][kk][ii][a] = 2.0*c[a] - bd_info_y[k_b][i_b][kk][ii][a]; // 1
#endif
#ifdef NOBOUNDZ
                    //bdz = 2.0*c[a] - bdz; // 1
                    bd_info_z[jj][ii][a] = 2.0*c[a] - bd_info_z[jj][ii][a];
#else
                    bd_info_z[j_b][i_b][jj][ii][a] = 2.0*c[a] - bd_info_z[j_b][i_b][jj][ii][a]; // 1
#endif
                    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~//


                    //========================================================//
                  } // ANGLES
                } // I
              } // J
            } // K
          } // BLOCKS
        } // PLANES
      } // OCTANTS
    } // GROUPS
  } // PARALLEL
  for (int i = 0; i < phi.size(); ++i)
    phi[i] = phi_t[i];
}
