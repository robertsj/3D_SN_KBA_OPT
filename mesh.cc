#include "mesh.hh"
#include <cstdio>
#include <cstdlib>
#include <time.h>
#include <algorithm>
#include "auxiliary_function.hh"

//----------------------------------------------------------------------------//
Mesh::Mesh(int cm_xy, int fm_xy, int cm_z, int fm_z, int bs_x, int bs_y, int bs_z, int nt_in)
: nx(fm_xy * cm_xy)
, ny(fm_xy * cm_xy)
, nz(fm_z  * cm_z)
, ncell(nx*ny*nz)
, dx(nx, 1.0/nx)
, dy(ny, 1.0/ny)
, dz(nz, 1.0/nz)
, nmat(cm_xy*cm_xy*cm_z)
, mat_id(ncell, 0)
, nbx(nx / bs_x + (nx % bs_x ? 1 : 0))
, nby(ny / bs_y + (ny % bs_y ? 1 : 0))
, nbz(nz / bs_z + (nz % bs_z ? 1 : 0))
, xbs(nbx, bs_x)
, ybs(nby, bs_y)
, zbs(nbz, bs_z)
, i_s(nbx)
, j_s(nby)
, k_s(nbz)
, planes(8)
, chunks(8)
, nt(nt_in)
{

  printf("mesh details: \n");
  printf("    nx = %i\n", nx);
  printf("    ny = %i\n", ny);
  printf("    nz = %i\n", nz);
  printf("  nmat = %i\n", nmat);


  // fix size of last blocks block sizes
  xbs[nbx-1] = nx - bs_x * (nbx - 1);
  ybs[nby-1] = ny - bs_y * (nby - 1);
  zbs[nbz-1] = nz - bs_z * (nbz - 1);

  // randomize coarse mesh materials
  std::srand(time(0));
  vec_int reg_mat(nmat);
  for(int i = 0; i < nmat; ++i)
    reg_mat[i] = i;
  std::random_shuffle(mat_id.begin(), mat_id.end());
  for (int z = 0; z < cm_z; ++z)
  {
    for(int y = 0; y < cm_xy; ++y)
    {
      for (int x = 0; x < cm_xy; ++x)
      {
        int m = reg_mat[z * cm_xy * cm_xy + y * cm_xy + x];
        // assign this material to all fine cells in this coarse cell
        for (int k = z * fm_z; k < (z + 1) * fm_z; ++k)
          for(int j = y * fm_xy; j < (y + 1) * fm_xy; ++j)
            for(int i = x * fm_xy; i < (x + 1) * fm_xy; ++i)
              mat_id[ijk_to_cell(i, j, k)] = m;
      }
    }
  }

  // define planes
  vec2_int column_planes(8, vec_int(nbx*nby));
  int np = nbx + nby + nbz - 2;
  for (int o = 0; o < 8; ++o)
  {
    planes[o].resize(np);
    for (int p = 0; p < np; ++p)
    {
      for (int ii = 0; ii < nbx; ++ii)
      {
        // counter-clock octant
        int i = (o % 4 == 0 or o % 4 == 3) ? ii : nbx - ii - 1;
        for (int jj = 0; jj < nby; ++jj)
        {
          int j = (o % 4 < 2) ? jj : nby - jj - 1;
          for (int kk = 0; kk < nbz; ++kk)
          {
            int k =  (o < 4) ? kk : nbz - kk - 1;
            if (ii + jj + kk == p)
            {
              vec_int ijk(3);
              ijk[0] = i;
              ijk[1] = j;
              ijk[2] = k;
              planes[o][p].push_back(ijk);
              if (kk == 0)
                column_planes[o][ii+jj*nbx] = p;
            }
          }
        }
      }
    }
  }

  // define chunks
  vec_int column_start(nbx*nby, 0);
  vec_int nct(nt, 0);

  vec2_int tmp_chunks;
  for (int c = 0; c < nbx*nby; ++c)
  {
    int i = c % nbx;
    int j = c / nby;
    int t = c % nt;
    int p = column_planes[0][i+j*nbx];
    int step = 0;
    if (nct[t] == 0)
      step = p;
    column_start[c] = nct[t] + step;
    if (p > column_start[c])
      column_start[c] = p;
    //printf("I=%i, J=%i, C=%i, T=%i, P=%i, CS=%i\n",i,j,c,t,p,column_start[c]);
    for (int k = 0; k < nbz; ++k)
    {
      int ci = nct[t]+step;
      vec_int chunk(4, 0);
      chunk[0] = i;
      chunk[1] = j;
      chunk[2] = k;
      chunk[3] = ci;
      tmp_chunks.push_back(chunk);
      nct[t]++;
    }
    nct[t] += step;
  } //

//  for (int i = 0; i < tmp_chunks.size(); ++i)
//  {
//    printf("(%i, %i, %i, %i) \n",
//        tmp_chunks[i][0],tmp_chunks[i][1],tmp_chunks[i][2],tmp_chunks[i][3]);
//  }

  int num_chunks = tmp_chunks[tmp_chunks.size()-1][3]+1;
  for (int o = 0; o < 8; ++o)
  {
    chunks[o].resize(num_chunks);
    for (int c = 0; c < tmp_chunks.size(); ++c)
    {
      int I = tmp_chunks[c][0];
      int J = tmp_chunks[c][1];
      int K = tmp_chunks[c][2];
      int C = tmp_chunks[c][3];
      if (o == 1 or o == 2 or o == 5 or o == 6)
        I = nbx - I - 1;
      if (o == 2 or o == 3 or o == 6 or o == 7)
        J = nby - J - 1;
      if (o > 3)
        K = nbz - K - 1;
      vec_int chunk(3);
      chunk[0] = I; chunk[1] = J; chunk[2] = K;
      //printf(" c = %i, C = %i\n", c, C);
      chunks[o][C].push_back(chunk);
    }
  }

#ifdef CHUNKS
  planes = chunks;
#endif

//  for (int o = 0; o < 8; ++o)
//  {
//    printf("octant %i\n", o);
//  for (int i = 0; i < chunks[0].size(); ++i)
//  {
//    printf("  chunk %i \n", i);
//    for (int j = 0; j < chunks[0][i].size(); ++j)
//    {
//      printf("   (%i, %i, %i) \n",
//        chunks[o][i][j][0], chunks[o][i][j][1], chunks[o][i][j][2]);
//    }
//  }
//  }

  // block indices
  i_s.resize(nbx, 0);
  j_s.resize(nby, 0);
  k_s.resize(nbz, 0);

  for (int i = 0; i < nbx; ++i)
    for (int ii = 0; ii < i; ++ii)
      i_s[i] += xbs[ii];

  for (int j = 0; j < nby; ++j)
    for (int jj = 0; jj < j; ++jj)
      j_s[j] += ybs[jj];

  for (int k = 0; k < nbz; ++k)
    for (int kk = 0; kk < k; ++kk)
      k_s[k] += zbs[kk];
}



//----------------------------------------------------------------------------//
void Mesh::print_blocks() const
{
  print_dline();
  for (int o = 0; o < 8; ++o)
  {
    print_dline();
    for (int p = 0; p < planes[o].size(); ++p)
    {
      vec2_int plane = planes[o][p];
      for (int b = 0; b < plane.size(); ++b)
      {
        vec_int block = plane[b];
        int i0=get_block_i(o, block[X], 0);
        int i1=get_block_i(o, block[X], xbs[block[X]]-1);
        int j0=get_block_j(o, block[Y], 0);
        int j1=get_block_j(o, block[Y], ybs[block[Y]]-1);
        int k0=get_block_k(o, block[Z], 0);
        int k1=get_block_k(o, block[Z], zbs[block[Z]]-1);
        printf("octant=%i, plane=%i, block=%i, i=%i-->%i, j=%i-->%i, k=%i-->%i\n",
            o, p, b, i0, i1, j0, j1, k0, k1);
      }
      print_line();
    }
    print_dline();
  }

  for (int i = 0; i < k_s.size(); ++i)
  {
    printf("%i %i\n", i, zbs[i]);
  }
}

//----------------------------------------------------------------------------//
void Mesh::print_summary() const
{
  print_dline();
  for (int o = 0; o < 1; ++o)
  {
    print_dline();
    printf("DECOMPOSITION HAS %i PLANES \n", planes[o].size());
    for (int p = 0; p < planes[o].size(); ++p)
    {
      printf("   PLANE %i HAS %i BLOCKS \n", p, planes[o][p].size());
    }
    print_dline();
  }
}

