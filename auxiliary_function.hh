/*
 * auxiliary_function.hh
 *
 *  Created on: Aug 12, 2015
 *      Author: kevin
 */

#ifndef AUXILIARY_FUNCTION_HH_
#define AUXILIARY_FUNCTION_HH_

#ifdef FLOAT
typedef float real;
#else
typedef double real;
#endif

void SetValue(real *array, int length, real value);
void DetermineDir(int o, bool &forward_z, bool &forward_x, bool &forward_y);
void copy(real *left, int lstart, real *right, int rstart, int length);
void Set_start_TID(int start_TID[], int N);
void Find_start_plane(int start_TID[], int N, int TID, int &start_plane);
void Set_block_ID_xy(int start_TID[], int N, int start_plane, int TID,
		int &block_x, int &block_y);
void print_line(const int n = 80);
void print_dline(const int n = 80);
#endif /* AUXILIARY_FUNCTION_HH_ */
