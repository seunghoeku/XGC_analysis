//============================================================================
//  Copyright (c) Kitware, Inc.
//  All rights reserved.
//  See LICENSE.txt for details.
//
//  This software is distributed WITHOUT ANY WARRANTY; without even
//  the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
//  PURPOSE.  See the above copyright notice for more information.
//============================================================================

#ifndef XGCParameters_h
#define XGCParameters_h

struct XGCParameters
{
  int numPlanes = -1;
  int numNodes = -1;
  int numTri = -1;

  double eq_axis_r = -1, eq_axis_z = -1, eq_x_psi = -1, eq_x_r = -1, eq_x_z = -1;
  double eq_min_r = -1, eq_max_r = -1;
  double eq_min_z = -1, eq_max_z = 1;
  double psi_min = -1, psi_max = -1;
  int eq_mr = -1, eq_mz = -1;
};

#endif
