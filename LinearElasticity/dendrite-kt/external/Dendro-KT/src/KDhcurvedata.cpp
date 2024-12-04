//
// @author Milinda Fernando
// @author Masado Ishii
// School of Computing, University of Utah

// NOTE: DO NOT CHANGE THIS FILE FOR ANY REASON.


// This header file contains all the rotation permutations + hilbert rotation table data hard corded to improve the performance of the hilbert curve.
// Note that: Rotations contains the concatenated strings of rot_perm and rot_index.
// Created by Milinda Fernando
// on 10/2/15.
//
// Regenerated by Masado Ishii
// on 01/18/2019
//

#include "../include/hcurvedata.h"
#include "../include/KDhcurvedata_decl.h"
#include <assert.h>

char* rotations;
int* HILBERT_TABLE;

void _InitializeHcurve(int pDim)
{
#ifdef HILBERT_ORDERING

  switch(pDim)
  {
    case 2:
      HilbertData<2>::copyData(rotations, HILBERT_TABLE);
      break;

    case 3:
      HilbertData<3>::copyData(rotations, HILBERT_TABLE);
      break;

    case 4:
      HilbertData<4>::copyData(rotations, HILBERT_TABLE);
      break;

    default:
      assert(false);
  }

#else

  const int num_orthant = (1u<<pDim);
  rotations = new char[2*num_orthant];
  HILBERT_TABLE = new int[num_orthant];
  for (int ort = 0; ort < num_orthant; ort++)
  {
    rotations[ort] = ort;
    rotations[ort + num_orthant] = ort;
    HILBERT_TABLE[ort] = 0;
  }

#endif
}

void _DestroyHcurve()
{
  if (rotations != NULL)
  {
    delete [] rotations;
    rotations = NULL;
  }
  if (HILBERT_TABLE != NULL)
  {
    delete [] HILBERT_TABLE;
    HILBERT_TABLE = NULL;
  }
}