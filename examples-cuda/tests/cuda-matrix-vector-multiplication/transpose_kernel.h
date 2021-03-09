/*
 * Copyright 1993-2009 NVIDIA Corporation.  All rights reserved.
 *
 * NVIDIA Corporation and its licensors retain all intellectual property and 
 * proprietary rights in and to this software and related documentation. 
 * Any use, reproduction, disclosure, or distribution of this software 
 * and related documentation without an express license agreement from
 * NVIDIA Corporation is strictly prohibited.
 *
 * Please refer to the applicable NVIDIA end user license agreement (EULA) 
 * associated with this source code for terms and conditions that govern 
 * your use of this NVIDIA software.
 * 
 */

/* Matrix transpose with Cuda
* Device code.
*/

#ifndef _TRANSPOSE_KERNEL_H_
#define _TRANSPOSE_KERNEL_H_

#define BLOCK_DIM 16

// This kernel is optimized to ensure all global reads and writes are coalesced,
// and to avoid bank conflicts in shared memory.  This kernel is up to 11x faster
// than the naive kernel below.  Note that the shared memory array is sized to 
// (BLOCK_DIM+1)*BLOCK_DIM.  This pads each row of the 2D block in shared memory 
// so that bank conflicts do not occur when threads address the array column-wise.
__global__ void transpose(float *odata, float *idata, int width, int height);


// This naive transpose kernel suffers from completely non-coalesced writes.
// It can be up to 10x slower than the kernel above for large matrices.
__global__ void transpose_naive(float *odata, float* idata, int width, int height);

#endif // _TRANSPOSE_KERNEL_H_
