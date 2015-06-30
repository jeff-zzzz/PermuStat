# -*- coding: utf-8 -*-
"""
Created on Fri Mar 13 10:09:07 2015
@author: jesong1126
"""


ExtractLFmatrix_angle (unsigned long int nE, double *originXYZ, double *ElecXYZ, unsigned long int nE_big, double *ElecXYZ_big, double *originXYZ_big, unsigned long int nS, double *LF, double *LF_big, 
						double cosAngle, double smallCosAngle, unsigned long int maxN, double exponent) {

// PLEASE NOTE: ElecXYZ will be altered (by centering on originXYZ)
// PLEASE NOTE: ElecXYZ_big will be altered (by centering on originXYZ_big)
// originXYZ = address of first element of "center" of subject's electrode cloud = [x y z] -> size: 1 x 3 {vector} 
// ElecXYZ = address of first element of electrode list for subject = [x1, y1, z1, x2, y2, z2, ...] -> size: nE x 3 {row major}
// originXYZ_big = address of first element of "center" of electrode list used to calculate LF_big = [x y z] -> size: 1 x 3 {vector}
// ElecXYZ_big = address of first element of electrode list used to calculate LF_big = [x1, y1, z1, x2, y2, z2, ...] -> size: nE_big x 3 {row major}
// LF = address of first element of forward model that will be filled by function -> size: nE x nS {row major}
// LF_big = address of first element of big forward model -> size: nE_big x nS {row major}

	unsigned long int i, i_big, size, k, pos, posE_max, pos2;
	unsigned long int *NPosE;
	unsigned long int *PosE;
	double dotP, dotP_max, norm, oneMinusCosAngle, oneMinusCosAngle_small, sum;
	double *Norm_big, *DotP;
	
	size = nE * maxN * sizeof(unsigned long int);
	PosE = (unsigned long int *) malloc(size);
	memset(PosE,0,size);
	
	size = nE * sizeof(unsigned long int);
	NPosE = (unsigned long int *) malloc(size);
	memset(NPosE,0,size);

	size = nE * maxN * sizeof(double);
	DotP = (double *) malloc(size);
	memset(DotP,0,size);
	
	Norm_big = (double *) malloc(nE_big * sizeof(double));
	if (originXYZ_big != NULL) {
		for (i_big = 0; i_big < nE_big; i_big++) {
// void cblas_daxpy(const int N, const double alpha, const double *X, const int incX, double *Y, const int incY);
			cblas_daxpy(3,-1.0,originXYZ_big,1,&ElecXYZ_big[i_big * 3],1);
// double cblas_dnrm2(const int N, const double *X, const int incX);
			Norm_big[i_big] = cblas_dnrm2(3,&ElecXYZ_big[i_big * 3],1);
		}
	}
	else {
		for (i_big = 0; i_big < nE_big; i_big++) {
// double cblas_dnrm2(const int N, const double *X, const int incX);
			Norm_big[i_big] = cblas_dnrm2(3,&ElecXYZ_big[i_big * 3],1);
		}
	}
	
	oneMinusCosAngle = 1.0 - fabs(cosAngle);
	oneMinusCosAngle_small = 1.0 - fabs(smallCosAngle);
	
	for (i = 0; i < nE; i++) {
		pos = i * maxN;
		k = 0;
		dotP_max = 0.0;
// void cblas_daxpy(const int N, const double alpha, const double *X, const int incX, double *Y, const int incY);
		if (originXYZ != NULL) cblas_daxpy(3,-1.0,originXYZ,1,&ElecXYZ[i * 3],1);
// double cblas_dnrm2(const int N, const double *X, const int incX);
		norm = cblas_dnrm2(3,&ElecXYZ[i * 3],1);
		for (i_big = 0; i_big < nE_big; i_big++) {
// double cblas_ddot(const int N, const double *X, const int incX, const double *Y, const int incY);
// double cblas_dnrm2(const int N, const double *X, const int incX);
			dotP = 1.0 - cblas_ddot(3,&ElecXYZ[i * 3],1,&ElecXYZ_big[i_big * 3],1) / Norm_big[i_big] / norm;
			if (dotP <= oneMinusCosAngle_small) {
				PosE[pos] = i_big;
				k = 1;
				break;
			}
			else if (dotP <= oneMinusCosAngle) {
				if (k < maxN) {
					PosE[pos + k] = i_big;
					DotP[pos + k] = dotP;
					if (dotP > dotP_max) {
						posE_max = pos + k;
						dotP_max = dotP;
					}
					k += 1;
				}
				else if (dotP < dotP_max) {
					PosE[posE_max] = i_big;
					DotP[posE_max] = dotP;
// CBLAS_INDEX cblas_idamax(const int N, const double *X, const int incX);
					posE_max = pos + (unsigned long int)cblas_idamax((int)maxN,&DotP[pos],1);
					dotP_max = DotP[posE_max];
				}
			}
		}
		NPosE[i] = k;
	}
	
	for (i = 0; i < nE; i++) {
		pos = i * nS;
		memset(&LF[pos],0,nS * sizeof(double));
		if (NPosE[i] == 0) break;
		pos2 = i * maxN;
// void cblas_daxpy(const int N, const double alpha, const double *X, const int incX, double *Y, const int incY);
		if (NPosE[i] == 1) cblas_daxpy((int)nS,1.0,&LF_big[PosE[pos2] * nS],1,&LF[pos],1);
		else {
			sum = DotP[pos2] = 1.0 / pow(DotP[pos2],exponent);
			for (k = 1; k < NPosE[i]; k++) {
				DotP[pos2 + k] = 1.0 / pow(DotP[pos2 + k],exponent);
				sum += DotP[pos2 + k];
			}
			for (k = 0; k < NPosE[i]; k++) {
// void cblas_daxpy(const int N, const double alpha, const double *X, const int incX, double *Y, const int incY);
				cblas_daxpy((int)nS,DotP[pos2 + k]/sum,&LF_big[PosE[pos2 + k] * nS],1,&LF[pos],1);
			}
		}
	}
	
	free(PosE);
	free(NPosE);
	free(DotP);
	free(Norm_big);
}


