#include <stdio.h>
#include "CKKS_CoeffToSlot.h"
#include "CKKS_poly.h"
#include "CKKS_EvalMod.h"
#include "CKKS_SlotToCoeff.h"
#include <iostream>
#include <typeinfo>
using namespace std;

#define logN 10
#define MDNUM 3
#define N (1<<logN)
#define L 25
#define DNUM 5
#define K (L/DNUM)
#define Delta (uint64_t) (1ULL << 58)


double norm_square_exp(const double zr[N / 2], const double zi[N / 2]) {
	double sum = 0;
	for (int i = 0; i < N / 2; i++) sum += zr[i] * zr[i] + zi[i] * zi[i];
	return (sum / (N / 2));
}

double norm_square(const double zr[N / 2], const double zi[N / 2]) {
	double sum = 0;
	for (int i = 0; i < N / 2; i++)
		sum += zr[i] * zr[i] + zi[i] * zi[i];
	return (sum);
}

double norm_max(const double zr[N / 2], const double zi[N / 2]) {
	double max = 0;
	for (int i = 0; i < N / 2; i++) {
		double abs = sqrt(zr[i] * zr[i] + zi[i] * zi[i]);
		if (abs > max) max = abs;
	}
	return max;
}

void main() {
	//---------------------------------------------------------------------
	// Initial setting
	//---------------------------------------------------------------------
	srand((unsigned int)time(NULL));

	uint64_t q[L];	uint64_t p[K];
	{
		uint64_t q58[30] = { 288230376155250689ULL, 288230376147386369ULL, 288230376138735617ULL, 288230376160755713ULL, 288230376158396417ULL, 288230376161280001ULL, 288230376161673217ULL, 288230376135196673ULL, 288230376162459649ULL, 288230376132182017ULL, 288230376131788801ULL, 288230376164294657ULL, 288230376167047169ULL, 288230376129691649ULL, 288230376173076481ULL, 288230376175828993ULL, 288230376176091137ULL, 288230376128643073ULL, 288230376176484353ULL, 288230376126545921ULL, 288230376122613761ULL, 288230376178057217ULL, 288230376180023297ULL, 288230376121434113ULL, 288230376121040897ULL, 288230376180416513ULL, 288230376181334017ULL, 288230376184086529ULL, 288230376118026241ULL, 288230376115535873ULL };
		for (int i = 0; i < K; i++)
			p[i] = q58[L+K-1-i];
		for (int i = 0; i < L; i++)
			q[i] = q58[L-1-i];

	}
	printf("58 -bit RNS prime set\n");
	printf("Scaling factor: 2^58 (=%lld)\n", Delta);

	int h = 192; int s[N];
	keygen<N>(h, s);

	//---------------------------------------------------------------------
	// swkgen
	//---------------------------------------------------------------------
//	SparseComplexMatrix<N / 2, 3> E[logN - 1];
//	splitU0R<logN>(E);
//	uint64_t rkey_hat[logN - 1][3][DNUM][2][DNUM * K + K][N];
//	for (int n = 0; n < logN - 1; n++)
//		rkey_gen<N, L, DNUM, K, 3>(E[n], q, p, s, rkey_hat[n]);

	uint64_t  evk_hat[DNUM][2][DNUM * K + K][N];
	//uint64_t ckey_hat[DNUM][2][DNUM * K + K][N];
	int ss[N]; conv<N>(s, s, ss);
	swkgen< N, L, DNUM>(ss, s, q, p, evk_hat);
	//	int sconj[N]; conj<N>(s, sconj);
	//	swkgen<N, L, DNUM>(sconj, s, q, p, ckey_hat);
	printf("KeyGen complete\n");

	// See CKKS_basic.h
	//---------------------------------------------------------------------
	// Ecd & Enc
	//---------------------------------------------------------------------
	uint64_t ct1[2][L][N];
	uint64_t ct2[2][L][N];
	uint64_t ct_tmp[2][L][N];
	uint64_t ct_out[2][L - 1][N];
	uint64_t pt[L][N];
	uint64_t pt_out[L][N];
	uint64_t Delta_out = Delta;
	double z1r[N / 2], z1i[N / 2], z2r[N / 2], z2i[N / 2];
	double wr[N / 2], wi[N / 2];
	double errr[N / 2], erri[N / 2];

	for (int i = 0; i < N / 2; i++) {
		z1r[i] = ((double)rand()) / RAND_MAX;
		z1i[i] = 0;
		z2r[i] = ((double)rand()) / RAND_MAX;
		z2i[i] = 0;
	}

	const int ITER = 100;

	/* Enc Dec */
	printf("\nEnc & Dec Error\n");
	long double max_err = 0;
	for (int iter = 0; iter < ITER; iter++) {

		encode<N, logN, L>(z1r, z1i, Delta, q, pt);
		enc<N, L>(pt, s, q, ct1);

		dec<N, L>(ct1, s, q, pt_out);
		decode<N, logN, L>(pt_out, Delta, q, wr, wi);

		for (int i = 0; i < N / 2; i++) {
			errr[i] = wr[i] - z1r[i];
			erri[i] = 0;
		}
		long double err = std::log2l(norm_max(errr, erri));
		if (err < max_err)
			max_err = err;
	}
	printf("maximal error in Linf norm: %lf\n", max_err);



	/* Add */
	printf("\nAdd error\n");
	max_err = 0;
	for (int iter = 0; iter < ITER; iter++) {

		encode<N, logN, L>(z1r, z1i, Delta, q, pt);
		enc<N, L>(pt, s, q, ct1);
		encode<N, logN, L>(z2r, z2i, Delta, q, pt);
		enc<N, L>(pt, s, q, ct2);

		for (int k = 0; k < 2; ++k) {
			for (int i = 0; i < L; i++) {
				for (int j = 0; j < N; j++) {
					ct_tmp[k][i][j] = (ct1[k][i][j] + ct2[k][i][j]) % q[i];
				}
			}
		}

		dec<N, L>(ct_tmp, s, q, pt_out);
		decode<N, logN, L>(pt_out, Delta, q, wr, wi);

		for (int i = 0; i < N / 2; i++) {
			errr[i] = wr[i] - (z1r[i] + z2r[i]);
			erri[i] = 0;
		}
		long double err = std::log2l(norm_max(errr, erri));
		if (err < max_err)
			max_err = err;
	}
	printf("maximal error in Linf norm: %lf\n", max_err);


	/* Rescale w/ a const scaling factor */
	printf("\nRescale (single) error\n");
	printf("* With a constant scaling factor\n");
	max_err = 0;
	for (int iter = 0; iter < ITER; iter++) {

		encode<N, logN, L>(z1r, z1i, Delta, q, pt);
		enc<N, L>(pt, s, q, ct1);

		for (int k = 0; k < 2; ++k) {
			for (int i = 0; i < L-1; i++) {
				for (int j = 0; j < N; j++) {
					ct_out[k][i][j] = ct1[k][i][j];
				}
			}
		}

		Delta_out = Delta;
		dec<N, L - 1>(ct_out, s, q, pt_out);
		decode<N, logN, L - 1>(pt_out, Delta_out, q, wr, wi);

		for (int i = 0; i < N / 2; i++) {
			errr[i] = wr[i] - z1r[i];
			erri[i] = 0;
		}
		long double err = std::log2l(norm_max(errr, erri));
		if (err < max_err)
			max_err = err;
	}
	printf("maximal error in Linf norm: %lf\n", max_err);



	/* Rescale w/ leveled scaling factors */
	printf("\nRescale (single) error\n");
	printf("* With leveled scaling factors\n");
	max_err = 0;
	for (int iter = 0; iter < ITER; iter++) {

		encode<N, logN, L>(z1r, z1i, Delta, q, pt);
		enc<N, L>(pt, s, q, ct1);
		for (int k = 0; k < 2; ++k) {
			for (int i = 0; i < L; i++) {
				for (int j = 0; j < N; j++) {
					ct_tmp[k][i][j] = mul_mod(Delta, ct1[k][i][j], q[i]);
				}
			}
		}
		RS_hat<N, L, L - 1>(q, ct_tmp, ct_out);
		Delta_out = Delta;
		Delta_out = (uint64_t)((long double)Delta_out * ((long double)Delta_out / (long double)q[L - 1]));
		dec<N, L - 1>(ct_out, s, q, pt_out);
		decode<N, logN, L - 1>(pt_out, Delta_out, q, wr, wi);

		for (int i = 0; i < N / 2; i++) {
			errr[i] = wr[i] - z1r[i];
			erri[i] = 0;
		}
		long double err = std::log2l(norm_max(errr, erri));
		if (err < max_err)
			max_err = err;
	}
	printf("maximal error in Linf norm: %lf\n", max_err);



	/* Mult w/ a const scaling factor */
	printf("\nMult once error\n");
	printf("* With a constant scaling factor\n");
	max_err = 0;
	for (int iter = 0; iter < ITER; iter++) {

		encode<N, logN, L>(z1r, z1i, Delta, q, pt);
		enc<N, L>(pt, s, q, ct1);
		encode<N, logN, L>(z2r, z2i, Delta, q, pt);
		enc<N, L>(pt, s, q, ct2);

		mul<N, L, DNUM, K>(q, p, evk_hat, ct1, ct2, ct_tmp);
		RS_hat<N, L, L - 1>(q, ct_tmp, ct_out);


		Delta_out = Delta;
		dec<N, L - 1>(ct_out, s, q, pt_out);
		decode<N, logN, L - 1>(pt_out, Delta_out, q, wr, wi);


		for (int i = 0; i < N / 2; i++) {
			errr[i] = wr[i] - z1r[i] * z2r[i];
			erri[i] = 0;
		}
		long double err = std::log2l(norm_max(errr, erri));
		if (err < max_err)
			max_err = err;
	}
	printf("maximal error in Linf norm: %lf\n", max_err);



	/* Mult w/ leveled scaling factors */
	printf("\nMult once error\n");
	printf("* With leveled scaling factors\n");
	max_err = 0;
	for (int iter = 0; iter < ITER; iter++) {

		encode<N, logN, L>(z1r, z1i, Delta, q, pt);
		enc<N, L>(pt, s, q, ct1);
		encode<N, logN, L>(z2r, z2i, Delta, q, pt);
		enc<N, L>(pt, s, q, ct2);

		mul<N, L, DNUM, K>(q, p, evk_hat, ct1, ct2, ct_tmp);
		RS_hat<N, L, L - 1>(q, ct_tmp, ct_out);


		Delta_out = Delta;
		Delta_out = (uint64_t)((long double)Delta_out * ((long double)Delta_out / (long double)q[L - 1]));
		dec<N, L - 1>(ct_out, s, q, pt_out);
		decode<N, logN, L - 1>(pt_out, Delta_out, q, wr, wi);


		for (int i = 0; i < N / 2; i++) {
			errr[i] = wr[i] - z1r[i] * z2r[i];
			erri[i] = 0;
		}
		long double err = std::log2l(norm_max(errr, erri));
		if (err < max_err)
			max_err = err;
	}
	printf("maximal error in Linf norm: %lf\n", max_err);

}






/* Mult twice *//*
encode<N, logN, L>(z1r, z1i, Delta, q, pt);
enc<N, L>(pt, s, q, ct1);
encode<N, logN, L>(z2r, z2i, Delta, q, pt);
enc<N, L>(pt, s, q, ct2);
printf("Encryption complete\n");

mul<N, L, DNUM, K>(q, p, evk_hat, ct1, ct2, ct_tmp);
RS_hat<N, L, L - 1>(q, ct_tmp, ct_out);

uint64_t ct_out2[2][L-2][N];
mul<N, L-1, DNUM, K>(q, p, evk_hat, ct_out, ct_out, ct_out);//TODO
RS_hat<N, L-1, L-2>(q, ct_out, ct_out2);

dec<N, L - 2>(ct_out2, s, q, pt_out);
uint64_t Delta_out = Delta;
Delta_out = (uint64_t)((long double)Delta_out * ((long double)Delta_out / (long double)q[L - 1]));
Delta_out = (uint64_t)((long double)Delta_out * ((long double)Delta_out / (long double)q[L - 2]));
decode<N, logN, L - 2>(pt_out, Delta_out, q, wr, wi);
printf("Decryption complete\n");

for (int i = 0; i < N / 2; i++) {
	errr[i] = wr[i] - (z1r[i] * z2r[i]) * (z1r[i] * z2r[i]);
	erri[i] = 0;
	//errr[i] = (wr[i] - ((z1r[i] * z2r[i]) - (z1i[i] * z2i[i]))) / wr[i];
	//erri[i] = (wi[i] - ((z1r[i] * z2i[i]) + (z1i[i] * z2r[i]))) / wi[i];
}
printf("error in Linf norm: %lf\n", -std::log2l(norm_max(errr, erri)));
printf("error in L2 norm: %lf\n", -std::log2l(sqrt(norm_square(errr, erri))));
printf("error in L2_exp norm: %lf\n", -std::log2l(norm_square_exp(errr, erri)));


/*
uint64_t ct3[2][L][N];
uint64_t ct4[2][L - 1][N];

encode<N, logN, L>(z1r, z1i, Delta, q, pt);
enc<N, L>(pt, s, q, ct3);
encode<N, logN, L-1>(z2r, z2i, (uint64_t)((long double)Delta* ((long double)Delta / (long double)q[L - 1])), q, pt);
enc<N, L - 1>(pt, s, q, ct4);
printf("Encryption complete\n");
printf("Delta %lld\n", Delta);
printf("Delta %lld\n", (uint64_t)((long double)Delta * ((long double)Delta / (long double)q[L - 1])));

uint64_t ct3_RS[2][L - 1][N];
uint64_t ct3_RS_cache[2][L][N];
for (int k = 0; k < 2; ++k) {
	for (int i = 0; i < L; i++) {
		for (int j = 0; j < N; j++) {
			ct3_RS_cache[k][i][j] = mul_mod(Delta, ct3[k][i][j], q[i]);
		}
	}
}
RS_hat<N, L, L - 1>(q, ct3_RS_cache, ct3_RS);

uint64_t ct_add[2][L - 1][N];
for (int i = 0; i < L - 1; i++) {
	for (int j = 0; j < N; j++) {
		ct_add[0][i][j] = ct3_RS[0][i][j] + ct4[0][i][j];
		ct_add[1][i][j] = ct3_RS[1][i][j] + ct4[1][i][j];
	}
}


dec<N, L - 1>(ct_add, s, q, pt_out);
decode<N, logN, L - 1>(pt_out, (uint64_t)((long double)Delta * ((long double)Delta / (long double)q[L - 1])), q, wr, wi);
//	dec<N, L>(ct3, s, q, pt_out);
//	decode<N, logN, L>(pt_out, (uint64_t)(Delta), q, wr, wi);
	printf("Decryption complete\n");
	for (int i = 0; i < N / 2; i++) {
		errr[i] = (wr[i] - (z1r[i] + z2r[i])) / wr[i];
		erri[i] = 0;
	}

	printf("error in Linf norm: %lf\n", -std::log2l(norm_max(errr, erri)));
	printf("error in L2 norm: %lf\n", -std::log2l(sqrt(norm_square(errr, erri))));
	printf("error in L2_exp norm: %lf\n", -std::log2l(norm_square_exp(errr, erri)));
	*/















