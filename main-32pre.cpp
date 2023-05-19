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
#define L 50
#define DNUM 5
#define K (L/DNUM)
#define Delta (uint64_t) (1ULL << 28)


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

		// for delta = 1<<31
		// uint64_t q32[60] = { 2147352577ULL, 2146959361ULL, 2148794369ULL, 2146041857ULL, 2151677953ULL, 2144468993ULL, 2152071169ULL, 2150760449ULL, 2142502913ULL, 2150891521ULL, 2135818241ULL, 2158231553ULL, 2153250817ULL, 2135162881ULL, 2154823681ULL, 2135031809ULL, 2161508353ULL, 2155216897ULL, 2130706433ULL, 2155610113ULL, 2132279297ULL, 2156265473ULL, 2126118913ULL, 2134638593ULL, 2175926273ULL, 2130444289ULL, 2178285569ULL, 2165702657ULL, 2112225281ULL, 2128740353ULL, 2182742017ULL, 2125725697ULL, 2185756673ULL, 2166620161ULL, 2109603841ULL, 2167013377ULL, 2113011713ULL, 2168455169ULL, 2106458113ULL, 2121793537ULL, 2194145281ULL, 2171207681ULL, 2102788097ULL, 2171600897ULL, 2099249153ULL, 2172518401ULL, 2100953089ULL, 2120482817ULL, 2204368897ULL, 2119303169ULL, 2203058177ULL, 2183135233ULL, 2082078721ULL, 2117468161ULL, 2211053569ULL, 2114977793ULL, 2214854657ULL, 2184183809ULL, 2080505857ULL, 2114191361ULL };

		// for delta = 1<<30
		//uint64_t q32[60] = { 1073872897ULL, 1073479681ULL, 1074266113ULL, 1071513601ULL, 1081212929ULL, 1062862849ULL, 1070727169ULL, 1083703297ULL, 1068236801ULL, 1085276161ULL, 1065484289ULL, 1086455809ULL, 1053818881ULL, 1064697857ULL, 1089601537ULL, 1062469633ULL, 1095761921ULL, 1087635457ULL, 1052508161ULL, 1060765697ULL, 1100873729ULL, 1088684033ULL, 1042415617ULL, 1089208321ULL, 1045430273ULL, 1091043329ULL, 1040056321ULL, 1091174401ULL, 1040842753ULL, 1056440321ULL, 1107296257ULL, 1056178177ULL, 1110048769ULL, 1055260673ULL, 1108738049ULL, 1092616193ULL, 1038745601ULL, 1093533697ULL, 1036779521ULL, 1054212097ULL, 1111883777ULL, 1051721729ULL, 1115815937ULL, 1093795841ULL, 1032192001ULL, 1093926913ULL, 1037303809ULL, 1094582273ULL, 1043464193ULL, 1049100289ULL, 1103626241ULL, 1099431937ULL, 1034813441ULL, 1102053377ULL, 1031667713ULL, 1113980929ULL, 1010565121ULL, 1048707073ULL, 1125646337ULL, 1116209153ULL };
		
		// for delta = 1<<29
		//uint64_t q32[60] = { 537133057ULL, 536215553ULL, 536608769ULL, 539754497ULL, 533463041ULL, 540540929ULL, 528351233ULL, 540672001ULL, 529924097ULL, 541327361ULL, 531628033ULL, 532283393ULL, 543293441ULL, 543031297ULL, 529268737ULL, 547749889ULL, 517472257ULL, 548012033ULL, 517079041ULL, 526123009ULL, 558235649ULL, 549978113ULL, 518520833ULL, 525729793ULL, 555220993ULL, 524943361ULL, 557187073ULL, 550371329ULL, 518914049ULL, 524812289ULL, 553254913ULL, 551288833ULL, 510001153ULL, 521011201ULL, 566886401ULL, 552861697ULL, 508690433ULL, 512360449ULL, 582746113ULL, 561774593ULL, 489422849ULL, 564658177ULL, 488243201ULL, 511967233ULL, 591265793ULL, 568066049ULL, 478937089ULL, 569638913ULL, 483655681ULL, 570163201ULL, 475267073ULL, 570949633ULL, 487063553ULL, 572915713ULL, 471072769ULL, 575275009ULL, 469762049ULL, 576454657ULL, 473694209ULL, 504496129ULL };

		// for delta = 1<<28
		uint64_t q32[60] = {268042241ULL, 269221889ULL, 270532609ULL, 264634369ULL, 265420801ULL, 272760833ULL, 270794753ULL, 263454721ULL, 274726913ULL, 263323649ULL, 276037633ULL, 256376833ULL, 261881857ULL, 279838721ULL, 261488641ULL, 277086209ULL, 276430849ULL, 256770049ULL, 260702209ULL, 281935873ULL, 260571137ULL, 285474817ULL, 258605057ULL, 292159489ULL, 283508737ULL, 244842497ULL, 257949697ULL, 284950529ULL, 254279681ULL, 302776321ULL, 253493249ULL, 302252033ULL, 253100033ULL, 297664513ULL, 288882689ULL, 234356737ULL, 290455553ULL, 228720641ULL, 291373057ULL, 230686721ULL, 249561089ULL, 309460993ULL, 295305217ULL, 223215617ULL, 246415361ULL, 319291393ULL, 245760001ULL, 319160321ULL, 299499521ULL, 218628097ULL, 245235713ULL, 323092481ULL, 244973569ULL, 316538881ULL, 304218113ULL, 211025921ULL, 241827841ULL, 330301441ULL, 240648193ULL, 332660737ULL};

		for (int i = 0; i < K; i++)
			p[i] = q32[L+K-1-i];
		for (int i = 0; i < L; i++)
			q[i] = q32[L-1-i];

	}
	printf("32 -bit RNS prime set\n");
	printf("Scaling factor: %lld\n", Delta);

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















