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

    uint64_t q[L];   uint64_t p[K];
    {
        uint64_t q30[60] = { 537133057ULL, 536608769ULL, 1094582273ULL, 263323649ULL, 1137049601ULL, 253493249ULL, 1302724609ULL, 221249537ULL, 903086081ULL, 319160321ULL, 1133510657ULL, 254279681ULL, 1172832257ULL, 245760001ULL, 1536688129ULL, 187564033ULL, 1480851457ULL, 194641921ULL, 1841692673ULL, 156499969ULL, 735969281ULL, 391643137ULL, 1122500609ULL, 256770049ULL, 778436609ULL, 370278401ULL, 801374209ULL, 359661569ULL, 1905524737ULL, 151257089ULL, 1520828417ULL, 189530113ULL, 611844097ULL, 471072769ULL, 1089208321ULL, 264634369ULL, 820248577ULL, 351404033ULL, 1386479617ULL, 207880193ULL, 935329793ULL, 308150273ULL, 1474953217ULL, 195428353ULL, 867434497ULL, 332267521ULL, 553254913ULL, 521011201ULL, 769130497ULL, 374734849ULL, 1997406209ULL, 144310273ULL, 1191837697ULL, 241827841ULL, 1260257281ULL, 228720641ULL, 1981022209ULL, 145489921ULL, 1677459457ULL, 171835393ULL};

        for (int i = 0; i < L; i++)
            q[i] = q30[L-1-i];
        for (int i = 0; i < K; i++)
            p[i] = q30[L+K-1-i];

    }
    printf("30 -bit RNS prime set\n");

    int h = 192; int s[N];
    keygen<N>(h, s);

    //---------------------------------------------------------------------
    // swkgen
    //---------------------------------------------------------------------
 //   SparseComplexMatrix<N / 2, 3> E[logN - 1];
 //   splitU0R<logN>(E);
 //   uint64_t rkey_hat[logN - 1][3][DNUM][2][DNUM * K + K][N];
 //   for (int n = 0; n < logN - 1; n++)
 //      rkey_gen<N, L, DNUM, K, 3>(E[n], q, p, s, rkey_hat[n]);

    uint64_t  evk_hat[DNUM][2][DNUM * K + K][N];
    //uint64_t ckey_hat[DNUM][2][DNUM * K + K][N];
    int ss[N]; conv<N>(s, s, ss);
    swkgen< N, L, DNUM>(ss, s, q, p, evk_hat);
    //   int sconj[N]; conj<N>(s, sconj);
    //   swkgen<N, L, DNUM>(sconj, s, q, p, ckey_hat);

    printf("KeyGen complete\n");


    // See CKKS_basic.h
    //---------------------------------------------------------------------
    // Ecd & Enc
    //---------------------------------------------------------------------
    uint64_t ct1[2][L][N];
    uint64_t ct2[2][L][N];
    uint64_t ct_tmp[2][L][N];
    uint64_t ct_tmp2[2][L - 1][N];
    uint64_t ct_out[2][L - 2][N];
    uint64_t ct6[2][6][N];
    uint64_t ct_bot[2][4][N];
    uint64_t pt[L][N];
    uint64_t pt_out[L][N];
    uint64_t Delta_out = Delta;
    double z1r[N / 2], z1i[N / 2], z2r[N / 2], z2i[N / 2];
    double wr[N / 2], wi[N / 2];
    double errr[N / 2], erri[N / 2];

    for (int i = 0; i < N / 2; i++) {
        z1r[i] = ((double)rand()) / RAND_MAX / sqrt(2);
        z1i[i] = 0;
        z2r[i] = ((double)rand()) / RAND_MAX / sqrt(2);
        //z2r[i] = 0;
        z2i[i] = 0;
    }



    /* Enc Dec */
    printf("\nEnc & Dec Error\n");

    encode<N, logN, L>(z1r, z1i, Delta, q, pt);
    enc<N, L>(pt, s, q, ct1);

    dec<N, L>(ct1, s, q, pt_out);
    decode<N, logN, L>(pt_out, Delta, q, wr, wi);

    for (int i = 0; i < N / 2; i++) {
        errr[i] = wr[i] - z1r[i];
        erri[i] = 0;
    }
    printf("error in Linf norm: %lf\n", -std::log2l(norm_max(errr, erri)));
    printf("error in L2 norm: %lf\n", -std::log2l(sqrt(norm_square(errr, erri))));
    printf("error in L2_exp norm: %lf\n", -std::log2l(norm_square_exp(errr, erri)));



    /* Add */
    printf("\nAdd error\n");
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
    printf("error in Linf norm: %lf\n", -std::log2l(norm_max(errr, erri)));
    printf("error in L2 norm: %lf\n", -std::log2l(sqrt(norm_square(errr, erri))));
    printf("error in L2_exp norm: %lf\n", -std::log2l(norm_square_exp(errr, erri)));



    /* Rescale w/ a const scaling factor */
    printf("\nRescale (single) error\n");
    printf("* With a constant scaling factor\n");
    encode<N, logN, L>(z1r, z1i, Delta, q, pt);
    enc<N, L>(pt, s, q, ct1);

    for (int k = 0; k < 2; ++k) {
        for (int i = 0; i < L - 2; i++) {
            for (int j = 0; j < N; j++) {
                ct_out[k][i][j] = ct1[k][i][j];
            }
        }
    }

    Delta_out = Delta;
    dec<N, L - 2>(ct_out, s, q, pt_out);
    decode<N, logN, L - 2>(pt_out, Delta_out, q, wr, wi);

    for (int i = 0; i < N / 2; i++) {
        errr[i] = wr[i] - z1r[i];
        erri[i] = 0;
    }
    printf("error in Linf norm: %lf\n", -std::log2l(norm_max(errr, erri)));
    printf("error in L2 norm: %lf\n", -std::log2l(sqrt(norm_square(errr, erri))));
    printf("error in L2_exp norm: %lf\n", -std::log2l(norm_square_exp(errr, erri)));



    /* Rescale w/ leveled scaling factors */
    printf("\nRescale (single) error\n");
    printf("* With leveled scaling factors\n");
    encode<N, logN, L>(z1r, z1i, Delta, q, pt);
    enc<N, L>(pt, s, q, ct1);
    for (int k = 0; k < 2; ++k) {
        for (int i = 0; i < L; i++) {
            for (int j = 0; j < N; j++) {
                ct_tmp[k][i][j] = mul_mod(Delta, ct1[k][i][j], q[i]);
            }
        }
    }
    RS_hat<N, L, L - 2>(q, ct_tmp, ct_out);
    Delta_out = Delta;
    Delta_out = (uint64_t)((long double)Delta_out * ((long double)Delta_out / (long double)q[L - 1] / (long double)q[L - 2]));
    dec<N, L - 2>(ct_out, s, q, pt_out);
    decode<N, logN, L - 2>(pt_out, Delta_out, q, wr, wi);

    for (int i = 0; i < N / 2; i++) {
        errr[i] = wr[i] - z1r[i];
        erri[i] = 0;
    }
    printf("error in Linf norm: %lf\n", -std::log2l(norm_max(errr, erri)));
    printf("error in L2 norm: %lf\n", -std::log2l(sqrt(norm_square(errr, erri))));
    printf("error in L2_exp norm: %lf\n", -std::log2l(norm_square_exp(errr, erri)));



    /* Mult w/ a const scaling factor */
    printf("\nMult once error\n");
    printf("* With a constant scaling factor\n");
    encode<N, logN, L>(z1r, z1i, Delta, q, pt);
    enc<N, L>(pt, s, q, ct1);
    encode<N, logN, L>(z2r, z2i, Delta, q, pt);
    enc<N, L>(pt, s, q, ct2);

    mul<N, L, DNUM, K>(q, p, evk_hat, ct1, ct2, ct_tmp);
    RS_hat<N, L, L - 2>(q, ct_tmp, ct_out);

    Delta_out = Delta;
    dec<N, L - 2>(ct_out, s, q, pt_out);
    decode<N, logN, L - 2>(pt_out, Delta_out, q, wr, wi);


    for (int i = 0; i < N / 2; i++) {
        errr[i] = wr[i] - z1r[i] * z2r[i];
        erri[i] = 0;
    }
    printf("error in Linf norm: %lf\n", -std::log2l(norm_max(errr, erri)));
    printf("error in L2 norm: %lf\n", -std::log2l(sqrt(norm_square(errr, erri))));
    printf("error in L2_exp norm: %lf\n", -std::log2l(norm_square_exp(errr, erri)));



    /* Mult w/ leveled scaling factors */
    printf("\nMult once error\n");
    printf("* With leveled scaling factors\n");
    encode<N, logN, L>(z1r, z1i, Delta, q, pt);
    enc<N, L>(pt, s, q, ct1);
    encode<N, logN, L>(z2r, z2i, Delta, q, pt);
    enc<N, L>(pt, s, q, ct2);

    mul<N, L, DNUM, K>(q, p, evk_hat, ct1, ct2, ct_tmp);
    RS_hat<N, L, L - 2>(q, ct_tmp, ct_out);

    Delta_out = Delta;
    Delta_out = (uint64_t)((long double)Delta_out * ((long double)Delta_out / (long double)q[L - 1] / (long double)q[L - 2]));
    dec<N, L - 2>(ct_out, s, q, pt_out);
    decode<N, logN, L - 2>(pt_out, Delta_out, q, wr, wi);


    for (int i = 0; i < N / 2; i++) {
        errr[i] = wr[i] - z1r[i] * z2r[i];
        erri[i] = 0;
    }
    printf("error in Linf norm: %lf\n", -std::log2l(norm_max(errr, erri)));
    printf("error in L2 norm: %lf\n", -std::log2l(sqrt(norm_square(errr, erri))));
    printf("error in L2_exp norm: %lf\n", -std::log2l(norm_square_exp(errr, erri)));



    /* Rescales w/ constant scaling factors */
    printf("\nRescale (multiple) error\n");
    printf("* With a constant scaling factor\n");
    encode<N, logN, L>(z1r, z1i, Delta, q, pt);
    enc<N, L>(pt, s, q, ct1);
    //printf("Encryption complete\n");

    for (int k = 0; k < 2; ++k) {
        for (int i = 0; i < 4; i++) {
            for (int j = 0; j < N; j++) {
                ct_bot[k][i][j] = ct1[k][i][j];
            }
        }
    }

    Delta_out = Delta;
    dec<N, 4>(ct_bot, s, q, pt_out);
    decode<N, logN, 4>(pt_out, Delta_out, q, wr, wi);

    for (int i = 0; i < N / 2; i++) {
        errr[i] = wr[i] - z1r[i];
        erri[i] = 0;
    }
    printf("error in Linf norm: %lf\n", -std::log2l(norm_max(errr, erri)));
    printf("error in L2 norm: %lf\n", -std::log2l(sqrt(norm_square(errr, erri))));
    printf("error in L2_exp norm: %lf\n", -std::log2l(norm_square_exp(errr, erri)));



    /* Rescales w/ leveled scaling factors */
    printf("\nRescale (multiple) error\n");
    printf("* With leveled scaling factors\n");
    encode<N, logN, L>(z1r, z1i, Delta, q, pt);
    enc<N, L>(pt, s, q, ct1);
    //printf("Encryption complete\n");


    Delta_out = Delta;
    for (int i = L - 1; i > 4; i -= 2) {
        Delta_out = (uint64_t)((long double)Delta_out * ((long double)Delta_out / (long double)q[i] / (long double)q[i - 1]));
    }
    uint64_t Delta_diff = (uint64_t)((long double)Delta_out * ((long double)q[4] * q[5] / (long double)Delta));

    for (int k = 0; k < 2; ++k) {
        for (int i = 0; i < 6; i++) {
            for (int j = 0; j < N; j++) {
                ct6[k][i][j] = mul_mod(Delta_diff, ct1[k][i][j], q[i]);
            }
        }
    }
    RS_hat<N, 6, 4>(q, ct6, ct_bot);

    dec<N, 4>(ct_bot, s, q, pt_out);
    decode<N, logN, 4>(pt_out, Delta_out, q, wr, wi);

    for (int i = 0; i < N / 2; i++) {
        errr[i] = wr[i] - z1r[i];
        erri[i] = 0;
    }
    printf("error in Linf norm: %lf\n", -std::log2l(norm_max(errr, erri)));
    printf("error in L2 norm: %lf\n", -std::log2l(sqrt(norm_square(errr, erri))));
    printf("error in L2_exp norm: %lf\n", -std::log2l(norm_square_exp(errr, erri)));




}














