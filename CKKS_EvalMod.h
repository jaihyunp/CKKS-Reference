#pragma once

#include "CKKS_ks_wih_dnum.h"
#include "CKKS_lineartransform.h"
#include "CKKS_poly.h"



template<int N, int L, int DNUM, int K>
void EvalMod_h_192_mul_q_over_Delta(const uint64_t q[L],
	const uint64_t p[K], uint64_t Delta,
	const uint64_t evk_hat[DNUM][2][DNUM * K + K][N],
	const uint64_t         ct_hat[2][L][N],
	uint64_t ct_evalmod_hat[2][L - 9][N], const int s[N], double tr[N / 2]) {
	
	double c[] = { -2.84529398747917162409e-18, 1.40829338704334061649e-01, -7.58941316265013770449e-17, -7.19566808032123528260e-02, -1.35521334298163801589e-16, 3.14734688318025090936e-01, -3.25449777990325958176e-17, -6.14678201819273239970e-01, -2.90528692129655971332e-17, 1.74769865183753869697e-01, 2.37412526119722306076e-17, 3.73430130070927512875e-02, 9.64574127417459471246e-17, -8.31909520306124417033e-02, -1.24518166394120163836e-16, 1.28405049044885422038e-01, 1.99648137023739073036e-18, 2.47970475294174041991e-01, -3.36700298340906705129e-17, -1.95489024287215262810e-01, -5.45348155061956594221e-17, 7.84513028613951113321e-01, -6.74946857968754315119e-17, -9.95925646437272993339e-01, 1.40467699145745675233e-17, 3.36175543040554858365e-01, -1.54650671797943756927e-16, -2.46866274710532740411e-01, 1.25269827167496108407e-16, 8.30328267314629081541e-02, 4.58101613769693512443e-16, -4.04302682927135115243e-02, 1.26155831044648452197e-17, 1.10368985456504647633e-03, -2.82748176353342681581e-17, -3.85726275603731595058e-04, -1.60581942887265385299e-16, 6.05169485509866731570e-05, -9.54506786731657775881e-17, -1.60357755525372959698e-05, -1.39830674553646880860e-17, 9.60231783399508091276e-07, 9.36089133798698824026e-17, -2.00419728640990657501e-07, -1.16196001303934638819e-17, 1.90463988603809574087e-08, -1.91246411704479019825e-16, -3.22204548202335592976e-09, 4.16190818855493052269e-17, 6.24875008270227891354e-11, 5.36960222705191952928e-17, -8.75809220492681329210e-12, -1.35534452846880670782e-16, 5.66502801344388376973e-13, -1.39400864208274666098e-16, -6.72225532112069847271e-14, -4.16946319642174047022e-17, 1.36132348366319765563e-15, -9.04945035868366645017e-17, 6.37864735960474548804e-16, 2.97281478745637994436e-16, -6.27223658499345852726e-16, 4.39184854924537499461e-16, 3.2709827308e-16 };
	double d[] = { 1.11967834533887089510e-01, 1.15836602658729409747e-15, 9.82830818362460556514e-01, -1.66309732631114582566e-15, 4.49882806538154877973e-01, 1.32216460780035379061e-15, -1.29917408307411719193e-02, 1.52915675194968074776e-16, -1.21563801272885629867e-01, 2.45324154487675598504e-15, 2.52472893007162380030e+00, -3.26912879873184984163e-15, -5.02140286840478733410e-01, -2.83592857616607464105e-16, -2.84037651199109086875e+00, -8.11028461236646722445e-17, -7.92066911484812508082e-02, 1.67480302441616736356e-15, 1.20504955203691865862e+00, -5.62085683739023608564e-16, -3.08348186415920988424e-02, 2.24600871629064016135e-15, -8.63537922882909869671e-01, -2.48994041178130976316e-15, 8.19094909701504780841e-01, 7.12733465748659975931e-16, -7.95885607817973017575e-01, 1.28478934022931736697e-16, 3.46134301386879117413e-01, 2.30076149981798539748e-15, -1.97261955893251239580e-01, -3.98734662154818075817e-15, 6.23485552935949295661e-03, -6.00272811362794957481e-16, -2.43177739966620462889e-03, 5.46074525951646090626e-18, 4.22796397012973625661e-04, 3.70965830166666639534e-16, -1.21980156047958579663e-04, 1.43615948318353411476e-15, 7.91333801993830504801e-06, -4.58548181077484460683e-16, -1.77111843399497543083e-06, -3.79507333298473781701e-17, 1.79809440840117585223e-07, -3.61135742443564421965e-15, -3.22822449759573795758e-08, 2.54690087496140642984e-15, 6.62507027153921759912e-10, -7.83840637892694580309e-16, -9.78188038965720494429e-11, -4.29183062916307089207e-16, 6.65024020537251935677e-12, 1.91363411765791569904e-17, -8.21376829120915831134e-13, 2.05047750643674072802e-15, 2.38522793010224241371e-14, -1.14982385495602525629e-15, -3.59146446060181189964e-15, 9.04493085252262006972e-16, -2.37375504008781032242e-15, -1.37029353660641094982e-17, 1.45096501114769884061e-16, -3.1651415092e-17 };

	uint64_t A1_hat[2][L - 6][N];
	uint64_t B1_hat[2][L - 6][N];

	uint64_t ct_T_32 [2][L - 5][N];
	EvalPoly<N, L, DNUM, K, 6>(q, p, evk_hat, ct_hat, Delta, Delta, c, A1_hat, ct_T_32, 1);
	EvalPoly<N, L, DNUM, K, 6>(q, p, evk_hat, ct_hat, Delta, Delta, d, B1_hat, ct_T_32, 2);

	
	

	uint64_t A2_hat[2][L - 7][N];
	uint64_t B2_hat[2][L - 7][N];
	mul_rs<N, L - 6, DNUM, K>(q, p, evk_hat, A1_hat, B1_hat, A2_hat);
	mul_rs<N, L - 6, DNUM, K>(q, p, evk_hat, B1_hat, B1_hat, B2_hat);


	for (int i = 0; i < 2; i++)
		for (int j = 0; j < L - 7; j++)
			for (int k = 0; k < N; k++) {
				A2_hat[i][j][k] = (2 * A2_hat[i][j][k]) % q[j];
				B2_hat[i][j][k] = (i == 0) ? ((2 * B2_hat[i][j][k] + q[j] - (Delta % q[j])) % q[j])
					: ((2 * B2_hat[i][j][k]) % q[j]);
			}
	
	uint64_t A3_hat[2][L - 8][N];
	uint64_t B3_hat[2][L - 8][N];
	mul_rs<N, L - 7, DNUM, K>(q, p, evk_hat, A2_hat, B2_hat, A3_hat);
	mul_rs<N, L - 7, DNUM, K>(q, p, evk_hat, B2_hat, B2_hat, B3_hat);
	for (int i = 0; i < 2; i++)
		for (int j = 0; j < L - 8; j++)
			for (int k = 0; k < N; k++) {
				A3_hat[i][j][k] = (4 * A3_hat[i][j][k]) % q[j];
				B3_hat[i][j][k] = (i == 0) ? ((2 * B3_hat[i][j][k] + q[j] - (Delta % q[j])) % q[j])
					: ((2 * B3_hat[i][j][k]) % q[j]);
			}

	

	mul_rs<N, L - 8, DNUM, K>(q, p, evk_hat, A3_hat, B3_hat, ct_evalmod_hat);

}
