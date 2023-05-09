#pragma once
#include <stdint.h>
#include <inttypes.h>
#include <cassert>
#include "DFT.h"
#include "NumberTheory.h"
#include "BigInt.h"
#include "CKKS_ks_wih_dnum.h"

//
template<int N, int logN, int L>
void encode( const double zr[N/2],
	         const double zi[N/2], uint64_t Delta,
			 const uint64_t q[L], uint64_t pt[L][N]){
	double m[N]; ifft<N,logN>(zr,zi,m);
	for(int i=0;i<N;i++){
		uint64_t pti = uint64_t((m[i]>0?m[i]:-m[i])*Delta+0.5);
		for(int j=0;j<L;j++){
			pt[j][i] = pti%q[j];
			if(m[i]<0) pt[j][i]=q[j]-pt[j][i];
		}
	}
}

//
template<int N, int logN, int L>
void decode( const uint64_t pt[L][N], uint64_t Delta,
	         const uint64_t  q[L], double zr[N/2], double zi[N/2]){
	//
	BigInt Q(L), Qhalf(L); Q[0]=1;
	for(int i=0;i<L;i++){ BigInt temp(Q); temp.mul(q[i],Q); }
	{uint64_t r; Q.div(2,Qhalf,r); }
	double m[N];
	for(int j=0;j<N;j++){ 
		uint64_t a[L]; 
		for(int i=0;i<L;i++) 
			a[i]=pt[i][j];
		BigInt A; 
		icrt<L>(q,a,A);
		if(A>=Qhalf){ 
			BigInt temp(A); 
			Q.sub(temp,A); 
			m[j]=-A.to_real()/Delta;
		}
		else                                         
			m[j]= A.to_real()/Delta;
	}
	fft<N,logN>(m,zr,zi);
}

//
template<int N, int logN, int L>
void decode( const uint64_t pt[L][N], double Delta,
	         const uint64_t  q[L], double zr[N/2], double zi[N/2]){
	//
	BigInt Q(L), Qhalf(L); Q[0]=1;
	for(int i=0;i<L;i++){ BigInt temp(Q); temp.mul(q[i],Q); }
	{uint64_t r; Q.div(2,Qhalf,r); }
	double m[N];
	for(int j=0;j<N;j++){ 
		uint64_t a[L]; 
		for(int i=0;i<L;i++) 
			a[i]=pt[i][j];
		BigInt A; 
		icrt<L>(q,a,A);
		if(A>=Qhalf){ 
			BigInt temp(A); 
			Q.sub(temp,A); 
			m[j]=-A.to_real()/Delta;
		}
		else                                         
			m[j]= A.to_real()/Delta;
	}
	fft<N,logN>(m,zr,zi);
}

//
template<int N>
void keygen(int h, int s[N]){
	for(int i=0;i<N;i++) s[i]=0;
	for(int i=0;i<h;i++){
		int j=rand()%N;
		while(s[j]!=0)
			j=rand()%N;
		s[j]=(rand()%2==0)?1:-1;
	}
}

//
template<int N, int L>
void enc( const uint64_t pt[L][N], const int s[N],
	      const uint64_t q[L], uint64_t cthat[2][L][N]){
	int e[N]; for(int i=0;i<N;i++) e[i]=(rand()%17)-8;
	for(int j=0;j<L;j++){
		uint64_t shat[N];
		for(int i=0;i<N;i++){
			cthat[1][j][i]=rand_uint64()%q[j];
			shat      [i]=(  q[j]+s[i])%q[j];
			cthat[0][j][i]=(q[j]+e[i]+pt[j][i])%q[j];
		}
		ntt<N>(cthat[0][j],q[j]);
		ntt<N>(cthat[1][j],q[j]);
		// debug
		{
			intt<N>(cthat[1][j],q[j]);
			 ntt<N>(cthat[1][j],q[j]);
		}
		ntt<N>( shat      ,q[j]);
		for(int i=0;i<N;i++)
			cthat[0][j][i]=(cthat[0][j][i]+q[j]-mul_mod(shat[i],cthat[1][j][i],q[j]))%q[j];
	}
}

//
template<int N, int L>
void dec( const uint64_t cthat[2][L][N], const int s[N],
	      const uint64_t q[L], uint64_t pt[L][N]){
	for(int j=0;j<L;j++){
		uint64_t shat[N];
		for(int i=0;i<N;i++)
			shat[i]=(q[j]+s[i])%q[j];
		ntt<N>(shat,q[j]);

		for(int i=0;i<N;i++)
			pt[j][i]=(cthat[0][j][i]+mul_mod(cthat[1][j][i],shat[i],q[j]))%q[j];
		intt<N>(pt[j],q[j]);
	}
}

//
template<int N, int L>
void RS( const uint64_t q[L], const uint64_t pt[L][N], uint64_t pt_[L-1][N] ){
	uint64_t qinv[L-1];
	for(int j=0;j<L-1;j++)
		qinv[j]=inv_mod(q[L-1],q[j]);
	for(int i=0;i<N;i++){
		uint64_t r = pt[L-1][i]>(q[L-1]/2)?1:0;
		for(int j=0;j<L-1;j++)
			pt_[j][i]=(mul_mod((pt[j][i]+q[j]-(pt[L-1][i]%q[j]))%q[j],qinv[j],q[j])+r)%q[j];
	}
}

// ADDED by KSH 1024
// RS pt Lfr -> Lto
template<int N, int Lfr, int Lto>
void RS( const uint64_t q[Lfr], const uint64_t pt_orig[Lfr][N], uint64_t pt_[Lto][N] ){
	assert(Lfr > Lto);
	uint64_t pt[Lfr][N];
	for(int i=0; i<Lfr; ++i)
	for(int j=0; j<N; ++j)
		pt[i][j] = pt_orig[i][j];
	for(int k=Lfr-1;k>=Lto;--k){
		uint64_t qinv[Lfr];
		for(int j=0;j<k;++j)
			qinv[j]=inv_mod(q[k],q[j]);
		for(int i=0;i<N;++i){
			uint64_t r = (pt[k][i]>q[k]/2)?1:0;
			for(int j=0;j<k;++j){
				pt[j][i]=mul_mod((pt[j][i]+q[j]-(pt[k][i]%q[j]))%q[j], qinv[j], q[j]);
				if(k == Lto)
					pt_[j][i]=(pt[j][i]+r)%q[j];
			}
		}
	}
}

// MOVED FROM NumberTheory.h 1024
template<int N, int L>
void RS_hat(const uint64_t q[L],
			const uint64_t pt_hat [L  ][N],
				  uint64_t pt_hat_[L-1][N]){
	uint64_t temp[L][N];
	for(int i=0;i<L;i++)
	for(int j=0;j<N;j++)
		temp[i][j] = pt_hat[i][j];

	intt<N,L>(q, temp);
	RS<N,L>(q, temp, pt_hat_);
	ntt<N,L-1>(q, pt_hat_);
}

// ADDED by KSH 1024
// RS_hat ct Lfr -> Lto
template<int N, int Lfr, int Lto>
void RS_hat( const uint64_t q[Lfr], const uint64_t ct_hat[2][Lfr][N], uint64_t ct_hat_[2][Lto][N] ){
	uint64_t ct[2][Lfr][N];
	for(int i = 0; i < 2; ++i) {
		for(int j = 0; j < Lfr; ++j) {
			for(int k = 0; k < N; ++k) {
				ct[i][j][k] = ct_hat[i][j][k];
			}
		}
	}
	intt<N, Lfr>(q, ct[0]);
	intt<N, Lfr>(q, ct[1]);
	RS<N, Lfr, Lto>(q, ct[0], ct_hat_[0]);
	RS<N, Lfr, Lto>(q, ct[1], ct_hat_[1]);
	ntt<N, Lto>(q, ct_hat_[0]);
	ntt<N, Lto>(q, ct_hat_[1]);
}

// 
template<int N>
void conv( const int A[N],
	       const int B[N], int C[N] ){

	for(int i=0; i<N; i++){
		C[i] = 0;
		for(int j=0  ; j<=i; j++) C[i] += A[j]*B[i-j];
		for(int j=i+1; j<N ; j++) C[i] -= A[j]*B[i-j+N];
	}
}

//
template<int N>
void conv(const double a[N],
		  const double b[N], double c[N]){
	for(int i=0; i<N; i++){
		c[i]=0;
		for(int j=0; j<i+1; j++)
			c[i] += a[j]*b[i-j];
		for(int j=i+1; j<N; j++)
			c[i] -= a[j]*b[i+N-j];
	}
}

//
template<int N, int L>
void keygen(const int s[N], const uint64_t q[L],
								  uint64_t pk_hat[2][L][N]){
	uint64_t zero[L][N];
	for(int i=0; i<L; i++)
	for(int j=0; j<N; j++)
		zero[i][j] = 0;
	enc<N,L>(zero, s, q, pk_hat);
}

//
template<int N, int L>
void enc(const uint64_t pt[L][N],
		 const uint64_t pk_hat[2][L][N], const uint64_t q[L],
			   uint64_t ct_hat[2][L][N]){
	int e0[N], e1[N], v[N];
	for(int i=0; i<N; i++){
		e0[i] = (rand() % 17) - 8;
		e1[i] = (rand() % 17) - 8;
		int r = rand()%4;
		 v[i] = (r==0) ? (-1) : ((r==1) ? 1 : 0);
	}

	for(int j=0; j<L; j++){
		uint64_t v_hat[N], pt_j_hat[N], e0_hat[N], e1_hat[N];
		for(int i=0; i<N; i++){
			v_hat[i]    = (q[j]+v[i]) % q[j];
			pt_j_hat[i] = (q[j]+pt[j][i]) % q[j];
			e0_hat[i]   = (q[j]+e0[i]) % q[j];
			e1_hat[i]   = (q[j]+e1[i]) % q[j];
		}

		ntt<N>(   v_hat,q[j]);
		ntt<N>(pt_j_hat,q[j]);
		ntt<N>(  e0_hat,q[j]);
		ntt<N>(  e1_hat,q[j]);

		for(int i=0; i<N; i++){
			ct_hat[0][j][i] = (mul_mod(v_hat[i], pk_hat[0][j][i], q[j]) + pt_j_hat[i] + e0_hat[i]) % q[j];
			ct_hat[1][j][i] = (mul_mod(v_hat[i], pk_hat[1][j][i], q[j])               + e1_hat[i]) % q[j];
		}
	}
}



/***********************************************/
/***** For precision of higher than 64 bits ****/
/***********************************************/

/*
//
template<int N, int logN, int L>
void encode(const double zr[N / 2],
	const double zi[N / 2], uint64_t Delta,
	const uint64_t q[L], uint64_t pt[L][N]) {
	double m[N]; ifft<N, logN>(zr, zi, m);
	for (int i = 0;i < N;i++) {
		uint64_t pti = uint64_t((m[i] > 0 ? m[i] : -m[i]) * Delta + 0.5);
		for (int j = 0;j < L;j++) {
			pt[j][i] = pti % q[j];
			if (m[i] < 0) pt[j][i] = q[j] - pt[j][i];
		}
	}
}

//
template<int N, int logN, int L>
void decode(const uint64_t pt[L][N], uint64_t Delta,
	const uint64_t  q[L], double zr[N / 2], double zi[N / 2]) {
	//
	BigInt Q(L), Qhalf(L); Q[0] = 1;
	for (int i = 0;i < L;i++) { BigInt temp(Q); temp.mul(q[i], Q); }
	{uint64_t r; Q.div(2, Qhalf, r); }
	double m[N];
	for (int j = 0;j < N;j++) {
		uint64_t a[L];
		for (int i = 0;i < L;i++)
			a[i] = pt[i][j];
		BigInt A;
		icrt<L>(q, a, A);
		if (A >= Qhalf) {
			BigInt temp(A);
			Q.sub(temp, A);
			m[j] = -A.to_real() / Delta;
		}
		else
			m[j] = A.to_real() / Delta;
	}
	fft<N, logN>(m, zr, zi);
}

//
template<int N, int logN, int L>
void decode(const uint64_t pt[L][N], double Delta,
	const uint64_t  q[L], double zr[N / 2], double zi[N / 2]) {
	//
	BigInt Q(L), Qhalf(L); Q[0] = 1;
	for (int i = 0;i < L;i++) { BigInt temp(Q); temp.mul(q[i], Q); }
	{uint64_t r; Q.div(2, Qhalf, r); }
	double m[N];
	for (int j = 0;j < N;j++) {
		uint64_t a[L];
		for (int i = 0;i < L;i++)
			a[i] = pt[i][j];
		BigInt A;
		icrt<L>(q, a, A);
		if (A >= Qhalf) {
			BigInt temp(A);
			Q.sub(temp, A);
			m[j] = -A.to_real() / Delta;
		}
		else
			m[j] = A.to_real() / Delta;
	}
	fft<N, logN>(m, zr, zi);
}

//
template<int N>
void keygen(int h, int s[N]) {
	for (int i = 0;i < N;i++) s[i] = 0;
	for (int i = 0;i < h;i++) {
		int j = rand() % N;
		while (s[j] != 0)
			j = rand() % N;
		s[j] = (rand() % 2 == 0) ? 1 : -1;
	}
}

//
template<int N, int L>
void enc(const uint64_t pt[L][N], const int s[N],
	const uint64_t q[L], uint64_t cthat[2][L][N]) {
	int e[N]; for (int i = 0;i < N;i++) e[i] = (rand() % 17) - 8;
	for (int j = 0;j < L;j++) {
		uint64_t shat[N];
		for (int i = 0;i < N;i++) {
			cthat[1][j][i] = rand_uint64() % q[j];
			shat[i] = (q[j] + s[i]) % q[j];
			cthat[0][j][i] = (q[j] + e[i] + pt[j][i]) % q[j];
		}
		ntt<N>(cthat[0][j], q[j]);
		ntt<N>(cthat[1][j], q[j]);
		// debug
		{
			intt<N>(cthat[1][j], q[j]);
			ntt<N>(cthat[1][j], q[j]);
		}
		ntt<N>(shat, q[j]);
		for (int i = 0;i < N;i++)
			cthat[0][j][i] = (cthat[0][j][i] + q[j] - mul_mod(shat[i], cthat[1][j][i], q[j])) % q[j];
	}
}

//
template<int N, int L>
void dec(const uint64_t cthat[2][L][N], const int s[N],
	const uint64_t q[L], uint64_t pt[L][N]) {
	for (int j = 0;j < L;j++) {
		uint64_t shat[N];
		for (int i = 0;i < N;i++)
			shat[i] = (q[j] + s[i]) % q[j];
		ntt<N>(shat, q[j]);

		for (int i = 0;i < N;i++)
			pt[j][i] = (cthat[0][j][i] + mul_mod(cthat[1][j][i], shat[i], q[j])) % q[j];
		intt<N>(pt[j], q[j]);
	}
}

//
template<int N, int L>
void RS(const uint64_t q[L], const uint64_t pt[L][N], uint64_t pt_[L - 1][N]) {
	uint64_t qinv[L - 1];
	for (int j = 0;j < L - 1;j++)
		qinv[j] = inv_mod(q[L - 1], q[j]);
	for (int i = 0;i < N;i++) {
		uint64_t r = pt[L - 1][i] > (q[L - 1] / 2) ? 1 : 0;
		for (int j = 0;j < L - 1;j++)
			pt_[j][i] = (mul_mod((pt[j][i] + q[j] - (pt[L - 1][i] % q[j])) % q[j], qinv[j], q[j]) + r) % q[j];
	}
}

// ADDED by KSH 1024
// RS pt Lfr -> Lto
template<int N, int Lfr, int Lto>
void RS(const uint64_t q[Lfr], const uint64_t pt_orig[Lfr][N], uint64_t pt_[Lto][N]) {
	assert(Lfr > Lto);
	uint64_t pt[Lfr][N];
	for (int i = 0; i < Lfr; ++i)
		for (int j = 0; j < N; ++j)
			pt[i][j] = pt_orig[i][j];
	for (int k = Lfr - 1;k >= Lto;--k) {
		uint64_t qinv[Lfr];
		for (int j = 0;j < k;++j)
			qinv[j] = inv_mod(q[k], q[j]);
		for (int i = 0;i < N;++i) {
			uint64_t r = (pt[k][i] > q[k] / 2) ? 1 : 0;
			for (int j = 0;j < k;++j) {
				pt[j][i] = mul_mod((pt[j][i] + q[j] - (pt[k][i] % q[j])) % q[j], qinv[j], q[j]);
				if (k == Lto)
					pt_[j][i] = (pt[j][i] + r) % q[j];
			}
		}
	}
}

// MOVED FROM NumberTheory.h 1024
template<int N, int L>
void RS_hat(const uint64_t q[L],
	const uint64_t pt_hat[L][N],
	uint64_t pt_hat_[L - 1][N]) {
	uint64_t temp[L][N];
	for (int i = 0;i < L;i++)
		for (int j = 0;j < N;j++)
			temp[i][j] = pt_hat[i][j];

	intt<N, L>(q, temp);
	RS<N, L>(q, temp, pt_hat_);
	ntt<N, L - 1>(q, pt_hat_);
}

// ADDED by KSH 1024
// RS_hat ct Lfr -> Lto
template<int N, int Lfr, int Lto>
void RS_hat(const uint64_t q[Lfr], const uint64_t ct_hat[2][Lfr][N], uint64_t ct_hat_[2][Lto][N]) {
	uint64_t ct[2][Lfr][N];
	for (int i = 0; i < 2; ++i) {
		for (int j = 0; j < Lfr; ++j) {
			for (int k = 0; k < N; ++k) {
				ct[i][j][k] = ct_hat[i][j][k];
			}
		}
	}
	intt<N, Lfr>(q, ct[0]);
	intt<N, Lfr>(q, ct[1]);
	RS<N, Lfr, Lto>(q, ct[0], ct_hat_[0]);
	RS<N, Lfr, Lto>(q, ct[1], ct_hat_[1]);
	ntt<N, Lto>(q, ct_hat_[0]);
	ntt<N, Lto>(q, ct_hat_[1]);
}

// 
template<int N>
void conv(const int A[N],
	const int B[N], int C[N]) {

	for (int i = 0; i < N; i++) {
		C[i] = 0;
		for (int j = 0; j <= i; j++) C[i] += A[j] * B[i - j];
		for (int j = i + 1; j < N; j++) C[i] -= A[j] * B[i - j + N];
	}
}

//
template<int N>
void conv(const double a[N],
	const double b[N], double c[N]) {
	for (int i = 0; i < N; i++) {
		c[i] = 0;
		for (int j = 0; j < i + 1; j++)
			c[i] += a[j] * b[i - j];
		for (int j = i + 1; j < N; j++)
			c[i] -= a[j] * b[i + N - j];
	}
}

//
template<int N, int L>
void keygen(const int s[N], const uint64_t q[L],
	uint64_t pk_hat[2][L][N]) {
	uint64_t zero[L][N];
	for (int i = 0; i < L; i++)
		for (int j = 0; j < N; j++)
			zero[i][j] = 0;
	enc<N, L>(zero, s, q, pk_hat);
}

//
template<int N, int L>
void enc(const uint64_t pt[L][N],
	const uint64_t pk_hat[2][L][N], const uint64_t q[L],
	uint64_t ct_hat[2][L][N]) {
	int e0[N], e1[N], v[N];
	for (int i = 0; i < N; i++) {
		e0[i] = (rand() % 17) - 8;
		e1[i] = (rand() % 17) - 8;
		int r = rand() % 4;
		v[i] = (r == 0) ? (-1) : ((r == 1) ? 1 : 0);
	}

	for (int j = 0; j < L; j++) {
		uint64_t v_hat[N], pt_j_hat[N], e0_hat[N], e1_hat[N];
		for (int i = 0; i < N; i++) {
			v_hat[i] = (q[j] + v[i]) % q[j];
			pt_j_hat[i] = (q[j] + pt[j][i]) % q[j];
			e0_hat[i] = (q[j] + e0[i]) % q[j];
			e1_hat[i] = (q[j] + e1[i]) % q[j];
		}

		ntt<N>(v_hat, q[j]);
		ntt<N>(pt_j_hat, q[j]);
		ntt<N>(e0_hat, q[j]);
		ntt<N>(e1_hat, q[j]);

		for (int i = 0; i < N; i++) {
			ct_hat[0][j][i] = (mul_mod(v_hat[i], pk_hat[0][j][i], q[j]) + pt_j_hat[i] + e0_hat[i]) % q[j];
			ct_hat[1][j][i] = (mul_mod(v_hat[i], pk_hat[1][j][i], q[j]) + e1_hat[i]) % q[j];
		}
	}
}

*/