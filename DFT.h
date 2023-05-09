#pragma once
#include <cassert>
#include <complex>

#define PI 3.1415926535897932384626433

/*
//-------------------------------------------------------------------------------
// dft and idft
//-------------------------------------------------------------------------------
template< int N >
void dft( const double m[N], double zr[N/2],
	                         double zi[N/2]){
	for(int i=0,poweri=1; i<N/2; i++,poweri=(5*poweri)%(2*N)){
		zr[i]=0; zi[i]=0;
		for(int j=0,powerij=0; j<N/2; j++,powerij=(powerij+poweri)%(2*N)){
			double c=cos(PI/N*powerij);
			double s=sin(PI/N*powerij);
			zr[i]+=c*m[j]-s*m[j+N/2];
			zi[i]+=c*m[j+N/2]+s*m[j];
		}
	}
}

template< int N >
void idft( const double zr[N/2],
	       const double zi[N/2], double m[N]){
	for(int i=0;i<N/2;i++){
		m[i]=0; m[i+N/2]=0;
		for(int j=0,powerij=i; j<N/2; j++,powerij=(powerij*5)%(2*N)){
			double c=cos(PI/N*powerij);
			double s=sin(PI/N*powerij);
			m[i    ]+=c*zr[j]+s*zi[j];
			m[i+N/2]+=-s*zr[j]+c*zi[j];
		}
	}
}
*/
//-------------------------------------------------------------------------------
// bitReverse
//-------------------------------------------------------------------------------
template< int N >
void bitReverse( const double a[N], double b[N] ){
	int logN=0;
	while(N>(1<<logN))
		logN++;
	for(int i=0;i<N;i++){
		int j=0;
		for(int k=0;k<logN;k++)
			j+=((i>>k)&1)<<(logN-1-k);
		b[j]=a[i];
	}
}

//-------------------------------------------------------------------------------
// SparseComplexMatrix
//-------------------------------------------------------------------------------
template<int N, int NumDiags>
struct SparseComplexMatrix{
	int shift[NumDiags];
	double diagr[NumDiags][N];
	double diagi[NumDiags][N];
	SparseComplexMatrix(){
		for(int k=0;k<NumDiags;k++){shift[k]=0;
		for(int i=0;i<N;i++) diagr[k][i]=diagi[k][i]=0;
		}
	}
	void applyA( const double xr[N],
				 const double xi[N], double Axr[N],
									 double Axi[N]) const{
		for(int i=0;i<N;i++){
			Axr[i]=Axi[i]=0;
			for(int k=0;k<NumDiags;k++){
				int j=(i+shift[k])%N;
				double aijr = diagr[k][i];
				double aiji = diagi[k][i];
				Axr[i]+=aijr*xr[j]-aiji*xi[j];
				Axi[i]+=aijr*xi[j]+aiji*xr[j];
			}
		}
	}
	void transpose(){
		for(int k=0;k<NumDiags;k++){
			int s=shift[k];
			shift[k]=(N-s)%N;
			double tempr[N], tempi[N];
			for(int i=0;i<N;i++){tempr[i]=diagr[k][i]; tempi[i]=diagi[k][i];}
			for(int i=0;i<N;i++){diagr[k][i]=tempr[(i+N-s)%N]; 
			                     diagi[k][i]=tempi[(i+N-s)%N];}
		}
	}
	void conjugate(){
		for(int k=0;k<NumDiags;k++)
			for(int i=0;i<N;i++) diagi[k][i]=-diagi[k][i];
	}
};

//-------------------------------------------------------------------------------
// splitU0R into FFT Sparse matrices
//-------------------------------------------------------------------------------
template< int logN >
void splitU0R( SparseComplexMatrix<(1<<(logN-1)),3> E[logN-1] ){
	int N = 1<<logN;
	for(int k=2,k_=0;k<=N/2;k*=2,k_++){
		int M=N/2/k; E[k_].shift[0]=0;
					 E[k_].shift[1]=M;
					 E[k_].shift[2]=N/2-M;
		double Wr[1<<(logN-2)];
		double Wi[1<<(logN-2)];
		for(int i=0,fiveik=k/2;i<M;i++,fiveik=(fiveik*5)%(2*N)){
			Wr[i]=cos(PI/N*fiveik);
			Wi[i]=sin(PI/N*fiveik);
		}
		for(int b=0;b<k/2;b++)
		for(int i=0;i<M;i++){
			E[k_].diagr[0][2*M*b  +i]=    1 ;
			E[k_].diagr[2][2*M*b+M+i]=    1 ;
			E[k_].diagr[1][2*M*b  +i]= Wr[i];
			E[k_].diagi[1][2*M*b  +i]= Wi[i];
			E[k_].diagr[0][2*M*b+M+i]=-Wr[i];
			E[k_].diagi[0][2*M*b+M+i]=-Wi[i];
		}
	}
}

//-------------------------------------------------------------------------------
// fft and ifft
//-------------------------------------------------------------------------------
template< int N, int logN >
void fft( const double m[N], double zr[N/2],
							 double zi[N/2]){
	SparseComplexMatrix<N/2,3> E[logN-1]; 
	splitU0R<logN>(E);
	bitReverse<N/2>(m    ,zr);
	bitReverse<N/2>(m+N/2,zi);
	double tempr[N/2], tempi[N/2];
	for(int k=logN-2;k>=0;k--){
		for(int i=0;i<N/2;i++){tempr[i]=zr[i]; tempi[i]=zi[i];}
		E[k].applyA(tempr,tempi,zr,zi);
	}
}

template< int N, int logN >
void ifft( const double zr[N/2],
		   const double zi[N/2], double m[N]){
	SparseComplexMatrix<(N/2),3> E[logN-1]; 
	splitU0R<logN>(E);
	double tempr[N/2], tempi[N/2];
	for(int i=0;i<N/2;i++){tempr[i]=zr[i]; tempi[i]=zi[i];}
	
	for(int k=0;k<logN-1;k++){
		E[k].transpose();
		E[k].conjugate();
		E[k].applyA(tempr,tempi,m,m+N/2);
		for(int i=0;i<N/2;i++){tempr[i]=m[i]; tempi[i]=m[N/2+i];}
	}
	bitReverse<N/2>(tempr,m); bitReverse<N/2>(tempi,m+N/2);
	for(int i=0;i<N;i++) m[i]*=2./N;
}