#include "FT.h"
#include <vector>
#include <complex>
FT::FT()
{
}

void FT::DiscreteFourierTransform(int** InputImage, int** OutputImage, double ** FreqReal, double ** FreqImag, int h, int w)
{
	int M = h;
	int N = w;
	
	double** pFreq = new double*[M];
	for (int newcnt = 0; newcnt<M; newcnt++)
	{
		pFreq[newcnt] = new double[N]; // 傅立葉頻率陣列
	}
	for (int forzero_i = 0; forzero_i<M; forzero_i++)
	{
		for (int forzero_j = 0; forzero_j<N; forzero_j++)
		{
			pFreq[forzero_i][forzero_j] = 0.0f;
		}
	}
	//-------------------------------------------
	for (int i = 0; i < M; i++)
	{
		for (int j = 0; j < N; j++)
		{
			DFT(FreqReal, FreqImag, InputImage,M, N, j, i);
		}
	}
	for (int i = 0; i < M; i++)
	{
		for (int j = 0; j < N; j++)
		{
			// 將計算好的傅立葉實數與虛數部分作結合 
			pFreq[i][j] = sqrt(pow(FreqReal[i][j], (double) 2.0) + pow(FreqImag[i][j], (double) 2.0));
			// 結合後之頻率域丟入影像陣列中顯示 
			OutputImage[i][j] = pFreq[i][j];
		}
	}
	//-------------------------------------------
	for (int delcnt = 0; delcnt < M; delcnt++)
	{
		delete[] pFreq[delcnt];
	}
	delete[] pFreq;
}

void FT::DFT(double ** pFreqReal, double ** pFreqImag, int ** InputImage, int h, int w, int u, int v)
{
	// M = N 必須是方陣
	int M = h;
	int N = w;

	for (int y = 0; y < M; y++)
	{
		for (int x = 0; x < N; x++)
		{
			// 可先計算Eular's equation e^{j*theta} = cos{theta}+j*sin{theta}			
			double angleDFT = (-1.0f * 2.0f * 3.14159 * (double)(u*x + v*y) / (double)M);
			double c = cos(angleDFT);
			double s = sin(angleDFT);

			// 利用Eular's equation計算傅立葉之實虛數部分
			pFreqReal[u][v] += (double)InputImage[y][x] * c;
			pFreqImag[u][v] -= (double)InputImage[y][x] * s;
		}
	}

	pFreqReal[u][v] = pFreqReal[u][v] / (double)(M);
	pFreqImag[u][v] = pFreqImag[u][v] / (double)(M);
}

void FT::InverseDiscreteFourierTransform(int ** InputImage, int ** OutputImage, double ** FreqReal, double ** FreqImag, int h, int w)
{
	int M = h;
	int N = w;

	double** InverseReal = new double*[M];
	double** InverseImag = new double*[M];
	double** pFreq = new double*[M];

	for (int i = 0; i<M; i++)
	{
		InverseReal[i] = new double[N];
		InverseImag[i] = new double[N];
		pFreq[i] = new double[N]; // 傅立葉頻率陣列
	}

	for (int i = 0; i<M; i++)
	{
		for (int j = 0; j<N; j++)
		{
			InverseReal[i][j] = 0.0f;
			InverseImag[i][j] = 0.0f;
			pFreq[i][j] = 0.0f;
		}
	}
	//-------------------------------------------
	for (int i = 0; i < M; i++)
	{
		for (int j = 0; j < N; j++)
		{
			InverseDFT(InverseReal, InverseImag,FreqReal, FreqImag, M, N, j, i);
		}
	}
	for (int i = 0; i < M; i++)
	{
		for (int j = 0; j < N; j++)
		{
			// 將計算好的傅立葉實數與虛數部分作結合 
			pFreq[i][j] = sqrt(pow(InverseReal[i][j], (double) 2.0) + pow(InverseImag[i][j], (double) 2.0));
			// 結合後之頻率域丟入影像陣列中顯示 
			OutputImage[i][j] = pFreq[i][j];
			//存下反傅立葉實數與虛數部分
			FreqReal[i][j] = InverseReal[i][j];
			FreqImag[i][j] = InverseImag[i][j];

		}
	}
	//-------------------------------------------
	for (int i = 0; i < M; i++)
	{
		delete[] pFreq[i];
		delete[] InverseReal[i];
		delete[] InverseImag[i];

	}
	delete[] pFreq;
	delete[] InverseReal;
	delete[] InverseImag;

}

void FT::InverseDFT(double ** InverseReal, double ** InverseImag, double ** pFreqReal, double ** pFreqImag, int h, int w, int x, int y)
{
	// M = N 必須是方陣
	int M = h;
	int N = w;

	for (int v = 0; v < M; v++)
	{
		for (int u = 0; u < N; u++)
		{
			// 可先計算Eular's equation e^{j*theta} = cos{theta}+j*sin{theta}			
			double angleIDFT = (2.0f * 3.14159 * (double)(u*x + v*y) / (double)M);
			double c = cos(angleIDFT);
			double s = sin(angleIDFT);

			// 利用Eular's equation計算傅立葉之實虛數部分
			InverseReal[x][y] += (pFreqReal[v][u] * c - pFreqImag[v][u] * s);
			InverseImag[x][y] += (pFreqReal[v][u] * s + pFreqImag[v][u] * c);
		}
	}
	InverseReal[x][y] = InverseReal[x][y] / (float)M;
	InverseImag[x][y] = InverseImag[x][y] / (float)M;
}

void FT::FastFourierTransform(int ** InputImage, int ** OutputImage, double ** FreqReal, double ** FreqImag, int h, int w)
{	
	int M = h;
	int N = w;
	std::vector<std::vector<double>> x;
	std::vector<std::vector<double>> X;
	std::vector<std::vector<double>> temp;
	double** pFreqReal = new double*[M];
	double** pFreqImag = new double*[M];
	double** pFreq = new double*[M];

	for (int i = 0; i< M; i++) {
		pFreqReal[i] = new double[N]; // 傅立葉實數部分
		pFreqImag[i] = new double[N]; // 傅立葉虛數部分
		pFreq[i] = new double[N]; // 傅立葉頻率陣列
	}
	for (int init_i = 0; init_i < M; init_i++) {
		for (int init_j = 0; init_j < N; init_j++) {
			pFreqReal[init_i][init_j] = 0.0f;
			pFreqImag[init_i][init_j] = 0.0f;
			pFreq[init_i][init_j] = 0.0f;
		}
	}

	for (int init = 0; init < M; init++) {
		std::vector<double> t(2,0);
		temp.push_back(t);
		x.push_back(t);
		X.push_back(t);
	}
	//-------------------------------------------
	for (int i = 0; i < M; i++) {
		for (int j = 0; j < N; j++) 
			x[j][0] = InputImage[j][i];
		FFT(X, x, temp, N, 0, 1);
		for (int k = 0; k<N; k++) {
			pFreqReal[k][i] += (X[k][0] / (double)M);
			pFreqImag[k][i] += (X[k][1] / (double)M);
			pFreq[k][i] = std::sqrt(std::pow(pFreqReal[k][i], (double) 2.0) + 
							        std::pow(pFreqImag[k][i], (double) 2.0));
			// 結合後之頻率域丟入影像陣列中顯示 
			OutputImage[k][i] = pFreq[k][i];
		}
	}
	temp.clear();
	x.clear();
	X.clear();

	for (int init = 0; init < M; init++) {
		std::vector<double> t(2, 0);
		temp.push_back(t);
		x.push_back(t);
		X.push_back(t);
	}

	for (int j = 0; j < N; j++) {
		for (int i = 0; i < M; i++) {
			x[i][0] = pFreqReal[j][i];
			x[i][1] = pFreqImag[j][i];
		}
		FFT(X, x, temp, N, 0, 1);
		for (int k = 0; k < M; k++) {
			pFreqReal[j][k] = (X[k][0]);
			pFreqImag[j][k] = (X[k][1]);
			pFreq[j][k] = std::sqrt(std::pow(pFreqReal[j][k], 2) + 
						  std::pow(pFreqImag[j][k],2));
			// 結合後之頻率域丟入影像陣列中顯示 
			OutputImage[j][k] = pFreq[j][k];
		}
	}

	//-------------------------------------------
	for (int i = 0; i < M; i++) {
		delete[] pFreqReal[i];
		delete[] pFreqImag[i];
		delete[] pFreq[i];
	}
	delete[] pFreqReal;
	delete[] pFreqImag;
	delete[] pFreq;
}

void FT::FFT(std::vector<std::vector<double>> &X, std::vector<std::vector<double>> &x,
			 std::vector<std::vector<double>> &temp,int N, int odd, int d){
	double pi = 3.14159;
	//deal with the value if the size is 2*2
	if (N == 2){
		int i=0, j=0;
		i = odd;
		j = i + d;
		X[i][0] = x[i][0] + x[j][0];
		X[i][1] = x[i][1] + x[j][1];
		X[j][0] = x[i][0] - x[j][0];
		X[j][1] = x[i][1] - x[j][1];
		return;
	}

	int i00=0, i01=0, i10=0, i11=0;
	//deal with even index
	FFT(temp, x, X, N / 2, odd, 2 * d);
	//deal with odd index
	FFT(temp, x, X, N / 2, odd + d, 2 * d);
	//compute omega into the xs
	//because odd index and even can be done at the same time
	for (int k = 0; k < N / 2; k++) {
		i00 = odd + k*d;
		i01 = i00 + (N / 2)*d;
		i10 = odd + 2 * k*d;
		i11 = i10 + d;
		double angle = -2.0 * pi * k / (double)N;
		double c = std::cos(angle);
		double s = std::sin(angle);
		double t1 = 0, t2 = 0;

		t1 = c * temp[i11][0] + s * temp[i11][1];
		t2 = c * temp[i11][1] - s * temp[i11][0];

		X[i00][0] = temp[i10][0] + t1;
		X[i00][1] = temp[i10][1] + t2;
		X[i01][0] = temp[i10][0] - t1;
		X[i01][1] = temp[i10][1] - t2;
	}
}

void FT::InverseFastFourierTransform(int ** InputImage, int ** OutputImage, double ** FreqReal, double ** FreqImag, int h, int w)
{
}

void FT::InverseFFT(double ** InverseReal, double ** InverseImag, double ** pFreqReal, double ** pFreqImag, int h, int w, int x, int y)
{
}
void FT::LowpassFilter(double** Real, double** Img, double** filter)
{
}

void FT::HighpassFilter(double** Real, double** Img, double** filter)
{

}
void FT::highpass(int **input,int **output,int h,int w) {
	int M = h, N = w;
	int mask[3][3] = { (-1,-1,-1),(-1,8,-1),(-1,-1,-1)};
	for (int i = 0; i < h; i++) {
		for (int j = 0; j < w; j++) 
			output[i][j] = input[i][j];
	}
	for (int i = 1; i < M-1;i++) {
		for (int j = 1; j < N - 1;j++) {
			int sum = 0.0;
			sum = -1 * input[i - 1][j] + -1 * input[i][j - 1] +
				   4 * input[i][j] + -1 * input[i][j + 1] +
				  -1 * input[i + 1][j];
			output[i][j] = sum;
		}
	}
}

void FT::lowpass(int**inputimage,int **outputimage,int h,int w) {
	int M = h, N = w;
	int **mask = new int*[3];

	for (int newi = 0; newi < 3; newi++) {
		mask[newi] = new int[3];
		for (int init = 0; init < 3; init++) mask[newi][init] = 0.0f;
	}

	for (int i = 0; i < h; i++) {
		for (int j = 0; j < w; j++) {
			outputimage[i][j] = inputimage[i][j];
		}
	}

	for (int i = 0; i < M - 2; i++) {
		for (int j = 1; j < N - 1; j++) {
			mask[0][0] = inputimage[i][j - 1];
			mask[0][1] = inputimage[i][j];
			mask[0][2] = inputimage[i][j + 1];
			mask[1][0] = inputimage[i + 1][j - 1];
			mask[1][1] = inputimage[i + 1][j];
			mask[1][2] = inputimage[i + 1][j + 1];
			mask[2][0] = inputimage[i + 2][j - 1];
			mask[2][1] = inputimage[i + 2][j];
			mask[2][2] = inputimage[i + 2][j + 1];
			int sum = 0;
			for (int r = 0; r < 3; r++) {
				for (int s = 0; s < 3; s++) sum += mask[r][s];
			}
			outputimage[i][j] = sum / 9;
		}
	}

	for (int i = 0; i < 3; i++) delete[] mask[i];
	delete mask;
}
