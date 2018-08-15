#pragma once
#include <iostream>
#include <vector>
#include<complex>
class FT
{
private:
	friend std::vector<std::vector<double>> operator*(std::vector<std::vector<double>>,std::vector<std::vector<double>>);
public:
	FT();
	void DiscreteFourierTransform(int** InputImage, int** OutputImage, double ** FreqReal, double ** FreqImag, int h, int w);
	void DFT(double** pFreqReal, double** pFreqImag, int** InputImage, int h, int w, int u, int v);

	void InverseDiscreteFourierTransform(int** InputImage, int** OutputImage, double ** FreqReal, double ** FreqImag, int h, int w);
	void InverseDFT(double ** InverseReal, double ** InverseImag, double ** pFreqReal, double ** pFreqImag, int h, int w, int x, int y);

	void FastFourierTransform(int** InputImage, int** OutputImage, double ** FreqReal, double ** FreqImag, int h, int w);
	void FFT(std::vector<std::vector<double>>&, std::vector<std::vector<double>>&, 
			 std::vector<std::vector<double>>&, int, int, int);

	void InverseFastFourierTransform(int** InputImage, int** OutputImage, double ** FreqReal, double ** FreqImag, int h, int w);
	void InverseFFT(std::vector<std::vector<double>>& ,std::vector<std::vector<double>>&,
					std::vector<std::vector<double>>&,int h);
	
	void lowpass(int **,int **,int,int);
	void LowpassFilter(int**,int**,double** Real, double** Img,int);
	void highpass(int **, int **, int, int);
	void HighpassFilter(int**, int**, double**, double**, int);

private:

};


