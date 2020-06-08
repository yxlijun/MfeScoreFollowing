#ifndef SPECTRALODF_H
#define SPECTRALODF_H

#include <string>
#include "kiss_fft.h"
#include <vector>
//#include <resample.h>
#include <algorithm>
#include <cmath>
#include <numeric>

using namespace std;
class SpectralODF
{
public:
	SpectralODF();
	~SpectralODF();


	void SpectralFlux();

	void spectrumLog();

	double getfps();

	vector<double> getSFresultOnline();

	double CalSFOnline(vector<double> &spec);

	void frequencies();

	vector<double> triang(int start, int mid, int stop);

	void Filter(double fs = 44100,double fminx = 27.5,double fmaxx=16000);
private:
	int SPECTRUMLENGTH;              // Æ×³¤¶È£¨fftSize/2+1£©
	int floors;
	int relaxation;

	double fps;
	int mul;
	vector<double> SFResultOnline;
	vector<vector<double>> spectrumOnline;
	vector<double> curspectrum;
	vector<double> frequency;
	vector<vector<double>> filterbank;
};


#endif