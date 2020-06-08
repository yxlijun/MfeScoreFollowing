#include "SpectralODF.h"
#include <iostream>
#include <fstream>
using namespace std;

SpectralODF::SpectralODF(){
	SPECTRUMLENGTH = 513;
	floors = 5;
	relaxation = 10;
	fps = static_cast<double>(44100) / 512;
	mul = 10;
	Filter();
}

SpectralODF::~SpectralODF(){

}

double SpectralODF::getfps(){
	return fps;
}


void SpectralODF::SpectralFlux(){
	if (SFResultOnline.empty()){
		SFResultOnline.push_back(0.0);
		return;
	}
	else{
		int SpecOnLineSize = spectrumOnline.size() - 1;
		vector<double> diff;
		for (int i = 0; i < spectrumOnline[SpecOnLineSize].size(); i++){
			diff.push_back(spectrumOnline[SpecOnLineSize][i] - spectrumOnline[SpecOnLineSize - 1][i]);
		}
		for (int j = 0; j < spectrumOnline[SpecOnLineSize].size(); j++){
			if (diff[j] < 0)
				diff[j] = 0;
		}
		double sum = accumulate(diff.begin(), diff.end(), 0.0);
		SFResultOnline.push_back(sum);
	}
}


void SpectralODF::spectrumLog(){
	for (int i = 0; i < curspectrum.size(); i++){
		curspectrum[i] = log10(mul*curspectrum[i] + 1);
	}
}

vector<double> SpectralODF::getSFresultOnline(){
	return SFResultOnline;
}

double SpectralODF::CalSFOnline(vector<double> &spec){
	curspectrum.resize(filterbank[0].size());
	for (int j = 0; j < filterbank[j].size(); j++){
		for (int k = 0; k <SPECTRUMLENGTH; k++){
			curspectrum[j] += spec[k] * filterbank[k][j];
		}
	}
	spectrumLog();
	spectrumOnline.push_back(curspectrum);
	SpectralFlux();
	return SFResultOnline.back();
}



void SpectralODF::frequencies(){
	double fminx = 27.5;
	double fmaxx = 16000;
	int bands = 12;
	double a = 440;
	double factor = pow(2,1.0/bands);
	double freq = a;
	frequency.push_back(freq);
	while (freq<=fmaxx){
		freq *= factor;
		frequency.push_back(freq);
	}
	freq = a;
	while (freq>=fminx){
		freq /= factor;
		frequency.push_back(freq);
	}
	sort(frequency.begin(), frequency.end());
}

vector<double> SpectralODF::triang(int start, int mid, int stop){
	double height = 1.0;
	vector<double> triang_filter(stop - start);
	for (int i = 0; i < mid-start; i++){
		triang_filter[i] = i * (height / (mid - start));
	}
	double ave = height / (stop-mid);
	for (int i = 0; i < stop-mid; i++){
		triang_filter[i+mid-start] = height - i*ave;
	}
	return triang_filter;
}

void SpectralODF::Filter(double fs, double fminx, double fmaxx){
	if (fmaxx > fs / 2)
		fmaxx = fs / 2;
	frequencies();
	double factor = (fs/2.0) / SPECTRUMLENGTH;
	for (int i = 0; i < frequency.size(); i++){
		frequency[i] = round(frequency[i] / factor);
	}
	vector<double>::iterator pos = unique(frequency.begin(), frequency.end());
	frequency.erase(pos, frequency.end());
	vector<int> frequencyTemp;
	for (int i = 0; i < frequency.size(); i++){
		if (frequency[i] < SPECTRUMLENGTH)
			frequencyTemp.push_back(frequency[i]);
	}
	int bands = frequencyTemp.size() - 2;
	filterbank.resize(SPECTRUMLENGTH);
	for (int i = 0; i < SPECTRUMLENGTH; i++){
		filterbank[i].resize(bands);
	}
	for (int i = 0; i < bands; i++){
		int start = frequencyTemp[i];
		int mid = frequencyTemp[i + 1];
		int stop = frequencyTemp[i + 2];
		vector<double> triang_filter = triang(start, mid, stop);
		for (int j = start; j <stop; j++){
			filterbank[j][i] = triang_filter[j - start];
		}
	}
}