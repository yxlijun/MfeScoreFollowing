#pragma once

#include "SpectralODF.h"
#include <vector>
#include <fstream>
#include "Timer.h"

class  Onset
{
public:
	Onset();
	~ Onset();

	bool detect(vector<double> &spectrum);

	void saveResult(string path);

	vector<double> getdecectionResult();

	vector<double> getSFresult();

private:
	SpectralODF spectralOdf;
	string audiopath;
	double fps;
	vector<double> SFResult;
	vector<double> dectectionResult;
	bool online;
	double cursfresult;
	int last_onset;
	bool LessCombine;

	double threshold;     
	double combine;       // w5
	double pre_avg;      // w3
	double pre_max;		// w1  
	double delay;          
};