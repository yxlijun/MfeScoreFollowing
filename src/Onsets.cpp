#include "Onsets.h"
#include <iostream>
using namespace std;

Onset::Onset(){
	pre_max = 6;
	pre_avg = 20;
	combine = 7;
	threshold = 1.8;
	last_onset = 0;
	fps = spectralOdf.getfps();
}
Onset::~Onset(){

}


bool Onset::detect(vector<double> &spectrum){
	cursfresult = spectralOdf.CalSFOnline(spectrum);
	SFResult.push_back(cursfresult);
	int index = SFResult.size() - 1;
	int start_max = (index - pre_max) >= 0 ? index - pre_max : 0;
	int start_avg = (index - pre_avg) >= 0 ? index - pre_avg : 0;
	vector<double> maxfilter;
	vector<double> avgfilter;
	for (int j = start_max; j <= index; j++){
		maxfilter.push_back(SFResult[j]);
	}
	for (int j = start_max; j <= index; j++){
		avgfilter.push_back(SFResult[j]);
	}
	double mov_max = *max_element(maxfilter.begin(), maxfilter.end());
	double mov_avg = accumulate(avgfilter.begin(), avgfilter.end(), 0.0) / static_cast<double>(avgfilter.size());
	bool curmax = cursfresult == mov_max;
	bool curavg = cursfresult >= (mov_avg + threshold);
	bool currange = index > (last_onset + combine);
	if (curmax && curavg && currange){
		double onset = static_cast<double>(index) / fps;
		dectectionResult.push_back(onset);
		last_onset = index;
		return true;
	}
	return false;
}

void Onset::saveResult(string path){
	ofstream fout;
	fout.open(path);
	for (int i = 0; i < dectectionResult.size(); i++){
		fout << dectectionResult[i] << endl;
	}
	fout.close();
}

vector<double> Onset::getdecectionResult(){
	return dectectionResult;
}


vector<double> Onset::getSFresult(){
	return SFResult;
}