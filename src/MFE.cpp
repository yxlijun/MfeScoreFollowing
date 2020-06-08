#include "MFE.h"
#include <iostream>
MFE::MFE() :originNTemplate(89)
{
	noteSize = 88;
	nTemplate = 89;
	maxPitches = 4;
	threshCoeff = 0.2;
	minThresh = 240;
	framecol = 4096;
	niter = 15;
	firstCheck = false;
}
MFE::~MFE()
{
}
void MFE::paramInit(){
	noteSize = 88;
	nTemplate = 89;
	maxPitches = 4;
	threshCoeff = 0.2;
	minThresh = 240;
	framecol = 4096;
	niter = 15;
}
void MFE::setresample(bool resample){
	isresample = resample;
	if (isresample){
		FFTSIZE = 1024;
		SPECTRUMLENGTH = 513;
		HOPSIZE = 110;
		FRAMESIZE = 1024;
	}
	else{
		FFTSIZE = 4096;
		SPECTRUMLENGTH = 2049;
		HOPSIZE = 441;
		FRAMESIZE = 4096;
	}
}

int MFE::getrow(){
	return framerow;
}
void  MFE::ReadFrameFile(const string &framePath, vector<vector<double>>& xFrame){
	ifstream HStream(framePath.c_str());
	framerow = 0;
	string tmp;
	while (getline(HStream, tmp, '\n')){
		framerow++;
	}
	HStream.close();

	xFrame.resize(framerow);
	for (int i = 0; i < framerow; i++){
		xFrame[i].resize(framecol);
	}
	FILE *xFrameFile;
	if ((xFrameFile = fopen(framePath.c_str(), "rb")) == NULL){
		printf("cannot open!\n");
		exit(1);
	}
	for (int i = 0; i < framerow; ++i){
		for (int j = 0; j < framecol; ++j){
			fscanf(xFrameFile, "%lf", &xFrame[i][j]);
		}
	}
	fclose(xFrameFile);

}


void MFE::ReadTemplateFile(const string &templatePath, vector<double>& templateData, vector<int> &notesInScore){
	templateData.resize(SPECTRUMLENGTH*nTemplate);
	FILE* templateFile;
	if ((templateFile = fopen(templatePath.c_str(), "rb")) == NULL){
		printf("cannot open template.dat!\n");
		exit(1);
	}
	for (int i = 0; i < SPECTRUMLENGTH*nTemplate; ++i){
		fscanf(templateFile, "%lf", &templateData[i]);
	}

	int newTemplateColumn;

	double *newTemplate = chooseTemplate(&templateData[0], &notesInScore[0], notesInScore.size(),&newTemplateColumn);
	vector<double>().swap(templateData);
	templateData.resize(SPECTRUMLENGTH*nTemplate);
	for (int i = 0; i <SPECTRUMLENGTH*nTemplate; i++){
		templateData[i] = newTemplate[i];
	}
	fclose(templateFile);
}

void MFE::GenerateH(vector<vector<double>> &H, vector<double> &error, vector<vector<double>> &ratios, vector<vector<double>> &xFrame, vector<double> &templateData, const string &path, vector<int> &notesInScore,vector<vector<double>> &Specurm,vector<vector<double>> &outframe){
	//getframeOnset();
	H.resize(framerow);
	vector<vector<double>> xH;
	xH.resize(framerow);
	Specurm.resize(framerow);
	outframe.resize(framerow);
	for (int i = 0; i < framerow; i++){
		H[i].resize(88);
		xH[i].resize(noteSize);
		Specurm[i].resize(SPECTRUMLENGTH);
		outframe[i].resize(FRAMESIZE);
	}
	int H0Flag = 0;
	vector<double> H0(nTemplate,1);
	//H0.resize(nTemplate);
	notesScore = notesInScore;
	for (int i = 0; i < framerow; i++){
		vector<double> tempratio;
		double *ratio = (double*)calloc(3, sizeof(double));
		double *errs = (double*)calloc(niter, sizeof(double));
		double *specrum = (double*)calloc(SPECTRUMLENGTH, sizeof(double));
		double *Outframe = (double*)calloc(FRAMESIZE, sizeof(double));

		memset(errs, 0, sizeof(double)*niter);
		multiF0Estimation(&xFrame[i][0], &templateData[0], H0Flag, &H0[0], &xH[i][0],errs,path,ratio,specrum,Outframe,onsetdec);
		for (int j = 0; j <notesInScore.size(); j++){
			H[i][notesInScore[j] - 1] = xH[i][j];
		}
		for (int j = 0; j < SPECTRUMLENGTH; j++){
			Specurm[i][j] = specrum[j];
		}
		for (int j = 0; j < FRAMESIZE; j++){
			outframe[i][j] = Outframe[j];
		}
		error.push_back(errs[niter-1]);
		for (int j = 0; j < 3; j++){
			tempratio.push_back(ratio[j]);
		}
		ratios.push_back(tempratio);
		free(ratio); ratio = NULL;
		free(errs); errs = NULL;
		free(specrum); specrum = NULL;
		free(Outframe); Outframe = NULL;

	}
	onsetdec.saveResult("../data/onset.dat");
}


int* MFE::MFEinterface(short* originWav, int wChannels, double* templateData, double* xH,double * error, const string &path,double *ratio){
	int i;
	// 得到的原始录音数据，先归一化
	double* normalizedData = (double*)malloc(FRAMESIZE * sizeof(double));
	for (i = 0; i < FRAMESIZE; ++i) {
		normalizedData[i] = (double)2*originWav[i] / 32767.0000;
	}

	// 读入频谱模板数据
	// double* templateData = (double *)malloc(sizeof(double)*SPECTRUMLENGTH*nTemplate);
	// getTemplateData(templateData);

	// 调用多音调检测算法

	int H0Flag = 0;
	double* H0 = (double*)calloc(nTemplate, sizeof(double));
	int* pianoRoll = (int*)calloc(noteSize, sizeof(int));
	//multiF0Estimation(normalizedData, templateData, H0Flag, H0, pianoRoll,xH);
	//multiF0Estimation(normalizedData, templateData, H0Flag, H0, xH,path);
	// 释放分配的内存
	free(normalizedData); normalizedData = NULL;
	free(H0); H0 = NULL;

	return pianoRoll;
}

void MFE::multiF0Estimation(double* xFrame, double* templateData, int &H0Flag, double* H0, double* xH, double *error, const string &path,double *ratio,double *specrum,double*OutFrame,Onset &onsetdec){
	vector<double> outframe;
	vector<double> xframe(4096);
	for (int i = 0; i < 4096; i++){
		xframe[i] = xFrame[i];
	}
	if (isresample)
		outframe = resample(xframe, 1024, 4096,path);
	else
		outframe = xframe;
	int i, j;
	// 先判断xFrame是否全为零，全为零就将pianoRoll置0返回
	double xFrameSum = 0.0;
	//for (i = 0; i < FRAMESIZE; ++i) xFrameSum += fabs(xFrame[i]);     
	for (i = 0; i < FRAMESIZE; ++i) xFrameSum += fabs(outframe[i]);

	if (xFrameSum < 1e-8) {
		for (i = 0; i < noteSize; ++i) {
			xH[i] = 0;
		}
		return;
	}

	// 计算信号频谱
	double* spectrum = (double*)malloc(sizeof(double)*SPECTRUMLENGTH);
	spectrumOfFrame(&outframe[0], spectrum);
	/*for (int i = 0; i < SPECTRUMLENGTH; i++){
		spectrum[i] = spectrum[i] * (1 + 9.0*i / 512.0);
	}*/

	vector<double> spec;
	for (int z = 0; z < SPECTRUMLENGTH; z++){
		specrum[z] = spectrum[z];
		spec.push_back(specrum[z]);
	}
	SpectrumOffLine.push_back(spec);
	for (int z = 0; z < FRAMESIZE; z++){
		OutFrame[z] = outframe[z];
	}
	// 判断spectrum 中是否含零，包含零就将pianoRoll置0返回
	for (i = 0; i < SPECTRUMLENGTH; ++i) {
		if ((spectrum[i] - 0) < 1e-8) {
			for (j = 0; j < noteSize; ++j) {
				xH[j] = 0;
			}
			// 先释放分配的内存，再返回，防止内存泄漏
			free(spectrum); spectrum = NULL;
			return;
		}
	}

	//频谱因式分解
	double* H = (double*)calloc(nTemplate, sizeof(double)); // 定义H并初始化为0
	Timer timer;
	timer.Tic();
	bool findonset = onsetdec.detect(spec);
	
	nmf_beta(spectrum, SPECTRUMLENGTH, nTemplate, templateData, H, H0Flag, H0, error,beta, niter,findonset);
	timer.Toc();
	//cout << "cost time:" << timer.Elasped() << endl;
	//nmf_beta_cblas(spectrum, SPECTRUMLENGTH, nTemplate, templateData, H, H0Flag, H0, beta, niter);

	double TopSpectrum = 0.0;
	double BottomSpectrum = 0.0;
	for (int i = 0; i < 300; i++){
		BottomSpectrum += spectrum[i];
	}
	for (int i = 300; i < SPECTRUMLENGTH; i++){
		TopSpectrum += spectrum[i];
	}
	double TopSpectrumMean = TopSpectrum / (double)(213);
	//cout << TopSpectrumMean << endl;
	int excessMeanCount = 0;
	for (int i = 300; i < SPECTRUMLENGTH; i++){
		if (spectrum[i] > TopSpectrumMean)
			excessMeanCount++;
	}
	if (BottomSpectrum!=0){
		ratio[0] = TopSpectrum / BottomSpectrum;
	}
	ratio[1] = TopSpectrum;
	ratio[2] = excessMeanCount;
	free(spectrum); spectrum = NULL;

	//有多套模板时，将对应于同一音调的pitch activity相加
	for (i = 1; i < nTemplate / noteSize; ++i) {
		for (j = 0; j < noteSize; ++j) {
			H[j] += H[j + i * noteSize];
		}
	}

	for (i = 0; i < noteSize; ++i) xH[i] = H[i];

	free(H); H = NULL;
}


void MFE::spectrumOfFrame(double* xFrame, double* spectrum){
	int i;
	double* win = (double*)calloc(FRAMESIZE, sizeof(double));
	hammingWindow(win, FRAMESIZE);
	//加窗
	double* xFft = (double*)calloc(FFTSIZE, sizeof(double));
	for (i = 0; i < FRAMESIZE; ++i)
	{
		xFft[i] = xFrame[i] * win[i];
	}
	if (FFTSIZE > FRAMESIZE) //超过的部分补0
	{
		for (i = FRAMESIZE; i < FFTSIZE; ++i)
		{
			xFft[i] = 0.0;
		}
	}
	free(win); win = NULL;

	//kiss-fft
	//input
	kiss_fft_cpx* kiss_xFft = (kiss_fft_cpx*)malloc(sizeof(kiss_fft_cpx)*FFTSIZE);
	for (i = 0; i < FFTSIZE; ++i)
	{
		kiss_xFft[i].r = xFft[i];
		kiss_xFft[i].i = 0.0;
	}
	free(xFft); xFft = NULL;

	//output
	kiss_fft_cpx* fftFrame = (kiss_fft_cpx*)malloc(sizeof(kiss_fft_cpx)*FFTSIZE);
	kiss_fft_cfg cfg = kiss_fft_alloc(FFTSIZE, 0, 0, 0);

	kiss_fft(cfg, kiss_xFft, fftFrame);

	//计算频谱振幅
	for (i = 0; i < SPECTRUMLENGTH; ++i)
	{
		spectrum[i] = sqrt(fftFrame[i].r * fftFrame[i].r + fftFrame[i].i * fftFrame[i].i); // 计算幅值
	}

	free(kiss_xFft); free(cfg);
	free(fftFrame); fftFrame = NULL;
}
void MFE::hammingWindow(double *w, int windowSize){
	int i;
	double pi = 3.14159265358979323846;
	for (i = 0; i < windowSize; ++i) {
		w[i] = 0.54 - 0.46 * cos(2 * pi * i / (windowSize - 1));
	}
}
void MFE::nmf_beta(double* data, int rows, int columns, double* W, double* H, int &H0Flag, double* H0, double *errs, double beta, int maxiter,bool findonset){
	int i, j, k;
	if (findonset) H0Flag = 0;
	else H0Flag = 1;
	// 初始化H
	if (H0Flag == 0/* || !firstCheck*/) { // 没有H0输入，则随机初始化
		for (i = 0; i < columns; ++i) {
			H[i] = 1;  // this is for test，the code below is for release
			//H[i] = rand() / ((double)RAND_MAX + 1);
		}
	}
	else if (H0Flag == 1/*&&firstCheck*/) { // 有H0输入，则将H初始化为H[0]
		for (i = 0; i < columns; ++i) {
			H[i] = H0[i];
		}
	}

	//    // 将模板标准化
	//    double* mrowsum = (double*)calloc(columns, sizeof(double));
	//    for (i = 0; i < columns; ++i) {
	//        for (j = 0; j < rows; ++j) {
	//            mrowsum[i] += w[j * columns + i]; // 每一列的总和
	//        }
	//    }
	//    for (i = 0; i < rows; ++i) {
	//        for (j = 0; j < columns; ++j) {
	//            w[i * columns + j] /= mrowsum[j];
	//        }
	//    }
	//    free(mrowsum); mrowsum = null;

	// 初始化R
	double* R = (double*)calloc(rows, sizeof(double));
	for (i = 0; i < rows; ++i) {
		for (k = 0; k < columns; ++k) {
			R[i] += W[i * columns + k] * H[k];  // 二维矩阵相乘
		}
	}

	double *temp1 = (double*)calloc(rows, sizeof(double));
	double *temp2 = (double*)calloc(rows, sizeof(double));
	double *temp3 = (double*)calloc(columns, sizeof(double));
	double *temp4 = (double*)calloc(columns, sizeof(double));

	//double *errs = (double*)calloc(maxiter, sizeof(double));
	//smemset(errs, 0, sizeof(double)*maxiter);
	double *tempV = (double*)calloc(rows, sizeof(double));
	double *tempR = (double*)calloc(rows, sizeof(double));
	double *tempV1 = (double*)calloc(rows, sizeof(double));
	double *tempR1 = (double*)calloc(rows, sizeof(double));
	double *tempVR = (double*)calloc(rows, sizeof(double));
	double *temp5 = (double*)calloc(rows, sizeof(double));

	vector<double> difErr;
	vector<double> Err;

	// 迭代
	double myeps = 1.0e-20; // 精度
	int it; // 迭代次数
	for (it = 0; it < maxiter; ++it) {
		//update R
		memset(R, 0, sizeof(double)*rows); // 将R清零
		for (i = 0; i < rows; ++i) {
			for (k = 0; k < columns; ++k) {
				R[i] += W[i * columns + k] * H[k];
			}
		}
		//update H
		// 初始化中间变量，内存清零
		memset(temp1, 0, sizeof(double)*rows);
		memset(temp2, 0, sizeof(double)*rows);
		memset(temp3, 0, sizeof(double)*columns);
		memset(temp4, 0, sizeof(double)*columns);
		for (i = 0; i < rows; ++i) {
			temp1[i] = pow(R[i], beta - 2) * data[i];
			temp2[i] = pow(R[i], beta - 1);
		}
		for (k = 0; k < rows; ++k) {
			for (i = 0; i < columns; ++i) {
				temp3[i] += W[k * columns + i] * temp1[k];
				temp4[i] += W[k * columns + i] * temp2[k];
			}
		}
		for (i = 0; i < columns; ++i) {
			if (temp4[i] < myeps) {
				temp4[i] = myeps;
			}
			H[i] = H[i] * (temp3[i] / temp4[i]);
		}
		//再一次更新update R
		memset(R, 0, sizeof(double)*rows); // 将R清零
		for (i = 0; i < rows; ++i) {
			for (k = 0; k < columns; ++k) {
				R[i] += W[i * columns + k] * H[k];
			}
		}
		double sum = 0.0;
		for (i = 0; i < rows; i++){
			tempR[i] = pow(R[i], beta);
			tempR[i] = tempR[i] * (beta - 1);
			tempV[i] = pow(data[i], beta);
			tempV1[i] = data[i] * beta;
			tempR1[i] = pow(R[i], (beta - 1));
			tempVR[i] = tempV1[i] * tempR1[i];
			temp5[i] = tempV[i] + tempR[i] - tempVR[i];
		}
		for (i = 0; i < rows; i++){
			sum += temp5[i];
		}
		errs[it] = static_cast<double>(sum) / static_cast<double>(beta*(beta - 1));
		if (it >= 1){
			double thresh = errs[it - 1] - errs[it];
			difErr.push_back(thresh);
			Err.push_back(errs[it]);
			if (thresh <= 1.0){
				niternum.push_back(it);
				break;
			}

		}
	}
	for (int i = 0; i < nTemplate; i++){
		H0[i] = H[i];
		firstCheck = true;
	}
	differr.push_back(difErr);
	ERRS.push_back(Err);
	free(temp1); temp1 = NULL;
	free(temp2); temp2 = NULL;
	free(temp3); temp3 = NULL;
	free(temp4); temp4 = NULL;
	free(tempV); tempV = NULL;
	free(tempR); tempR = NULL;
	free(tempV1); tempV1 = NULL;
	free(tempR1); tempR1 = NULL;
	free(tempVR); tempVR = NULL;
	free(temp5); temp5 = NULL;
	free(R); R = NULL;
}
double* MFE::chooseTemplate(double* templateData, int* noteInScore, int noteLength, int* newTemplateColumn){
	int nTemplateSet = originNTemplate / nPitch; // 有几套模板
	int silenceTemplate = originNTemplate % nPitch; // 是否有静音模板

	double *newTemplate;
	if (silenceTemplate == 0) { // 没有静音模板
		*newTemplateColumn = noteLength * nTemplateSet;
		newTemplate = (double*)malloc(sizeof(double) * SPECTRUMLENGTH * (*newTemplateColumn));
	}
	else { // 有静音模板
		*newTemplateColumn = noteLength * nTemplateSet + silenceTemplate;
		newTemplate = (double*)malloc(sizeof(double) * SPECTRUMLENGTH * (*newTemplateColumn));
	}

	int i, j, k;
	for (i = 0; i < nTemplateSet; ++i) {
		for (j = 0; j < noteLength; ++j) {
			int column = noteInScore[j] - 1;
			for (k = 0; k < SPECTRUMLENGTH; ++k) {
				newTemplate[i*nTemplateSet + j + k * (*newTemplateColumn)] = templateData[i*nPitch + column + k * originNTemplate];
			}
		}
	}

	// 添加静音模板
	if (silenceTemplate != 0) {
		for (i = 0; i < SPECTRUMLENGTH; ++i) {
			newTemplate[*newTemplateColumn - 1 + i*(*newTemplateColumn)] = templateData[originNTemplate - 1 + i*originNTemplate];
		}
	}

	// 更新模板列数
	nTemplate = *newTemplateColumn;
	// 更新音符长度
	noteSize = noteLength;

	return newTemplate;

}
double* MFE::getTemplateData(const char *inputPath){
	int i;
	double* templateData = (double *)malloc(sizeof(double)*SPECTRUMLENGTH*nTemplate);
	FILE* templateFile;
	//从template.dat中读入频谱模板template (spectrumLength*nTemplate)
	if ((templateFile = fopen(inputPath, "rb")) == NULL) {
		printf("Cannot open template data file!\n");
		exit(1);
	}
	for (i = 0; i < SPECTRUMLENGTH*nTemplate; ++i) {
		fscanf(templateFile, "%lf", &templateData[i]);
	}
	fclose(templateFile);
	return templateData;
}


void MFE::setscorePitches(vector<vector<double>> pitches){
	scorePitches = pitches;
}


void MFE::setscoreOnset(vector<double> onset){
	scoreOnset = onset;
}

void MFE::getframeOnset(){
	for (int i = 0; i <scoreOnset.size() ; i++){
		frameonset.push_back(round(scoreOnset[i]*44100 /static_cast<double>(512)));
	}
	sort(frameonset.begin(), frameonset.end());
	vector<int>::iterator pos = unique(frameonset.begin(), frameonset.end());
	frameonset.erase(pos, frameonset.end());
}

vector<vector<double>> MFE::getdifferr(){
	return differr;
}

vector<int> MFE::getniternum(){
	return niternum;
}
vector<vector<double>> MFE::geterr(){
	return ERRS;
}

void MFE::GenerateHOffLine(vector<vector<double>> &H, vector<vector<vector<double>>> &sfResult, vector<vector<vector<double>>> &scoreEvent, const string &templatePath, vector<Correctness> &correctness){
	onsetResult = onsetdec.getdecectionResult();
	sfResult = excessLoc(sfResult,scoreEvent);
	changeSfResult(sfResult);
	int sfResultSize = sfResult[0].size();
	H.resize(framerow);
	for (int i = 0; i < framerow; i++){
		H[i].resize(88);
	}
	
	for (int i = 0; i < sfResultSize; i++){
		vector<double> timeNmf;
		for (int j = 0; j < sfResult[0][i].size(); j+=2){
			timeNmf.push_back(sfResult[0][i][j]);
		}
		paramInit();
		vector<double> H0;
		H0.resize(nTemplate);
		vector<int> excession = correctness[i].excess;
		vector<int> NoteInScore = GetNoteInScore(sfResult[1][i][0], sfResult[2][i][0], scoreEvent);
		NoteInScore.insert(NoteInScore.end(), excession.begin(), excession.end());
		double sum = accumulate(timeNmf.begin(), timeNmf.end(), 0.0);
		double mintime = sum / (static_cast<double>(timeNmf.size()));
		vector<vector<double>> OnsetSpectrum = GetSpectrumRange(mintime);
		vector<double> templateData;
		ReadTemplateFile(templatePath, templateData, NoteInScore);
		bool fisrtNmf = true;
		for (int j = 0; j < OnsetSpectrum.size(); j++){
			vector<double> xH = multiDection(&OnsetSpectrum[j][0], &templateData[0], fisrtNmf, &H0[0]);
			for (int k = 0; k <NoteInScore.size(); k++){
				H[startframe+j][NoteInScore[k] - 1] = xH[k];
			}
			fisrtNmf = false;
		}
	}
}


vector<double> MFE::multiDection(double *spectrum, double* templateData, bool fistNmf, double* H0){
	vector<double> xH(noteSize);
	double* H = (double*)calloc(nTemplate, sizeof(double)); // 定义H并初始化为0
	nmf_beta_test(spectrum, SPECTRUMLENGTH, nTemplate, templateData, H, beta, niter,fistNmf,H0);
	for (int i = 1; i < nTemplate / noteSize; ++i) {
		for (int j = 0; j < noteSize; ++j) {
			H[j] += H[j + i * noteSize];
		}
	}
	for (int i = 0; i < noteSize; ++i) xH[i] = H[i];
	free(H);
	H = NULL;
	return xH;
}


vector<vector<vector<double>>> MFE::excessLoc(vector<vector<vector<double>>> &sfResult, vector<vector<vector<double>>> &scoreEvent){
	vector<vector<vector<double>>> ModisfResult(3);
	vector<int> InsertLoc;
	vector<double> diffOnset;
	for (int i = 0; i < onsetResult.size(); i++){
		bool findOnset = false;
		int loc = 0;
		double minmargin = 10000;
		double mindiff = 0.0;
		for (int j = 0; j < sfResult[0].size(); j++){
			double sum = 0.0;
			for (int k = 0; k < sfResult[0][j].size(); k += 2){
				sum += sfResult[0][j][k];
				if (fabs(onsetResult[i] - sfResult[0][j][k]) <= 0.06)
					findOnset = true;
			}
			double difftime = abs(sum / (static_cast<double>(sfResult[0][j].size() / 2)) - onsetResult[i]);
			if (minmargin > difftime){
				loc = j;
				minmargin = difftime;
				mindiff = sum / (static_cast<double>(sfResult[0][j].size() / 2));
			}
		}
		if (!findOnset){
			if (fabs(mindiff - onsetResult[i])>0.2 && loc>0 && loc<sfResult[0].size() - 1){
				int nowloc = static_cast<int>(sfResult[2][loc][0]);
				int lastloc = static_cast<int>(sfResult[2][loc - 1][0]);
				int nextloc = static_cast<int>(sfResult[2][loc + 1][0]);
				if (nowloc - lastloc == 2){
					InsertLoc.push_back(lastloc+1);
					diffOnset.push_back(onsetResult[i]);
				}
				else if (nextloc - nowloc == 2){
					InsertLoc.push_back(nowloc+1);
					diffOnset.push_back(onsetResult[i]);
				}
			}
		}
	}
	for (int i = 0; i < sfResult[0].size(); i++){
		vector<int> lastpitchOctive;
		vector<int> pitchOctive;
		int lastindex = i - 1 >= 0 ? i - 1 : 0;
		for (int j = 1; j < sfResult[0][i].size(); j+=2){
			pitchOctive.push_back(static_cast<int>(sfResult[0][i][j]) % 12);
		}
		for (int j = 0; j < sfResult[0][lastindex].size(); j+=2){
			lastpitchOctive.push_back(static_cast<int>(sfResult[0][lastindex][j]) % 12);
		}
		vector<int> intersection;
		sort(lastpitchOctive.begin(),lastpitchOctive.end());
		sort(pitchOctive.begin(), pitchOctive.end());
		set_intersection(lastpitchOctive.begin(), lastpitchOctive.end(), pitchOctive.begin(), pitchOctive.end(), back_inserter(intersection));
		bool octivequal = (intersection.size() == pitchOctive.size());
		bool findOnset = false;
		vector<double> timeonset;
		for (int j = 0; j < sfResult[0][i].size(); j += 2){
			timeonset.push_back(sfResult[0][i][j]);
		}
		double sum = accumulate(timeonset.begin(), timeonset.end(), 0.0);
		double avgtime = sum / (static_cast<double>(timeonset.size()));
		for (size_t j = 0; j < onsetResult.size(); j++){
			if (fabs(onsetResult[j] - avgtime) <= 0.06){
				findOnset = true;
			}
		}
		if (octivequal&&!findOnset && i!=0) continue;
		else{
			ModisfResult[0].push_back(sfResult[0][i]);
			ModisfResult[1].push_back(sfResult[1][i]);
			ModisfResult[2].push_back(sfResult[2][i]);
		}
		/*vector<int>::iterator pos = find(InsertLoc.begin(), InsertLoc.end(), sfResult[2][i][0] + 1);
		if (pos != InsertLoc.end()){
			vector<double> pitch;
			for (int k = 0; k < scoreEvent[sfResult[2][i][0]][3].size(); k++){
				pitch.push_back(diffOnset[pos - InsertLoc.begin()]);
				pitch.push_back(scoreEvent[sfResult[2][i][0]][3][k]);
			}
			ModisfResult[0].push_back(pitch);
			ModisfResult[1].push_back(vector<double>{static_cast<double>(1)});
			ModisfResult[2].push_back(vector<double>{static_cast<double>(sfResult[2][i][0]+1)});
		}*/
	}
	return ModisfResult;
}


vector<vector<double>> MFE::GetSpectrumRange(double onsetime){
	vector<vector<double>> OnsetSpectrum;
	int frameindex = static_cast<int>(onsetime * 44100 / static_cast<double>(512));

	startframe = frameindex >= 0 ? frameindex : 0;
	endframe = frameindex + 9 <= SpectrumOffLine.size() - 1 ? frameindex + 9 : SpectrumOffLine.size() - 1;
	
	for (int i = startframe; i <=endframe; i++){
		OnsetSpectrum.push_back(SpectrumOffLine[i]);
	}
	return OnsetSpectrum;
}


void MFE::changeSfResult(vector<vector<vector<double>>> &sfResult){
	for (int i = 0; i < sfResult[1].size(); i++){
		if (sfResult[1][i][0] == 0){
			int last = i - 1 >= 0 ? i - 1 : 0;
			int next = i + 1 <= sfResult[1].size() - 1 ? i + 1 : sfResult[1].size() - 1;
			int lastloc = static_cast<int>(sfResult[2][last][0]);
			int nextloc =static_cast<int>(sfResult[2][next][0]);
			int nowloc =static_cast<int>(sfResult[2][i][0]);
			if ((nowloc - lastloc == 1) && (nextloc - nowloc == 1)){
				sfResult[1][i][0] = 1;
			}
		}
	}
}


void MFE::nmf_beta_test(double* data, int rows, int columns, double* W, double* H, double beta, int maxiter, bool fistNmf,double* H0){
	int i, j, k;
	// 初始化H
	//if (fistNmf){
		for (i = 0; i < columns; ++i) {
			H[i] = 1;  // this is for test，the code below is for release
		}
	//}
	//else{
	//	for (i = 0; i < columns; ++i) {
	//		H[i] = H0[i];
	//	}
	//}
	
	// 初始化R
	double* R = (double*)calloc(rows, sizeof(double));
	for (i = 0; i < rows; ++i) {
		for (k = 0; k < columns; ++k) {
			R[i] += W[i * columns + k] * H[k];  // 二维矩阵相乘
		}
	}

	double *temp1 = (double*)calloc(rows, sizeof(double));
	double *temp2 = (double*)calloc(rows, sizeof(double));
	double *temp3 = (double*)calloc(columns, sizeof(double));
	double *temp4 = (double*)calloc(columns, sizeof(double));

	double *tempV = (double*)calloc(rows, sizeof(double));
	double *tempR = (double*)calloc(rows, sizeof(double));
	double *tempV1 = (double*)calloc(rows, sizeof(double));
	double *tempR1 = (double*)calloc(rows, sizeof(double));
	double *tempVR = (double*)calloc(rows, sizeof(double));
	double *temp5 = (double*)calloc(rows, sizeof(double));


	// 迭代
	double myeps = 1.0e-20; // 精度
	int it; // 迭代次数
	for (it = 0; it < maxiter; ++it) {
		//update R
		memset(R, 0, sizeof(double)*rows); // 将R清零
		for (i = 0; i < rows; ++i) {
			for (k = 0; k < columns; ++k) {
				R[i] += W[i * columns + k] * H[k];
			}
		}
		//update H
		// 初始化中间变量，内存清零
		memset(temp1, 0, sizeof(double)*rows);
		memset(temp2, 0, sizeof(double)*rows);
		memset(temp3, 0, sizeof(double)*columns);
		memset(temp4, 0, sizeof(double)*columns);
		for (i = 0; i < rows; ++i) {
			temp1[i] = pow(R[i], beta - 2) * data[i];
			temp2[i] = pow(R[i], beta - 1);
		}
		for (k = 0; k < rows; ++k) {
			for (i = 0; i < columns; ++i) {
				temp3[i] += W[k * columns + i] * temp1[k];
				temp4[i] += W[k * columns + i] * temp2[k];
			}
		}
		for (i = 0; i < columns; ++i) {
			if (temp4[i] < myeps) {
				temp4[i] = myeps;
			}
			H[i] = H[i] * (temp3[i] / temp4[i]);
		}
		//再一次更新update R
		memset(R, 0, sizeof(double)*rows); // 将R清零
		for (i = 0; i < rows; ++i) {
			for (k = 0; k < columns; ++k) {
				R[i] += W[i * columns + k] * H[k];
			}
		}
		for (i = 0; i < rows; i++){
			tempR[i] = pow(R[i], beta);
			tempR[i] = tempR[i] * (beta - 1);
			tempV[i] = pow(data[i], beta);
			tempV1[i] = data[i] * beta;
			tempR1[i] = pow(R[i], (beta - 1));
			tempVR[i] = tempV1[i] * tempR1[i];
			temp5[i] = tempV[i] + tempR[i] - tempVR[i];
		}
		if (it >= 1){
			if (thresh <= 1.0){
				break;
			}

		}
	}
	for (int i = 0; i < nTemplate; i++){
		H0[i] = H[i];
	}
	free(temp1); temp1 = NULL;
	free(temp2); temp2 = NULL;
	free(temp3); temp3 = NULL;
	free(temp4); temp4 = NULL;
	free(tempV); tempV = NULL;
	free(tempR); tempR = NULL;
	free(tempV1); tempV1 = NULL;
	free(tempR1); tempR1 = NULL;
	free(tempVR); tempVR = NULL;
	free(temp5); temp5 = NULL;
	free(R); R = NULL;
}
vector<int> MFE::GetNoteInScore(int SureFlag, int location, vector<vector<vector<double>>> &scoreEvent){
	vector<int> NoteInScore;
	int start = 0;
	int end = 0;

	if (SureFlag == 1){
		start = end = location - 1;
	}
	else{
		start = location - 1 >= 0 ? location - 1 : 0;
		end = location + 1 <= scoreEvent.size() - 1 ? location + 1 : scoreEvent.size() - 1;
	}
	int bitmap[88] = { 0 };
	for (vector<int>::size_type i = start; i <= end; ++i) {
		for (vector<int>::size_type j = 0; j < scoreEvent[i][3].size(); ++j) {
			int temp = static_cast<int>(scoreEvent[i][3][j]);
			++bitmap[temp];
		}
	}
	for (int i = 0; i < 88; ++i) {
		if (bitmap[i] != 0) {
			NoteInScore.push_back(i);
		}
	}
	set<int> tempresult;
	vector<int> result;
	if (!NoteInScore.empty() && SureFlag==0){
		for (int i = 0; i < NoteInScore.size(); i++){
			int begin = NoteInScore[i];
			int end = NoteInScore[i];
			begin = begin > 0 ? begin : 1;
			end = end < 89 ? end : 88;
			for (int j = begin; j <= end; j++){
				tempresult.insert(j);
			}
		}
	}
	else if (!NoteInScore.empty() && SureFlag == 1){
		for (int i = 0; i < NoteInScore.size(); i++){
			int begin = NoteInScore[i];
			int end = NoteInScore[i];
			begin = begin > 0 ? begin : 1;
			end = end < 89 ? end : 88;
			for (int j = begin; j <= end; j++){
				tempresult.insert(j);
			}
		}
	}
	for (set<int>::iterator it = tempresult.begin(); it != tempresult.end(); it++){
		result.push_back(*it);
	}
	sort(result.begin(), result.end());
	NoteInScore = result;
	return NoteInScore;
}

vector<vector<double>> MFE::GetSpectrum(){
	return SpectrumOffLine;
}
void MFE::saveOnset(const string &path){
	onsetdec.saveResult(path);
}
