#pragma once
#ifndef MFE_H_
#define MFE_H_
#include <vector>
#include "kiss_fft.h"
#include "resample.h"
#include "Timer.h"
#include <fstream>
#include <string>
#include <algorithm>
#include "Onsets.h"
#include <set>
#include "EvaluateSfResult.h"
using namespace std;
class MFE
{
public:
	/**
	* @author zhangqianyi
	* @date 九月 2016
	* @brief 多音调检测接口函数，输入原始录音数据，读入频谱模板处理之后调用多音调检测算法返回pianoRoll
	* @param  short* originWav 输入原始录音数据，长度固定为4410，但是使用只取前4096个数据.
	* @param  int wChannels 输入原始录音数据的声道数，默认为单声道数据.
	* @return int* pianoRoll 输出的结果.
	*/
	int* MFEinterface(short* originWav, int wChannels, double* templateData, double* xH, double * error, const string &path,double *ratio);   //对每一帧进行多音调检测
	/**
	* @author zhangqianyi
	* @date 七月 2016
	* @brief 对一帧数据进行处理，得出一帧数据的pianoRoll值.
	* @param  xFrame 输入的一帧数据.
	* @param  templateData 输入的模板数据.
	* @param  pianoRoll 输出的结果.
	* @return 返回value.
	*/
	void multiF0Estimation(double* xFrame, double* templateData, int &H0Flag, double* H0, double* xH, double * error, const string &path, double *ratio, double *specrum, double*OutFrame, Onset &onsetdec);
	/**
	* @author zhangqianyi
	* @date 2016/7/12
	* @brief 对一帧信号进行FFT，计算振幅谱。 窗函数为Hamming窗，DFT偶对称；若帧长小于fft长度，信号加窗后末尾补0。
	* @param  xFrame 一帧信号（长度为帧长。若帧长为奇数，窗函数不满足DFT偶对称）
	* @return 振幅谱地址（前半段, fftSize/2+1）
	*/
	void spectrumOfFrame(double* xFrame, double* spectrum);
	/**
	* @author zhangqianyi
	* @date 七月 2016
	* @brief 获得一个指定大小的hamming窗.
	* @param  w 返回的窗函数值.
	* @param  windowSize 窗函数的大小。
	*/
	void hammingWindow(double *w, int windowSize);

	/**
	* @author zhangqianyi
	* @date 七月 2016
	* @brief nmf_beta非负矩阵分解函数（data = W * H，已知data和W来得到H）.
	* @param data：要分解的矩阵
	* @param rows：W的行数
	* @param classes：W的列数
	* @param W：矩阵分解的频谱模板
	* @param H：矩阵分解的结果
	* @param beta：用于矩阵分解的beta参数
	* @param maxiter：循环计算的次数
	* @param H0Flag：是否有H0输入，1为有H0，0为没有H0
	* @param H0：H0数据
	*/
	void nmf_beta(double* data, int rows, int column, double* W, double* H, int &H0Flag, double* H0, double* errs, double beta, int maxiter,bool findonset);
	/**
	* @author zhangqianyi
	* @date 四月 2017
	*
	* @param templateData  完整的模板。要求：一套模板对应NPITCH个音符，列号为音符序号；多套模板拼接，各套对应的音符相同；若有静音模板，置于末尾，静音模板数<NPITCH
	* @param noteInScore   乐谱包含的音符的序号
	* @param noteLength      noteInScore的长度
	* @param newTemplateColumn 返回模板的列数
	* @return newTemplate  乐谱包含的音符对应的模板。包含静音模板
	* @brief 根据乐谱包含的音符选择模板
	*/
	double* chooseTemplate(double* templateData, int* noteInScore, int noteLength, int* newTemplateColumn);
	/**
	* @author zhangqianyi
	* @date   2016/12/6
	* @param  inputPath 表示模板文件路径
	* @return 返回保存模板数据的动态分配double型数组指针
	* @brief  读取模板数据
	*/
	double* getTemplateData(const char *inputPath);

	void QSort(double* a, int low, int high);

	/**
	设置是否进行重采样操作
	*/
	void setresample(bool resample);

	/**
	@date  2017/7/11
	@param  framePath  表示乐谱frame文件路径
	@param	xFrame 将文件数据保存到xFrame
	@brief 在纯离线时读取整个乐谱frame数据
	*/
	void  ReadFrameFile(const string &framePath, vector<vector<double>>& xFrame);
	/*
	@date 2017/7/11
	@param  templatePath 表示模板文件路径
	@param templateData 模板数据
	@brief 读取模板文件
	*/
	void ReadTemplateFile(const string &templatePath, vector<double>& templateData, vector<int> &notesInScore);

	/*
	@param  H   定位H ，列数为88
	@param  xFrame  乐谱frame数据
	@param   templateData  模板数据
	@brief 生成用于定位的H
	*/
	void GenerateH(vector<vector<double>> &H, vector<double> &error, vector<vector<double>> &ratios, vector<vector<double>> &xFrame, vector<double> &templateData, const string &path, vector<int> &notesInScore, vector<vector<double>> &Specurm, vector<vector<double>> &outframe);



	int getrow();

	void setscorePitches(vector<vector<double>> pitches);

	void setscoreOnset(vector<double> onset);

	void getframeOnset();

	vector<vector<double>> getdifferr();
	vector<int> getniternum();
	vector<vector<double>> geterr();

	vector<vector<double>> GetSpectrum();


	/*
	以下函数是测试函数，在弹奏完成后进行二次NMF分解
	*/
	void GenerateHOffLine(vector<vector<double>> &H, vector<vector<vector<double>>> &sfResult, vector<vector<vector<double>>> &scoreEvent, const string &templatePath, vector<Correctness>& correctness);
	vector<int> GetNoteInScore(int SureFlag, int location,vector<vector<vector<double>>> &scoreEvent);
	vector<vector<double>> GetSpectrumRange(double onsetime);
	void nmf_beta_test(double* data, int rows, int columns, double* W, double* H, double beta, int maxiter, bool fistNmf, double* H0);
	vector<double> multiDection(double *spectrum, double* templateData, bool fistNmf,double* H0);
	vector<vector<vector<double>>> excessLoc(vector<vector<vector<double>>> &sfResult, vector<vector<vector<double>>> &scoreEvent);
	void changeSfResult(vector<vector<vector<double>>> &sfResult);

	void saveOnset(const string &path);


	MFE();
	~MFE();
	void paramInit();

private:
	const int originNTemplate;	// 初始template模板的列数
	int nTemplate;               // template模板的列数
	int noteSize;                // 乐谱中音符个数，读取完乐谱之后确定
	int maxPitches;				// 一帧数据中最多出现的音符个数
	const int nPitch = 88;            // 音符个数，这个是常数88

	int FFTSIZE;                      // fft长度
	int SPECTRUMLENGTH;              // 谱长度（fftSize/2+1）
	int HOPSIZE;                    // 帧移
	int FRAMESIZE;                  // 帧长（samples）

	double minThresh;     // 阈值最小值(真钢琴)
	double threshCoeff;   // 阈值系数(真钢琴)
	const double beta = 0.6;           //nmf_beta算法中beta
	int niter;          //nmf_beta算法中迭代次数
	bool isresample;

	int framerow;       // frame文件的行数
	int framecol;       //frame文件的列数

	double thresh;   //nmf跟新过程的loss阈值

	vector<vector<double>> scorePitches;   // 乐谱音符
	vector<double> scoreOnset;            // 乐谱音符onset
	vector<int> frameonset;
	
	vector<vector<double>> differr;
	vector<vector<double>> ERRS;
	vector<int> niternum;
	vector<int> notesScore;


	vector<vector<double>> SpectrumOffLine;
	Onset onsetdec;


	int startframe;
	int endframe;
	bool firstCheck;
	vector<double> onsetResult;

};


#endif