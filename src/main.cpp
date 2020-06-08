#include <iostream>
#include "ScoreFollowing.h"
#include "EvaluateSfResult.h"
#include "MFE.h"
#include "Timer.h"
#include <io.h>
#include "tool.h"
using namespace std;

vector<string> datfiles;
vector<string> audiofiles;
vector<string> datpaths;

void saveH(vector<vector<double>> &xH,string Hpath);     // 保存H
void saveRatio(vector<vector<double>> ratio, string RatioPath);        //保存滤掉人声的比例
void saveOutFrame(vector<vector<double>> &outframe, int row, string OutresamplePath);    //保存输出帧
void saveSpectrum(vector<vector<double>> &Spectrum, string SpectrumPath);    //保存频谱
void GetAllFileFromPath(string folderPath);   //遍历文件夹下的所有音频数据
void ReadOnset(const string& path, vector<vector<double>> &scorePitches, vector<double> &scoreOnset);
void saveDifferr(vector<vector<double>> &differr, int row, string differpath);
void saveErr(vector<vector<double>> &err, int row, string errpath);
void saveNiter(vector<int> &niternum, string niterpath);
void ReadOnsetandPitch(const string &path, vector<vector<double>> &scorePitches, vector<double> &scoreOnset);
//void findReplay(vector<vector<vector<double>>> &sfResult, vector<vector<vector<double>>> &scoreEvent);

const string windowpath = "../data/config/window4096_1024.txt";
const string FilePath = "../data/datasets/2017";
const string OnsetPath = "../data/midi/Fly/KuDong.txt";
const string pitchFreFile = "../data/config/pitchFre.txt";


int main()
{
	vector<vector<double>> scorePitches;
	vector<double> scoreOnset;
	//ReadOnsetandPitch(OnsetPath,scorePitches,scoreOnset);
	GetAllFileFromPath(FilePath);
	initFrame(FilePath);
	string arraysavefile[7] = {"/sfResult.txt","/sfResultOrigin.txt","/evaluateResult.txt","/evaluateResultOrigin.txt","/beatRhythm.json","/correctness.json","/playback.json"};
	for (int i = 0; i < datfiles.size(); i++){
		const string sfResultPath = datpaths[i] + arraysavefile[0];
		const string sfResultOriginPath = datpaths[i] + arraysavefile[1];
		const string evaluateResultPath = datpaths[i] + arraysavefile[2];
		const string evaluateResultOriginPath = datpaths[i] + arraysavefile[3];
		const string beatRhythmPath = datpaths[i] + arraysavefile[4];
		const string correctnessPath = datpaths[i] + arraysavefile[5];
		const string playbackPath = datpaths[i] + arraysavefile[6];
		cout << datfiles[i] << endl;
		bool isres = true;
		MFE mfe;
		mfe.setscorePitches(scorePitches);
		mfe.setscoreOnset(scoreOnset);
		ScoreFollowing scoreFollowing;
		if (scoreFollowing.Init(datfiles[i]) == -1) {
			cout << "Read scoreEvent error!" << endl;
			exit(-1);
		}
		scoreFollowing.processNMF.SetThreshParams(0.2, 240);
		double timeResolution = static_cast<double>(512) / 44100;

		vector<int> notesInScore = scoreFollowing.GetNoteInScore();
		int maxPitchesInEvent = scoreFollowing.GetMaxPitchesInFrame();
		mfe.setresample(isres);
		string templatePath;    // 模板路径
		if (isres){
			templatePath = "../data/config/template11025.txt";
		}
		else{
			templatePath = "../data/config/template44100.txt";
		}
		vector<vector<double>> xFrame;
		vector<double> templateData;
		vector<vector<double> > xH;
		vector<double> Error;
		vector<vector<double>> Ratios;
		vector<vector<double>> Specurm;
		vector<vector<double>> outframe;
		
		string framefile = datpaths[i] + "/frame.txt";
		mfe.ReadFrameFile(framefile, xFrame);
		mfe.ReadTemplateFile(templatePath, templateData, notesInScore);
		
		mfe.GenerateH(xH, Error, Ratios, xFrame, templateData, windowpath, notesInScore, Specurm, outframe);
		Timer tim;
		tim.Tic();
		scoreFollowing.ScoreFollowingOffline(xH, Error, Ratios, timeResolution, maxPitchesInEvent);
		tim.Toc();
		cout << "total time:" << tim.Elasped() << "ms" << endl;

		vector<int> niternum = mfe.getniternum();
		string niterPath = datpaths[i] + "/niter.txt";
		cout << accumulate(niternum.begin(), niternum.end(), 0) << endl;
		saveNiter(niternum, niterPath);
		string Hpath = datpaths[i] + "/H.txt";
		saveH(xH,Hpath);
		vector<vector<double>> SpectrumOffLine = mfe.GetSpectrum();
		string spectrumPath = datpaths[i] + "/spectrum.txt";
		saveSpectrum(SpectrumOffLine, spectrumPath);
		/*vector<vector<double>> differr = mfe.getdifferr();
		vector<vector<double>> errs = mfe.geterr();
		vector<int> niternum = mfe.getniternum();
		string DiffPath = datpaths[i] + "/thresh.txt";
		string errpath = datpaths[i] + "/err.txt";
		string niterpath = datpaths[i] + "/niter.txt";
		saveDifferr(differr, mfe.getrow(), DiffPath);
		saveErr(errs, mfe.getrow(), errpath);
		saveNiter(niternum, niterpath); */

		
		SaveSfResult(sfResultPath, scoreFollowing.GetSfResult());
		SaveSfResult(sfResultOriginPath, scoreFollowing.GetSfResultOrigin());
		// 乐谱结果评价
		//scoreFollowing.findReplay();

	/*	string sfResultPath1 = datpaths[i] + "/sfResult1.txt";
		string sfResultOriginPath1 = datpaths[i] + "/sfResultOrigin1.txt";

		SaveSfResult(sfResultPath1, scoreFollowing.GetSfResult());
		SaveSfResult(sfResultOriginPath1, scoreFollowing.GetSfResultOrigin());*/

		EvaluateSfResult evaluateResult(scoreFollowing);
		evaluateResult.Init();
		evaluateResult.SetSpectrum(SpectrumOffLine);
		evaluateResult.ReadpitchFre(pitchFreFile);
		vector<Correctness> correctness = evaluateResult.EvaluateCorrectnessModify();
		vector<BeatRhythm> beatRhythm = evaluateResult.EvaluateBeatRhythm();
		//vector<PitchRhythm> PitchRhythm = evaluateResult.EvaluatePerNoteRhythm();
		vector<PitchRhythm> PitchRhythm;
		double score = evaluateResult.CountScore(correctness, PitchRhythm);
		int starsNum = evaluateResult.GiveStars(score);
		cout << "score = " << score << " stars = " << starsNum << endl;
		evaluateResult.SaveEvaluateResult(evaluateResultPath, correctness, beatRhythm);
		vector<Correctness> correctnessOrigin = evaluateResult.EvaluateCorrectnessOrigin(maxPitchesInEvent);

		evaluateResult.SaveEvaluateResult(evaluateResultOriginPath, correctnessOrigin, beatRhythm);
		evaluateResult.SaveBeatRhythm(beatRhythmPath, beatRhythm);
		evaluateResult.SaveCorrectness(correctnessPath, correctness);
		evaluateResult.SavePlayBack(playbackPath);
		string starFile = datpaths[i] + "/" + to_string(starsNum) + "star.txt";
		cout << starFile << endl;
		ofstream ofs;
		ofs.open(starFile);
		ofs.close();
		
		string onsetPath = datpaths[i] + "/onset.txt";
		mfe.saveOnset(onsetPath);
		/*string pitchDurationPath = datpaths[i] + "/durationResult.txt";
		evaluateResult.SaveDuratinTime(pitchDurationPath);*/
		/*vector<vector<vector<double>>> sfResult = scoreFollowing.GetSfResult();
		vector<vector<vector<double>>> scoreEvent = scoreFollowing.GetScoreEvent();
		vector<vector<double>> OfflineH;
		string onsetPath = datpaths[i] + "/onset.txt";
		mfe.saveOnset(onsetPath);

		Timer time1;
		time1.Tic();
		mfe.GenerateHOffLine(OfflineH, sfResult, scoreEvent, templatePath, correctnessOrigin);
		scoreFollowing.Init(datfiles[i]);
		scoreFollowing.processNMF.SetThreshParams(0.2, 240);
		scoreFollowing.ScoreFollowingOffline(OfflineH, Error, Ratios, timeResolution, maxPitchesInEvent);
		scoreFollowing.findReplay();


		SaveSfResult(sfResultPath, scoreFollowing.GetSfResult());
		SaveSfResult(sfResultOriginPath, scoreFollowing.GetSfResultOrigin());

		EvaluateSfResult evaluateResult1(scoreFollowing);
		evaluateResult1.Init();
		vector<Correctness> correctness1 = evaluateResult1.EvaluateCorrectnessModify();
		vector<BeatRhythm> beatRhythm1 = evaluateResult1.EvaluateBeatRhythm();
		vector<int> noteRhythm1 = evaluateResult1.EvaluateNoteRhythm();
		double score1 = evaluateResult1.CountScore(correctness1, beatRhythm1);
		int starsNum1 = evaluateResult1.GiveStars(score1);
		cout << "scoreModify = " << score1 << " starsModify = " << starsNum1 << endl;


		evaluateResult1.SaveEvaluateResult(evaluateResultPath, correctness1, beatRhythm1);

		vector<Correctness> correctnessOrigin1 = evaluateResult1.EvaluateCorrectnessOrigin(maxPitchesInEvent);
		evaluateResult1.SaveEvaluateResult(evaluateResultOriginPath, correctnessOrigin1, beatRhythm1);

		time1.Toc();
		cout << "second total time:" << time1.Elasped() << "ms" << endl;

		string Hpath1 = datpaths[i] + "/H2.txt";
		saveH(OfflineH, Hpath1);
		*/
		remove(framefile.c_str());
		Genexcel(datpaths[i]);
	}
	return 0;
}




void saveH(vector<vector<double>> &xH,string Hpath){
	FILE *HFILE;
	HFILE = fopen(Hpath.c_str(), "w+");
	for (int i = 0; i < xH.size(); i++){
		for (int j = 0; j < 88; j++){
			fprintf(HFILE, "%lf\t", xH[i][j]);
		}
		fprintf(HFILE, "\n");
	}
	fclose(HFILE);
}

void saveDifferr(vector<vector<double>> &differr, int row, string differpath){
	FILE *HFILE;
	HFILE = fopen(differpath.c_str(), "w+");
	for (int i = 0; i < row; i++){
		for (int j = 0; j <differr[i].size() ; j++){
			fprintf(HFILE, "%lf\t", differr[i][j]);
		}
		fprintf(HFILE, "\n");
	}
	fclose(HFILE);
}

void saveErr(vector<vector<double>> &err, int row, string errpath){
	FILE *HFILE;
	HFILE = fopen(errpath.c_str(), "w+");
	for (int i = 0; i < row; i++){
		for (int j = 0; j <err[i].size(); j++){
			fprintf(HFILE, "%lf\t", err[i][j]);
		}
		fprintf(HFILE, "\n");
	}
	fclose(HFILE);
}

void saveNiter(vector<int> &niternum,string niterpath){
	FILE *HFILE;
	int row = niternum.size();
	HFILE = fopen(niterpath.c_str(), "w+");
	for (int i = 0; i < row; i++){
		fprintf(HFILE, "%d\n", niternum[i]);
	}
	fclose(HFILE);
}


void saveRatio(vector<vector<double>> ratio,string RatioPath){
	FILE *RFILE;
	RFILE = fopen(RatioPath.c_str(), "w+");
	for (int i = 0; i < ratio.size(); i++){
		fprintf(RFILE, "%f\t%f\t%f\n", ratio[i][0],ratio[i][1],ratio[i][2]);
	}
	fclose(RFILE);
}

void saveSpectrum(vector<vector<double>> &Spectrum, string SpectrumPath){
	FILE *SpecFILE;
	SpecFILE = fopen(SpectrumPath.c_str(), "w+");
	for (int i = 0; i < Spectrum.size(); i++){
		for (int j = 0; j < 513; j++){
			fprintf(SpecFILE, "%lf\t", Spectrum[i][j]);
		}
		fprintf(SpecFILE, "\n");
	}
	fclose(SpecFILE);
}


void saveOutFrame(vector<vector<double>> &outframe, int row,string OutresamplePath){
	FILE *OutFILE;
	OutFILE = fopen(OutresamplePath.c_str(), "w+");
	for (int i = 0; i < row; i++){
		for (int j = 0; j < 1024; j++){
			fprintf(OutFILE, "%lf\t", outframe[i][j]);
		}
		fprintf(OutFILE, "\n");
	}
	fclose(OutFILE);
}



void GetAllFileFromPath(string folderPath)
{
	_finddata_t FileInfo;
	string strfind = folderPath + "/*";
	intptr_t  Handle = _findfirst(strfind.c_str(), &FileInfo);

	if (Handle == -1L){
		_findclose(Handle);
		return;
	}
	do{
		if (FileInfo.attrib & _A_SUBDIR){
			if ((strcmp(FileInfo.name, ".") != 0) && (strcmp(FileInfo.name, "..") != 0)){
				string newPath = folderPath + "/" + FileInfo.name;
				GetAllFileFromPath(newPath);
			}
		}
		else{
			string newPath = folderPath + "/" + FileInfo.name;
			string path = FileInfo.name;
			int index = path.find(".dat", 0);
			int index2 = path.find("mp3", 0);
			if (index != -1){
				datfiles.push_back(newPath);
				datpaths.push_back(folderPath);
			}
			if (index2 != -1){
				audiofiles.push_back(newPath);
			}
		}
	} while (_findnext(Handle, &FileInfo) == 0);

	_findclose(Handle);
}


void ReadOnset(const string &path, vector<vector<double>> &scorePitches,vector<double> &scoreOnset){
	vector<double> scorePitch;
	ifstream OnsetFin(path);
	string line;
	int count = 0;
	while (OnsetFin>>line){
		if (count % 3 == 0){
			scorePitch.push_back(stod(line));
		}
		else if (count % 3 == 1){
			scoreOnset.push_back(stod(line));
		}
		count++;
	}
	for (int i = 0; i < scoreOnset.size(); i++){
		int j;
		vector<double> pitch;
		for (j = i; j < scoreOnset.size(); j++){
			if (scoreOnset[i] == scoreOnset[j]){
				pitch.push_back(scorePitch[j]);
			}
			else break;
		}
		scorePitches.push_back(pitch);
		i = j;
		i--;
	}
}

void ReadOnsetandPitch(const string &path, vector<vector<double>> &scorePitches, vector<double> &scoreOnset){
	vector<double> scorePitch;
	ifstream OnsetFin(path);
	string line;
	int count = 0;
	while (getline(OnsetFin,line)){
		vector<string> info = split(line, "\t");
		vector<double> newinfo;
		vector<double> pitch;
		vector<double> onset;
		for (int i = 0; i < info.size(); i++){
			if (info[i] != ""){
				newinfo.push_back(stod(info[i]));
			}
		}
		for (int i = 0; i < newinfo.size()/2; i++){
			pitch.push_back(newinfo[i]);
			onset.push_back(newinfo[i + newinfo.size() / 2]);
		}
		if (onset.empty())
			break;
		double minonset = *min_element(onset.begin(), onset.end());
		scoreOnset.push_back(minonset);
 		scorePitches.push_back(pitch);
	}
	return;
}

//void findReplay (vector<vector<vector<double>>> &sfResult, vector<vector<vector<double>>> &scoreEvent){
//	vector<int> loctionOffline;
//	vector<int> ReplayLoc;
//	for (int i = 0; i < sfResult[2].size()-1; i++){
//		int diff = sfResult[2][i + 1][0] - sfResult[2][i][0];
//		if (diff == 0 || diff < 0){
//			ReplayLoc.push_back(i + 1);
//		}
//	}
//	for (int i = 0; i < ReplayLoc.size(); i++){
//		int temp = 0;
//		int temp1 = 0;
//		for (int j = i + 1; j < ReplayLoc.size(); j++){
//			if (ReplayLoc[j] - ReplayLoc[i] <= (j - i + 1)){
//				temp = ReplayLoc[j];
//				temp1 = ReplayLoc[i];
//			}
//			else{
//				if (j - i <= 2) temp = 0;
//				i = j;
//				i--;
//				break;
//			}
//		}
//		if (temp != 0){
//			loctionOffline.push_back(temp1 - 1);
//			loctionOffline.push_back(temp);
//		}
//	}
//	for (int i = 0; i < loctionOffline.size(); i+=2){
//		vector<vector<int>> ReplayPitch;
//		vector<vector<double>> ScoreReplayPitch;
//		for (int j = loctionOffline[i]; j <=loctionOffline[i+1]; j++){
//			vector<int> pitch;
//			for (int k = 1; k < sfResult[0][j].size(); k+=2){
//				pitch.push_back(sfResult[0][j][k]);
//			}
//			ReplayPitch.push_back(pitch);
//		}
//		int lastLoc = sfResult[2][loctionOffline[i]][0];
//		int from = lastLoc - 12 >= 0 ? lastLoc - 12 : 0;
//		int to = lastLoc + 1 <= scoreEvent.size() - 1 ? lastLoc + 1 : scoreEvent.size() - 1;
//		for (int j = from; j <=to; j++){
//			ScoreReplayPitch.push_back(scoreEvent[j][g_Pitches]);
//		}
//		vector<double> matchSum;
//		for (int k = 0; k < ScoreReplayPitch.size(); k++){
//			if (k + ReplayPitch.size() <= ScoreReplayPitch.size()){
//				vector<double> matchDegree;
//				for (int j = 0; j < ReplayPitch.size(); j++){
//					matchDegree.push_back(matchIEvent(ReplayPitch[j], scoreEvent, from+k+j));
//				}
//				matchSum.push_back(accumulate(matchDegree.begin(), matchDegree.end(), 0.0));
//			}
//		}
//		int maxMatch = *max_element(matchSum.begin(), matchSum.end());
//		for (int j = matchSum.size()-1; j >=0; j--){
//			if (matchSum[j] == maxMatch){
//				maxMatch = from + j + 1;
//				break;
//			}
//		}
//		for (int j = loctionOffline[i]; j <= loctionOffline[i + 1]; j++){
//			sfResult[2][j][0] = maxMatch + j - loctionOffline[i];
//			//cout << maxMatch + j - loctionOffline[i] << endl;
//		}
//		//copy(matchSum.begin(), matchSum.end(), ostream_iterator<double>(cout,"\n"));
//	}
//}