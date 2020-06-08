//
//  readWav.c
//  XYAudioFollowing
//
//  Created by zxy on 16/12/9.
//  Copyright © 2016年 zxy. All rights reserved.
//

#include "readWav.h"
#include <stdlib.h>

// 下面结构体的定义存放在实现函数文件中，是读取音频文件独有的结构，不用存放在.h文件中
struct RIFF_HEADER
{
    //RIFF_HEADER
    unsigned   char     szRiffID[4]; //资源交换文件标志'R','I','F','F'
    unsigned   int     dwRiffSize; //从下个地址开始到文件尾的总字节数
    unsigned   char     szRiffFormat[4]; //wav文件标志'W','A','V','E'
}riff_header;

struct FMT_BLOCK
{
    //FMT_BLOCK
    unsigned   char     szFmtID[4]; //波形格式标志'f','m','t',' '
    unsigned   int     dwFmtSize; //过滤字节一般为10h
    unsigned   short    wFormatTag; //格式种类（值为1时，表示数据为线性pcm编码）
    unsigned   short    wChannels; //通道数，单声道为1，双声音为2
    unsigned   int     dwSamplesPerSec; //采样频率
    unsigned   int     dwAvgBytesPerSec; //数据传输率 (每秒字节＝采样频率×每个样本字节数)
    unsigned   short    wBlockAlign; // 块对齐字节数 = channles * bit_samp / 8
    unsigned   short    wBitsPerSample; //样本数据位数(量化位数)
    
}fmt_block;

struct FACT_BLOCK
{
    unsigned   char     szFactID[4];   //'f','a','c','t'
    unsigned   int     dwFactSize;
}fact_block;

struct DATA_BLOCK
{
    //DATA_BLOCK
    unsigned   char   szDataID[4];   //'d','a','t','a'
    unsigned   int   dwDataSize;
}data_block;


dataInfo readWav(char *inputPath)
{
    FILE *fp;
    if ((fp = fopen(inputPath, "rb")) == NULL) // 为读，打开一个wav文件,若打开文件失败，退出
    {
        printf("can't open this file\n");
        exit(0);
    }
    
    fread(&riff_header, sizeof(struct RIFF_HEADER), 1, fp); //riff_header的读取示
    
    fread(&fmt_block, sizeof(struct FMT_BLOCK), 1, fp); //fmt_block的读取
    int fs = fmt_block.dwSamplesPerSec;
    
    if (fmt_block.wFormatTag != 1) //如果读取的wav文件不是pcm编码，还需进行读取fact块的操作
    {
        unsigned short cbsize;
        fread(&cbsize, sizeof(unsigned short), 1, fp); //先要读入两个字节的cbSize
        fread(&fact_block, sizeof(struct FACT_BLOCK), 1, fp); //fact_block的读取
        unsigned int dwSampleLength;
        fread(&dwSampleLength, sizeof(unsigned int), 1, fp); //最后还要读入四个字节的dwSampleLength
    }
    int wBytesPerSample = fmt_block.wBitsPerSample / 8; //样本每个数据是多少字节
    int wChannels = fmt_block.wChannels; //声道数
    fread(&data_block, sizeof(struct DATA_BLOCK), 1, fp); //data_block的读取
    int length = data_block.dwDataSize / wBytesPerSample; //采样的总数据个数（16位一个）
    
    short *dataContent = (short *)malloc(length * sizeof(short));
    fread(dataContent, sizeof(short), length, fp);
    double *normalizedData = (double*)malloc(length * sizeof(double));
    int i = 0;
    for (i = 0; i < length; i++)
    {
        normalizedData[i] = (double)dataContent[i] / 32767.0000;
    }
    free(dataContent);
    dataContent = NULL;
    
    fclose(fp);
    
    dataInfo data_info;
    data_info.fs = fs;
    data_info.length = length;
    data_info.wChannels = wChannels;
    data_info.normalizedData = normalizedData;
    
    return data_info;
}

/**
 * @author zhangqianyi
 * @date 八月 2016
 
 * @param  data 多声道数据值.
 * @param  len 多声道数据长度（多个声道算一个数据）.
 * @param  wChannels 声道数.
 * @return 返回处理之后的数据.
 * @brief 对readWav函数读取的音频数据进行处理，将多声道音频平均化为一个声道或者只取一个声道
 */
double* wavDataProcess(double* data, int len, int wChannels) {
    
    int i;
    double* processedData = (double*)calloc(len, sizeof(double));
    
    //多声道平均化一个声道
    for (i = 0; i < len; ++i) {
        int j = 0;
    	for (j = 0; j < wChannels; ++j) {
    		processedData[i] += data[i * wChannels + j];
    	}
    	processedData[i] /= wChannels;
    }
    
    // //只取一个声道
    // for (i = 0; i < len;i++)
    // {
    //     processedData[i] = data[i * wChannels];
    // }
    
    return processedData;
}


