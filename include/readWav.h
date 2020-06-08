//
//  readWav.h
//  XYAudioFollowing
//
//  Created by zxy on 16/12/9.
//  Copyright © 2016年 zxy. All rights reserved.
//

#ifndef readWav_h
#define readWav_h

#include <stdio.h>
typedef struct
{
    double *normalizedData; // 存放数据的指针
    int length; // 数据长度
    int fs; // 采样频率
    int wChannels; // 声道数
}dataInfo;


/**
 * @author zhangqianyi
 * @date 2016/7/12
 * @brief 读取.wav文件，根据wav文件的存储格式读取wav文件中保存的数据.
 * @param  inputPath 文件的存放路径.
 * @return 返回dataInfo的结构体类型，包括存放数据的指针，数据长度，采样频率以及声道数（如果为双声道的话数据存放形式是一个数据左声道，接着下一个是右声道数据）
 */

dataInfo readWav(char *inputPath);


double* wavDataProcess(double* data, int len, int wChannels);




#endif /* readWav_h */
