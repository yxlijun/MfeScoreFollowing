#coding:utf--8

import os 
import argparse 
import shutil

parser = argparse.ArgumentParser()
parser.add_argument('--folderpath','-f',type=str,help='folderpath')
args = parser.parse_args()

filelist = [os.path.join(args.folderpath,x) for x in os.listdir(args.folderpath) if x.endswith(".wav")]

for file in filelist:
	filepath,_ext = os.path.splitext(file)
	if not os.path.exists(filepath):
		os.makedirs(filepath)
	shutil.move(file,filepath)

print 'done!'