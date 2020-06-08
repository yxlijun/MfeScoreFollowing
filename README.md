#本工程将多音调检测和音调定位合并为一个工程

用cmake编译便可以直接运行该工程，在工程目录下新建一个build目录   
mkdir build  
cd build    
cmake ..     

工程->属性->c/c++->预处理器->预处理器定义  添加      
_SCL_SECURE_NO_WARNINGS     
_CRT_SECURE_NO_WARNINGS     


/data 提供该工程需要的数据文件   
/src  工程所有源文件    
/include 工程所有头文件     



