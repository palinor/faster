#!/bin/zsh
logging=0
debugMode=1
compiler='/usr/bin/clang++'
workingDir='/Users/aionfeehan/Github/analytics/'
buildDir=${workingDir}debug/
logFile='builddebug.log'
libraryName='analytics.so'
frameworksDir="/Applications/Xcode.app/Contents/Developer/Platforms/MacOSX.platform/Developer/SDKs/MacOSX.sdk/System/Library/Frameworks/"

if [ $debugMode -gt 0]; then
optionList=('-std=gnu++20')
optionList+='-fcolor-diagnostics'
optionList+='-fansi-escape-codes'
optionList+='-Wall'
optionList+='-Wextra'
optionList+='-pedantic'
optionList+='-DEBUG'
optionList+='-O0'
optionList+='-g'
else
optionList=('-fcolor-diagnostics')
optionList+='-fansi-escape-codes'
optionList+='-Ofast'
optionList+='-std=gnu++20'
fi

#baseFiles=('linalg')
#baseFiles+='matrix'
baseFiles+=('basic_math')
baseFiles+='threadpool_better'
#baseFiles+='matrix_benchmarks'
baseFiles+='tests'

imguiDir='imgui/'
imguiBaseFiles=('imgui')
imguiBaseFiles+='imgui_demo'
imguiBaseFiles+='imgui_draw'
imguiBaseFiles+='imgui_tables'
imguiBaseFiles+='imgui_widgets'
imguiBaseFiles+='backends/imgui_impl_glut'
imguiBaseFiles+='backends/imgui_impl_opengl2'

imguiOptions=('-L/usr/local/lib')
imguiOptions+='-L/opt/local/lib'
imguiOptions+='-I/usr/local/include'
imguiOptions+='-I/opt/local/include'
imguiOptions+='-Iimgui'
imguiOptions+='-Iimgui/backends'
imguiOptions+='-F '${frameworksDir}
imguiOptions+='-framework'
imguiOptions+='OpenGL'
imguiOptions+='-framework'
imguiOptions+='GLUT'
imguiOptions+='-L/Users/aionfeehan/opt/anaconda3/lib/'

buildFiles=()
for file in ${baseFiles[@]}
do buildFiles+=${file}.cpp
done


cd ${workingDir}
rm ${logFile}
touch ${logFile}

objectFiles=()
for file in ${baseFiles[@]}
do
objectFile=${buildDir}${file}.o
if [ $logging -gt 0 ]; then
buildCommand="${compiler} ${file}.cpp -c -o ${objectFile} ${optionList[@]} &>> ${logFile}"
echo "Running ${buildCommand}" &>> ${logFile}
${compiler} ${file}.cpp -c -o ${objectFile} ${optionList[@]} &>> ${logFile} 
else
buildCommand="${compiler} ${file}.cpp -c -o ${objectFile} ${optionList[@]}"
echo "Running ${buildCommand}"
${compiler} ${file}.cpp -c -o ${objectFile} ${optionList[@]}
fi
objectFiles+=${objectFile}
done

mkdir ${buildDir}${imguiDir}
mkdir ${buildDir}${imguiDir}/backends
for file in ${imguiBaseFiles[@]}
do
objectFile=${buildDir}${imguiDir}${file}.o
if [ $logging -gt 0 ]; then
buildCommand="${compiler} ${imguiDir}${file}.cpp -c -o ${objectFile} ${optionList[@]} ${imguiOptions[@]} &>> ${logFile}"
echo "Running ${buildCommand}" &>> ${logFile}
${compiler} ${imguiDir}${file}.cpp -c -o ${objectFile} ${optionList[@]} ${imguiOptions[@]} &>> ${logFile} 
else
buildCommand="${compiler} ${imguiDir}${file}.cpp -c -o ${objectFile} ${optionList[@]} ${imguiOptions[@]}"
echo "Running ${buildCommand}"
${compiler} ${imguiDir}${file}.cpp -c -o ${objectFile} ${optionList[@]} ${imguiOptions[@]}
fi
objectFiles+=${objectFile}
done
# ${compiler} ${buildFiles[@]} -c -o ${libraryName} ${optionList[@]} &>> ${logFile}

# if [ $logging -gt 0 ]; then
# ${compiler} matrix_benchmarks.cpp -o ${buildDir}matrix_benchmarks.out ${objectFiles} -DIS_MAIN ${optionList[@]} &>> ${logFile}
# else
# ${compiler} matrix_benchmarks.cpp -o ${buildDir}matrix_benchmarks.out ${objectFiles} -DIS_MAIN ${optionList[@]}
# fi

if [ $logging -gt 0 ]; then
buildCommand="${compiler} tests.cpp -o ${buildDir}tests.out ${objectFiles} ${optionList[@]} -DIS_MAIN &>> ${logFile}"
echo ${buildCommand} $>> ${logFile}
${compiler} tests.cpp -o ${buildDir}tests.out ${objectFiles} ${optionList[@]} -DIS_MAIN &>> ${logFile}
else
buildCommand="${compiler} tests.cpp -o ${buildDir}tests.out ${objectFiles} ${optionList[@]} -DIS_MAIN"
echo ${buildCommand}
${compiler} tests.cpp -o ${buildDir}tests.out ${objectFiles} ${optionList[@]} -DIS_MAIN
fi

if [ $logging -gt 0 ]; then
buildCommand="${compiler} main.cpp -o ${buildDir}main ${objectFiles} ${optionList[@]} ${imguiOptions[@]} &>> ${logFile}"
echo ${buildCommand} $>> ${logFile}
${compiler} main.cpp -o ${buildDir}main ${objectFiles} ${optionList[@]} ${imguiOptions[@]} &>> ${logFile}
else
buildCommand="${compiler} main.cpp -o ${buildDir}main ${objectFiles} ${optionList[@]} ${imguiOptions[@]}"
echo ${buildCommand}
${compiler} main.cpp -o ${buildDir}main ${objectFiles} ${optionList[@]} ${imguiOptions[@]}
fi