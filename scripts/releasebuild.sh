#!/bin/zsh
logging=-1
compiler='/usr/bin/clang++'
workingDir='/Users/aionfeehan/Github/analytics/'
buildDir=${workingDir}release/
logFile='releasebuild.log'

optionList=('-fcolor-diagnostics')
optionList+='-fansi-escape-codes'
optionList+='-Ofast'
optionList+='-std=gnu++20'

baseFiles=('linalg')
baseFiles+='matrix'
baseFiles+='basic_math'

buildFiles=()
for file in ${baseFiles[@]}
do buildFiles+=${file}.cpp
done

objectFiles=()
for file in ${baseFiles[@]}
do objectFiles+=${buildDir}${file}.o
done

cd ${workingDir}
rm ${logFile}
touch ${logFile}

for file in ${baseFiles[@]}
do
if [ $logging -gt 0 ]; then
${compiler} ${file}.cpp -c -o ${buildDir}${file}.o ${optionList[@]} &>> ${logFile} 
else
${compiler} ${file}.cpp -c -o ${buildDir}${file}.o ${optionList[@]}
fi
done

# if [ $logging -gt 0 ]; then
# ${compiler} matrix_benchmarks.cpp -o ${buildDir}matrix_benchmarks.out ${objectFiles} -DIS_MAIN ${optionList[@]} &>> ${logFile}
# else
# ${compiler} matrix_benchmarks.cpp -o ${buildDir}matrix_benchmarks.out ${objectFiles} -DIS_MAIN ${optionList[@]}
# fi

if [ $logging -gt 0 ]; then
${compiler} tests.cpp -o ${buildDir}tests.out ${objectFiles} ${optionList[@]} -DIS_MAIN &>> ${logFile}
else
${compiler} tests.cpp -o ${buildDir}tests.out ${objectFiles} ${optionList[@]} -DIS_MAIN
fi
