# This script runs the stability analysis for Hether and Hohenlohe (2013)
# Created by TDH
# Aug. 13th, 2013

# Be sure to cd into Stability directory prior to running this script

#echo 'starting simulations'

R CMD BATCH ParameterSpace.R
cd ./Runs/

INPUT_FILE=../ParamSpace.txt

while read line
    do
    set 'echo $line'
    F1=$(echo $line | awk '{print substr($1, index($7, $9))}')
    F2=$(echo $line | awk '{print substr($2, index($7, $9))}')
    F3=$(echo $line | awk '{print substr($3, index($7, $9))}')
    F4=$(echo $line | awk '{print substr($4, index($7, $9))}')
    F5=$(echo $line | awk '{print substr($5, index($7, $9))}')
    F6=$(echo $line | awk '{print substr($6, index($7, $9))}')
    F7=$(echo $line | awk '{print substr($7, index($7, $9))}')
    F8=$(echo $line | awk '{print substr($8, index($8, $9))}')   
    F9=$(echo $line | awk '{print substr($9, index($9, $9))}')
    F10=$(echo $line | awk '{print substr($10, index($10, $10))}')        
    qsub -v MU=$F1,REGMU=$F2,OM11=$F3,OM12=$F4,N=$F5,OPTIMUM=$F6,SWITCHFREQ=$F7,MOTIF=$F8,REP=$F9,M=$F10 ../Rsub.pbs
# echo "qsub -v MU="$F1",REGMU="$F2",OM11="$F3",OM12="$F4",N="$F5",OPTIMUM="$F6",SWITCHFREQ="$F7",MOTIF="$F8",REP="$F9",M="$F10" ../Rsub.pbs" >> tmp.txt

done < $INPUT_FILE

# Email when finished:












