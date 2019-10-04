#!/bin/bash

echo "#!/bin/bash" >> run.sh
count=0
for Med in 195 200 295 300 495 500 995 1000 1995 2000 2495 2500 2995 3000
do
    for DM in 50 100 150 250 300 500 750 1000 1250 1500 2000
    do  
        targetDir=Vector_MonoTop_NLO_Mphi-${Med}_Mchi-${DM}_gSM-0p25_gDM-1p0_13TeV-madgraph

        if [ -e cards/MadGraph5Cards_2016/${targetDir} ]; then

            count=$((count+1))

            echo "nohup ./submit_cmsconnect_gridpack_generation.sh ${targetDir} cards/MadGraph5Cards_2016/${targetDir} > myLog_$count.debug 2>& 1 &" >> run.sh                                                            

        fi  
    done
done

