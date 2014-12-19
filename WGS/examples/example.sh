# bash script

reference=1kg
projectID=Project_Example

# generate support file
bash makesupport.sh ../data/ .fastq.gz > support.txt

# alignment
bash ../WGS_pipeline.sh --mode align\
                        --supportFrame support.txt
                        --reference ${reference}
                        --tparam 320
                        --inputFormat STDFQ 
                        --projectID ${projectID}
                        


