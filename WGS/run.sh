set -x
set -u

DIR=$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )
script=${DIR}/WGS_pipeline.sh 

function align() {
    bash $script --mode align --supportFrame ${supportFrame} --reference ${reference} --aligner-tparam 320 --inputFormat STDFQ  --extraID ${extraID} --projectID ${projectID}
}

function call() {
    bash $script --mode gvcf --supportFrame ${supportFrame} --reference ${reference}  --projectID ${projectID}
}


source $1 
$2




