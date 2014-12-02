

The two main scripts are:
- WGS_submission_script.sh
- WGS_pipeline.sh

*WGS_submission_script* defines the command line arguments to call
*WGS_pipeline*
with.

`
sh ${pipeline} --supportFrame ${supportFrame} --reference ${reference} --align ${align}  --tparam 320  --inputFormat STDFQ  --extraID $extraID --makeVCF ${makeVCF} --makegVCF ${makegVCF}  --projectID ${projectID}
`


The pipeline submission script assumes that you have written a "support" data frame.
The format of the first row is:
> code f1 f2 
where f1 and f2 are the paired end fastq files and code is an identifier for the exome of interest.

The running options are:
- makegVCF
- makeVCF
- align

Depending on the arguments *WGS_pipeline* will generate the following bash
scripts:

- cluster/submission/align.sh
- cluster/submission/align_table.sh

or

- cluster/submission/makeVCF.sh
- cluster/submission/makeVCF_table.sh


The default output folder is:
> aligned

A very useful parameter to set is "extraID" to know the batch.


