pipeline=/cluster/project8/vyp/vincent/Software/pipeline/msample_calling/msample_calling.sh


genotype=yes  #step1
recal=no  #step2
annovar=no ##step 3
convertToR=no  ##step 4

bash $pipeline --currentUCLex January2015b --gVCFlist support/gVCF.tab  --genotype ${genotype} --recal ${recal} --annovar ${annovar} --convertToR ${convertToR}

