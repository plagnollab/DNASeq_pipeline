pipeline=/cluster/project8/vyp/vincent/Software/pipeline/msample_calling/msample_calling.sh


genotype=no  #step1
recal=yes  #step2
annovar=yes ##step 3

bash $pipeline --currentUCLex January2015 --gVCFlist support/gVCF.tab  --genotype ${genotype} --recal ${recal} --annovar ${annovar}