#!/bin/bash

set -u
set -x

java17=/share/apps/jdk1.7.0_45/jre/bin/java
inFile=J6_sorted_unique.bam
outFile=J6_sorted_unique-out.bam
software=/cluster/project8/vyp/vincent/Software
AddOrReplaceReadGroups=${software}/picard-tools-1.100/AddOrReplaceReadGroups.jar
tempFolder=${SCRATCH2}/vincent/temp/novoalign
code=Hardcastle_Czech_J6

# INPUT (String)  Input file (bam or sam or a GA4GH url). Required.
# OUTPUT (File)   Output file (bam or sam). Required.
# SORT_ORDER (SortOrder)  Optional sort order to output in. If not supplied OUTPUT is in the same order as INPUT. Default value: null. Possible values: {unsorted, queryname, coordinate, duplicate}

# RGID (String)   Read Group ID Default value: 1. This option can be set to 'null' to clear the default value.
# RGLB (String)   Read Group Library Required.
# RGPL (String)   Read Group platform (e.g. illumina, solid) Required.
# RGPU (String)   Read Group platform unit (eg. run barcode) Required.
# RGSM (String)   Read Group sample name Required.
# RGCN (String)   Read Group sequencing center name Default value: null.
# RGDS (String)   Read Group description Default value: null.
# RGDT (Iso8601Date)  Read Group run date Default value: null.
# RGPI (Integer)  Read Group predicted insert size Default value: null.
# RGPG (String)   Read Group program group Default value: null.
# RGPM (String)   Read Group platform model Default value: null.

echo java -Djava.io.tmpdir=${tempFolder} -Xmx4g -jar $AddOrReplaceReadGroups I=$inFile O=$outFile RGLB=$code RGPL=illumina RGPU=run RGSM=$code RGID=$code


#VALIDATION_STRINGENCY=LENIENT 


