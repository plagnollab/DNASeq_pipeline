
py=/share/apps/python-2.7.3-static/bin/python

cd /cluster/project8/vyp/AdamLevine/annotations/barcode

##Download datasets and extract list of samples

#HGU133a (Human) v3
#AffyHGU133A.ann
wget -O HGU133aHumanV3.csv http://barcode.luhs.org/transcriptomes/abc-gpl96-formatted_v3.csv
$py scripts/work.py HGU133aHumanV3.csv > samplesHGU133aHumanV3.csv
wget -O HGU133aHumanV3.ann http://barcode.luhs.org/annotations/tabGPL96v3.csv

#HGU133plus2 (Human) tissues v3
#AffyHG133Plus2.ann
wget -O HGU133plus2HumanTissuesV3.csv http://barcode.luhs.org/transcriptomes/abc-tis-gpl570-formatted_v3.csv
$py scripts/work.py HGU133plus2HumanTissuesV3.csv > samplesHGU133plus2HumanTissuesV3.csv
wget -O HGU133plus2HumanTissuesV3.ann http://barcode.luhs.org/annotations/tabGPL570tissuev3.csv

#HGU133plus2 (Human) cells v3
#AffyHG133Plus2.ann
wget -O HGU133plus2HumanCellsV3.csv http://barcode.luhs.org/transcriptomes/abc-cell-gpl570-formatted_v3.csv
$py scripts/work.py HGU133plus2HumanCellsV3.csv > samplesHGU133plus2HumanCellsV3.csv
wget -O HGU133plus2HumanCellsV3.ann http://barcode.luhs.org/annotations/tabGPL570celltypev3.csv

#HGU133a2 (Human) v3
#AffyHGU133A_2.ann
wget -O HGU133a2HumanV3.csv http://barcode.luhs.org/transcriptomes/abc-gpl571-formatted_v3.csv
$py scripts/work.py HGU133a2HumanV3.csv > samplesHGU133a2HumanV3.csv
wget -O HGU133a2HumanV3.ann http://barcode.luhs.org/annotations/tabGPL571_v3.csv

#Human1.0ST (Human) tissues v3
#AffyHuGene1ST.ann
wget -O Human1STTissuesV3.csv http://barcode.luhs.org/transcriptomes/abc-tis-core-gpl6244-formatted.csv
$py scripts/work.py Human1STTissuesV3.csv > samplesHuman1STTissuesV3.csv
wget -O Human1STTissuesV3.ann http://barcode.luhs.org/annotations/tabGPL6244celltype.csv

#Human1.0ST (Human) cells v3
#AffyHuGene1ST.ann
wget -O Human1STCellsV3.csv http://barcode.luhs.org/transcriptomes/abc-cell-core-gpl6244-formatted.csv
$py scripts/work.py Human1STCellsV3.csv > samplesHuman1STCellsV3.csv
wget -O Human1STCellsV3.ann http://barcode.luhs.org/annotations/tabGPL6244celltype.csv

##Annotation files from BioMart
#AffyHG133Plus2.ann
#AffyHGU133A_2.ann
#AffyHGU133A.ann
#AffyHuGene1ST.ann

echo HGU133aHumanV3.csv > files.txt
echo HGU133plus2HumanTissuesV3.csv >> files.txt
echo HGU133plus2HumanCellsV3.csv >> files.txt
echo HGU133a2HumanV3.csv >> files.txt
echo Human1STTissuesV3.csv >> files.txt
echo Human1STCellsV3.csv >> files.txt


touch temp
rm temp
for i in `cat files.txt`;
do echo $i
awk -v source=${i} '{print source "\t" $0}' samples${i} >> temp
done
awk '$2!=""' temp > samples.txt
rm temp

#Manually create samples_annotated.txt

#Extract the bowel/immune expression information
${py} scripts/run.py samples_annotated.txt files.txt

#Prepare the ENSG to probe conversion
python scripts/to_ensg.py AffyHG133Plus2.ann 1 3 extracted_HGU133plus2HumanTissuesV3.csv > converted_HGU133plus2HumanTissuesV3.csv
python scripts/to_ensg.py AffyHG133Plus2.ann 1 3 extracted_HGU133plus2HumanCellsV3.csv > converted_HGU133plus2HumanCellsV3.csv
python scripts/to_ensg.py AffyHGU133A_2.ann 1 9 extracted_HGU133a2HumanV3.csv > converted_HGU133a2HumanV3.csv
python scripts/to_ensg.py AffyHGU133A.ann 1 9 extracted_HGU133aHumanV3.csv > converted_HGU133aHumanV3.csv
python scripts/to_ensg.py AffyHuGene1ST.ann 1 9 extracted_Human1STTissuesV3.csv > converted_Human1STTissuesV3.csv
python scripts/to_ensg.py AffyHuGene1ST.ann 1 9 extracted_Human1STCellsV3.csv > converted_Human1STCellsV3.csv

cat converted* | grep -v missing | sort > all.csv
python scripts/amalgamate.py  all.csv  > processed.txt

#Manually retrieve biomart.ann
# python scripts/gene.py | sort | uniq | sort -k1n -k2n > biomart.map
# python scripts/final.py  processed.txt biomart.map   > barcode.bed




