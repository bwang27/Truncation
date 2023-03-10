source ../../setGlobalVars.sh
echo $pipeline_root


downloadDir=$pipeline_root/ReferenceDB/geneOntology/01.download/
mkdir -p $downloadDir
cd $downloadDir

wget -O go.obo http://purl.obolibrary.org/obo/go.obo 
wget ftp://ftp.ncbi.nlm.nih.gov/gene/DATA/gene2go.gz


#need to install go-perl



##create GO to GO table
cd $pipeline_root/ReferenceDB/geneOntology/
perl 02.getGO_assn.all.pl

#create GO to gene ID table
perl 03.getGO2Gene.pl -s hs -t 9606  &
perl 03.getGO2Gene.pl -s mm -t 10090 &
perl 03.getGO2Gene.pl -s rn -t 10116 &
wait

