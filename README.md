![banner](https://raw.githubusercontent.com/microgenomics/tutorials/master/img/microgenomics.png)

#analysisMethods
----------------

Analysis Methods Module is a final module from SEPA (Simulation, Execution, Parse and Analysis,  a multi pipe to test pathogen software detection), that analyze a cvs file provided by Parse Module (or your own cvs file in particular format), the analysis aviable are R<sup>2</sup>, RRMSE (relative error) and ROC Curve

##Requirements

* Bash version 4 (comes with Linux and MacOSX)
* Internet connection if your realdata file (maybe metasim\_conf_file.mprf or same format [abundance 'gi or ti' ginumber]) contain Genomes ID (gi)

##Usage
analysisMethods provide three methods for the csv file:

* R<sup>2</sup> (determination coefficient)
* RRMSE (relative root mean square error)
* ROC Curve (receiver operating characteristic)

you need to specify what analysis you want in the configuration file like this:

	# parseMethods configuration file (config.conf)#
	# add comments using the pound character
	
	ANALYSISTYPE=ROC_CURVE,R2,RMS_RE
	TOTALGENOMES=426

	#ANALYSISTYPE is a flag that specify the analysis that you want to do.
	#TOTALGENOMES is a flag that you need to specify if ROC CURVE analysis is specified. its mean the universe of genomes (data base) that you used in mapping step (all genomes)
	
in example, three methods are required, can be only RMS_RE too or whatever.

Next, make sure you have a csv file provided by Parse Module that have the following format:

	Superkingdom,Phylum,Class,Order,Family,Genus,Specie,Name,ti,filename1,filename2,filenameN
	Bacteria,Firmicutes,Bacilli,Bacillales,Staphylococcaceae,Staphylococcus,Staphylococcus epidermidis,Staphylococcus epidermidis RP62A,"176279",20286,19568,19784,19562
	Bacteria,Firmicutes,Bacilli,Bacillales,Bacillaceae,Bacillus,Bacillus anthracis,Bacillus anthracis str. CDC 684,"568206",86380,85306,85925,87788


	
Where the first eight column names are the complete lineage of the pathogen, ti is a NCBI taxonomy identifier, and filenameX is the number of reads that software mapped on corresponding ti and saved in file called filenameX.

Also you can provide your own csv file without lineage (or within) like this:
	
	ti,filename1,filename2,filenameN
	176279,20286,19568,19784,19562
	568206",86380,85306,85925,87788
	etc,etc,etc,etc,etc

Then you can execute the module typing the following command:
	
	bash analysisMethods.bash --workpath [files directory] --cfile [config file] --simdata [simulation table csv format (',' separated)] --realdata [mconf from metasim] --local

####Options explanation:
* --workpath path where directory tree or your files are (included requirement files)
* --cfile is the configuration file
* --local is if you don't use the other modules of SEPA
* --simdata is the csv file that contain the simulated data
* --realdata a file that contain the real values of reads distribution (.mprf of metasim)

## Example output

R<sup>2</sup> (R2.dat):
	
	Analysis	ID1000-sam	ID150-sam	ID300-sam	ID75-sam
	R2			0.5906		0.6408		0.6379		0.6317
	
RRMSE (rms.dat):
	
	Analysis ID1000-sam	ID150-sam	ID300-sam	ID75-sam
	RRMSE	0.766182	0.580528	0.583302	0.593402
	AVGRE	0.662465	0.468544	0.469835	0.480457
	
ROC Curve (ROC.dat):

	file		fpr		  tpr
	ID1000-sam	0.0311751 0.714286
	ID150-sam	0.0311751 0.8
	ID300-sam	0.0311751 0.8
	ID75-sam	0.0311751 0.8
