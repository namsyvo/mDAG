import os
import sys
import re
import gc

from datetime import datetime

import sqlite3
#import MySQLdb

from Methods import *

#Reading the arguments
Argument = []
Argument = sys.argv[1:]
app_folder = Argument[0]
result_id = Argument[1]


def Executer():

    """
    Reading information
    This code reads dataset information, gene expression data file and other information
    """
    print "Reading information for result id: " + result_id

    db_path = os.path.join(app_folder, 'databases', 'storage.db')
    try:
	#Connecting to database and reading information
        #con = MySQLdb.connect(host = "localhost", user = "mdat", passwd = "mdat", db = "mdat")
	con = sqlite3.connect(db_path)
	con.row_factory = sqlite3.Row
	c = con.cursor()
	c.execute('SELECT * FROM results WHERE id = "' + result_id + '"')
	result = c.fetchone()
        dataset_id = str(result['dataset_id'])
        c.execute('SELECT * FROM datasets WHERE id = "' + dataset_id + '"')
        dataset = c.fetchone()
	c.close()
	con.close()
    except:
	print "Error for reading input parameters in Executer:", sys.exc_info()[0]
	return 0

    FDRrate = result['sign_level'] # False Discovery Rate for analysis
    permutation = result['perm_number'] # permutation number for Monte Carlo method

    nameoftreat = []
    numofrep = []
    rep_number = result['rep_number']
    temparr = rep_number.split("\n")
    for i in range(len(temparr)):
        temp = temparr[i].strip()
        temp = temp.split(" ")
        nameoftreat.append(temp[0])
        numofrep.append(int(temp[1]))

    numoftreat = (len(numofrep)) # Finding the number of treatments involved

    Indexarray = range(sum(numofrep))

    #Reading dataset file
    Header = [] #Header
    probeset = [] #Probe ID
    genesymbol = [] #Gene Symbol
    unigeneid = [] #UniGene ID
    reppubid = [] #Representative Public ID
    otherinfo = [] #Chromosomal Location
    Expressionmatrix = [] #Expression data

    #Taking path to dataset file
    file_path = os.path.join(app_folder, 'uploads', dataset['file_name'])
    try:
	d = file(file_path).read()
    except IOError:                     
	print "Error for reading input file in Executer:", sys.exc_info()[0]
	return 0

    rowspre = d.split('\n')
    rows = []
    for l in rowspre:
	if l.strip():
	    rows.append(l)
    NRows = len(rows) # number of rows in the dataset file
    Header = rows[0].split('\t')
    for i in range(1, NRows):
	rowlist = []
	rowlist = rows[i].split('\t')
	probeset.append(rowlist[0])
	genesymbol.append(rowlist[1])
	unigeneid.append(rowlist[2])
	reppubid.append(rowlist[3])
	otherinfo.append(rowlist[4])
	Expressionmatrix.append(map(float, rowlist[5:]))   # Storing Expression data as 2-D array


    print "Finish reading information for result id: " + result_id


    """
    Analyzing data
    This code uses other modules to perform Kruskal Wallis Test, False discovery rate, calculate Pvalue using MonteCarlo Permutation, perform PostHocComparison and generates patterns.
    It also uses Graphviz software to generate .svg and .png files to show clusters in the form of graphs.
    """
    print "Start selecting significant genes for result id: " + result_id


    """
    Calculating Hstatistics for the original or unpermuted data
    """
    Hstatistics = [] 
    for j in range(0,len(Expressionmatrix)):  # Finding and storing Hstatistics, KWPvalue 
        Hstat = 0; 
	Hstat = KruskalWallisTest(Expressionmatrix[j], numofrep, Indexarray) 
	Hstatistics.append(Hstat)


    """
    Calculating Hstatistics for Permuted data 
    """
    MonteCarloCounter = [0]*(len(Expressionmatrix)) 

    from datetime import timedelta
    start = datetime.now()
	
    for p in range(permutation):

        print "Calculating Hstats for permutation " + str(p) + " of result id " + str(result_id)

	Indexarray = Permutation(Indexarray)

	for k in range(0,len(Expressionmatrix)): # Finding and storing Hstatistics, KWPvalue 
	    Hstatp = 0
	    Hstatp = KruskalWallisTest(Expressionmatrix[k], numofrep, Indexarray)
		
	    if Hstatp >= Hstatistics[k]:
		MonteCarloCounter[k] = MonteCarloCounter[k] + 1

	"""
	Code for calculating approximate running time for Monte Carlo method using the time for 5 loops
	"""
	if p == 50:
	    try:
		#Connect to database again
		#con = MySQLdb.connect(host = "localhost", user = "mdat", passwd = "mdat", db = "mdat")
		con = sqlite3.connect(db_path)
		c = con.cursor()
		unit_duration = datetime.now() - start
		unit_seconds = unit_duration.seconds
		total_seconds = int (unit_seconds * permutation / 50)
		c.execute('UPDATE results SET description = ? WHERE id = ?', (total_seconds, result_id))
		con.commit()
		c.close()
		con.close()
	    except:
		print "Error for set running time in Executer:", sys.exc_info()[0]
                return 0

    """
    Calculating MonteCarloPvalue 
    """
    MonteCarlopvalue = []
    for counter in MonteCarloCounter:
        MonteCarlopval = 0
	MonteCarlopval = float(1 + counter)/float(permutation + 1)
	MonteCarlopvalue.append(MonteCarlopval)


    """
    Peforming FDR check for significant differentially expressed genes
    """

    Significantpvalues = []
    Significantpvalues = FDRAnalysis(MonteCarlopvalue, FDRrate)


    """
    If there is no gene significantly differentially expressed, return
    """

    if Significantpvalues == []:

	print "No gene is significantly differentially expressed for result id: " + result_id

	try:
            #Connect to database again
	    #con = MySQLdb.connect(host = "localhost", user = "mdat", passwd = "mdat", db = "mdat")
	    con = sqlite3.connect(db_path)
	    c = con.cursor()

	    timestamp = datetime.now()

	    status_update_sql='UPDATE results SET finish_time = ?, status = ? WHERE id = ?'

	    c.execute(status_update_sql, (timestamp, "No gene is significantly differentially expressed", result_id))
	    con.commit()

	    c.close()
	    con.close()

	except:
	    print "Error for updating database in Executer:", sys.exc_info()[0]

	return 0

    else:
        print "Finish selecting significant genes for result id: " + result_id



    """
    Peforming PostHoc testing and printing the Ternary Pattern for significant differentially expressed gene
    """
    print "Start peforming post hoc comparisons for result id: " + result_id

    """
    Taking significant gene indexes and their pvalues, fold-change values
    """
    Uniquepvalues = []
    Uniquepvalues = unique(Significantpvalues)

    SignificantGenesindex = []
    Pvaluedict = {}
    FoldChangedict = {}

    for pvalue in Uniquepvalues:
        genearray = []
	genearray  = get_index_new(MonteCarlopvalue, pvalue)
	for genekey in genearray:
	    SignificantGenesindex.append(genekey)
	    Pvaluedict[genekey] = MonteCarlopvalue[genekey]

	    exp_mean = []
	    left = right = 0
	    count = 0
	    for m in numofrep:
	        right += m
		exp_mean.append(math.fsum(Expressionmatrix[genekey][left:right])/m)
		if count > 0:
		    if genekey not in FoldChangedict:
		        FoldChangedict[genekey] = []
		    FoldChangedict[genekey].append(exp_mean[count]/exp_mean[0])
		left = right
		count += 1


    """
    Taking patterns and gene indexes in each pattern
    """
    CombinationList = [] 
    for y in combinations(range(1,(numoftreat+1)),2): CombinationList.append(tuple(y))


    Patterndict = {}
    IndexPatterndict = {}

    for l in SignificantGenesindex:
		
        print "Performing post hoc comparison for gene index " + str(l) + " of result id " + str(result_id)

        Pattern = ''
	Treatmentdic = {}
	left = right = 0
	count = 0

	for m in numofrep:
	    right = right + m
	    count = count+1
	    Treatmentdic[count] = Expressionmatrix[l][left:right]
	    left = right
		
	for (a,b) in CombinationList:
	    Pvalue,PvalueC = RanksumComparison(Treatmentdic[a],Treatmentdic[b])

	    if Pvalue < 0.05:
	        Pattern = Pattern+str(2)
	    elif PvalueC < 0.05:
		Pattern = Pattern+str(0)
	    else:
		Pattern = Pattern+str(1)

	if Pattern in Patterndict:
	    Patterndict[Pattern].append(probeset[l])
	    IndexPatterndict[Pattern].append(l)
	else:
	    Patterndict[Pattern] = []
	    IndexPatterndict[Pattern] = []

	    Patterndict[Pattern].append(probeset[l])
	    IndexPatterndict[Pattern].append(l)


    Binarydict = {}

    for keys in Patterndict.keys():
        Binarydict[int(keys)] = keys


    print "Finish peforming post hoc analysis for result id: " + result_id


    """
    Creating results and store them in database
    """

    print "Start creating results for result id: " + result_id

    #Generating directories where result files will be stored
    result_dir_path = os.path.join(app_folder, 'static', 'results', str(result_id))
    os.mkdir(result_dir_path)

    first_dir_path = os.path.join(result_dir_path, 'Firstlevel')
    os.mkdir(first_dir_path)

    second_dir_path = os.path.join(first_dir_path, 'Secondlevel')
    os.mkdir(second_dir_path)

    """
    Creating overview graph
    """
    Overviewclusters(SignificantGenesindex, Patterndict, numoftreat, numofrep, nameoftreat, FDRrate, result_dir_path)

    """
    Creating all graphical patterns for the first level cluster
    """
    for pattern in Patterndict.keys():
	FirstLevelClusters(pattern, len(Patterndict[pattern]), numoftreat, nameoftreat, first_dir_path)


    """
    Creating all graphical patterns for the second level cluster
    """
    Resultdict, Symboldict = SecondLevelClusterDictionary(nameoftreat, numoftreat, Patterndict)

    Sumdict = {}
    Totalsum = 0

    for a in Resultdict.keys():
	Sum = 0
	for b in Resultdict[a]:
	    Sum = Sum + len(Patterndict[b])
	    Totalsum = Totalsum + len(Patterndict[b])
	
	Sumdict[a] = Sum

		
    Contractratiodict = {}

    for second_cluster in Resultdict.keys():

	allpatterns = []
	Allpattern = {}
	allpatterns = Resultdict[second_cluster]
	contractiblecount = 0.0;

	for eachpattern in allpatterns:
	    Allpattern[eachpattern] = len(Patterndict[eachpattern])
	    if Contractibility(numoftreat,str(eachpattern)):
	        contractiblecount = contractiblecount + len(Patterndict[eachpattern])

	Contractratio = float(contractiblecount)/float(Sumdict[second_cluster])
	Contractratiodict[second_cluster] = round(Contractratio,3)

	SecondLevelClusters(second_cluster, Allpattern, Contractratio, Resultdict, Symboldict, Sumdict, Patterndict, numoftreat, nameoftreat, second_dir_path)


    """
    Creates gene files for all significant genes and genes in each pattern
    """
    SignificantGeneListFiles(SignificantGenesindex, Pvaluedict, FoldChangedict, Header, probeset, genesymbol, unigeneid, reppubid, otherinfo, numoftreat, nameoftreat, numofrep, Expressionmatrix, result_dir_path)
    
    #geneIndex = range(NRows-1)
    #AllGeneListFiles(geneIndex, Pvaluedict, FoldChangedict, Header, probeset, genesymbol, unigeneid, reppubid, otherinfo, numoftreat, nameoftreat, numofrep, Expressionmatrix, result_dir_path)

    for pattern in IndexPatterndict.keys():
	PatternGeneListFiles(pattern, IndexPatterndict, Pvaluedict, FoldChangedict, Header, probeset, genesymbol, unigeneid, reppubid, otherinfo, numoftreat, nameoftreat, numofrep, Expressionmatrix, first_dir_path)

    for symbol in Symboldict.keys():
        MetaPatternGeneListFiles(symbol, Symboldict, Resultdict, IndexPatterndict, Pvaluedict, FoldChangedict, Header, probeset, genesymbol, unigeneid, reppubid, otherinfo, numoftreat, nameoftreat, numofrep, Expressionmatrix, second_dir_path)


    """
    Storing all results in database
    """

    import cPickle
    import base64

    Binarydict_string = base64.b64encode(cPickle.dumps(Binarydict, cPickle.HIGHEST_PROTOCOL))
    Patterndict_string = base64.b64encode(cPickle.dumps(Patterndict, cPickle.HIGHEST_PROTOCOL))
    IndexPatterndict_string = base64.b64encode(cPickle.dumps(IndexPatterndict, cPickle.HIGHEST_PROTOCOL))
    MonteCarlopvalue_string = base64.b64encode(cPickle.dumps(MonteCarlopvalue, cPickle.HIGHEST_PROTOCOL))

    Pvaluedict_string = base64.b64encode(cPickle.dumps(Pvaluedict, cPickle.HIGHEST_PROTOCOL))
    FoldChangedict_string = base64.b64encode(cPickle.dumps(FoldChangedict, cPickle.HIGHEST_PROTOCOL))

    Resultdict_string = base64.b64encode(cPickle.dumps(Resultdict, cPickle.HIGHEST_PROTOCOL))
    Symboldict_string = base64.b64encode(cPickle.dumps(Symboldict, cPickle.HIGHEST_PROTOCOL))
    Sumdict_string = base64.b64encode(cPickle.dumps(Sumdict, cPickle.HIGHEST_PROTOCOL))
    Contractratiodict_string = base64.b64encode(cPickle.dumps(Contractratiodict, cPickle.HIGHEST_PROTOCOL))

    try:
	#Connect to database again to store results
	#con = MySQLdb.connect(host = "localhost", user = "mdat", passwd = "mdat", db = "mdat")
	con = sqlite3.connect(db_path)
	con.text_factory = str
	c = con.cursor()

	#Set status to Completed
	timestamp = datetime.now()

	status_update_sql='UPDATE results SET finish_time = ?, status = ?, binary_dict = ?, pattern_dict = ?, index_pattern_dict = ?, monte_carlo_pvalue = ?, result_dict = ?, symbol_dict = ?, sum_dict = ?, contract_dict = ?, pvalue_dict = ?, fold_change_dict = ? WHERE id = ?'

	c.execute(status_update_sql, (timestamp, "Completed", Binarydict_string, Patterndict_string, IndexPatterndict_string, MonteCarlopvalue_string, Resultdict_string, Symboldict_string, Sumdict_string, Contractratiodict_string, Pvaluedict_string, FoldChangedict_string, result_id))

	con.commit()

	c.close()
	con.close()

    except:
	print "Error for updating database in Executer:", sys.exc_info()[0]
	return 0


    print "Finish creating results for result id: " + result_id
    return 1



"""
Main program
"""
if __name__ == "__main__":

    print '############################################'
    print 'Start analyzing for result id: ' + result_id

    r = Executer()
    
    if r == 1:
        print 'Finish analyzing for result id: ' + result_id
        print '#############################################'
    else:
        print 'Fail to analyze for result id: ' + result_id
        print '############################################'
        

    """
    Call Scheduler again after finish this analysis
    """
    scheduler_path = os.path.join(app_folder, 'modules', 'Scheduler.py')
    cmd = 'python ' + scheduler_path + ' ' + app_folder + ' &'
    os.system(cmd)

    gc.collect()
    sys.exit()
