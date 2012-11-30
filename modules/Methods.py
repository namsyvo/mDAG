import os
import gc
import math


"""
Module takes list containing duplicated elements and returns list with unique elements
"""
def unique(List):
    checked = []
    for z in List:
        if z not in checked:
            checked.append(z)
    return checked


"""
Module returns the index of a particular element in list
"""
def get_index_new(array,item):
    count = 0
    List = []
    for i in array:
	if i == item:
            List.append(count)
	count = count + 1
    return List


"""
Module rounds off the decimal values to 2 places after decimal. For eg. it will round off 10.13240494 to 10.13
"""
def roundlist(List, num):
    RList = []
    for i in List:
        j = str(round(i, num))
        RList.append(j)
    
    return RList


"""
Module returns the Mean of numbers in a list
"""
def Mean(List):
    return (math.fsum(List))/len(List)


"""
Module returns the sample standard deviation of numbers in a list
"""
def SampStdDev(List):
    mean = Mean(List)
    sd = 0.0
    for i in List:
	sd += (float(i)-mean)**2
    return math.sqrt(sd/(len(List)-1))


"""
Module sorts the dictionary by value
"""
def sortbyvalue(d, reverse = False):
    return sorted(d.iteritems(), key = lambda x: x[1] , reverse = reverse)


"""
Module generates the complementary pattern for a pattern
"""
def comppattern(pattern):
    comppattern = ''
    for i in pattern:
        if i == '0':
            comppattern = comppattern + '2'
        elif i == '2':
            comppattern = comppattern + '0'
        elif i == '1':
            comppattern = comppattern + '1'
    return comppattern


"""
Module creates all possible combinations of the elements in a list.
"""
def combinations(items, n):
    if n == 0: yield []
    else:
        for i in range(len(items)):
            for cc in combinations(items[i + 1:], n - 1):
                yield [items[i]] + cc


#####################################################################################################
"""
Module implements False Discovery rate (FDR) procedure developed by Yoav Benjamini and Yosef Hochberg.
Input: List containing the P-values , significance level(Alpha) 
Output: List containing P-values qualifying as significant P-values after FDR correction
Reference: Controlling the False Discovery Rate: A Practical and Powerful Approach to Multiple Testing by Yoav Benjamini and Yosef Hochberg.
"""
def FDRAnalysis(Pvaluelist,signlevel):
    Pvaluelistsort = Pvaluelist[:]
    Pvaluelistsort.sort(lambda a,b: cmp(float(a), float(b))) 

    for i in range(len(Pvaluelistsort) - 1, -1, -1):
	if Pvaluelistsort[i] <= ((i*signlevel)/(len(Pvaluelistsort))):
            return Pvaluelistsort[:i]
    return []


"""
Module takes list as an argument and shuffles the elements using random number generator and returns Permuted list.
Input: List containing elements to be shuffled
Output: List containing shuffled elements from the original array
"""
import random

def Permutation(List):
    for i in range(len(List)):
	j = random.randint(0, i)
	List[i], List[j] = List[j], List[i]
    return List


""" 
Module SortnRank takes a list and returns sorted list and ranked list containing ranks of elements 
preserving order of that of original list.
Input: List containing data
Output: Sorted List , Ranked List 
"""
import math

def SortnRank(Array):
	
    Arraysort = Array[:]
    Arraysort.sort(lambda a,b: cmp(float(a), float(b))) 
    Arraysortc = Arraysort[:] 
    RankedArray = [0]*len(Arraysortc)
	
    for m in Arraysortc:
	RankedArray[Arraysortc.index(m)] = Array.index(m)
	Arraysortc[Arraysortc.index(m)] = -1
	Array[Array.index(m)] = -1
	
    return Arraysort,RankedArray


"""
Module takes List and resolves ties in ranks and substitutes them with average value.
This module uses module SortnRank for sorting and ranking.
Input: List containing data.
Output: Ranked list with ties resolved.
"""

def RanknTieresolver(inputlist):
    
    n = len(inputlist)
    sortlist, ranked = SortnRank(inputlist)
    sumranks = 0 
    dupcount = 0
    
    newlist = [0]*n
    for i in range(n):
        sumranks = sumranks + i
        dupcount = dupcount + 1
        if i == n - 1 or sortlist[i] != sortlist[i + 1]:
            averank = sumranks / float(dupcount) + 1
            for j in range(i - dupcount + 1, i + 1):
                newlist[ranked[j]] = averank
            sumranks = 0
            dupcount = 0
    return newlist


"""
Module calculates Hstatistics for a given data using KruskalWallis Procedure. 
Input: List containing Gene Expression data ,List containing indexes of permuted array,List containing number of replicates for each treatment.
Output: Hstatistics
"""
def KruskalWallisTest(GeneExpression, NumofReplicates, Indexarray):
        
    Listgene = GeneExpression[:]  	
    RankedList = []
    Totalelements = len(Listgene)
    RankedList = RanknTieresolver(Listgene)

    TotalSumSqDev = 0
    start = 0
    end = 0
	
    for m in NumofReplicates:
        GroupList = []
        SumSqDev = 0
		
        end = end + m
        for n in range(start,end):
            GroupList.append(RankedList[Indexarray.index(n)])
            start = end
        Ranksum = sum(GroupList)
	SumSqDev = (math.pow(Ranksum, 2))/m
	TotalSumSqDev = TotalSumSqDev + SumSqDev
			
    Hstatistics = float(12.0/(Totalelements*(Totalelements + 1))*TotalSumSqDev ) - float(3*(Totalelements + 1))
    return float(Hstatistics)


"""
Module prepares input to be sent to RecursionPvalue program. This function is based on Wilcoxon Rank Test procedure but uses Recursive method for obtaining P-values.
Input: Two Lists to be compared (List1,List2)
Output: Pvalue for two cases Pvalue1,Pvalue2
"""

def RanksumComparison(List1, List2):

    CommonList = []
    RankedList = []
	
    Numsample1 = len(List1)
    Numsample2 = len(List2)
    CommonList = List1 + List2
    RankedList = RanknTieresolver(CommonList)
    RanksumList1 = sum(RankedList[0:Numsample1])

    """
    This part implements Recursive equation given by P(n,m,r) = n/(n+m)*P(n-1,m,r-n-m)+ m/(n+m)*P(n,m-1,r)
    Memoization technique
    Base cases: P(1,0,k<=0) = 0, P(1,0,k>0) = 1, P(0,1,k<0) = 0, P(0,1,k>=0) = 1
    Argument n,m are positive integer types and r can be integer or float.
    Input: Elements in Sample 1, Elements in Sample 2, Sum of ranks of elements in Sample 1
    Reference: M.S. Ross Simulation 3rd edition, 2002
    """

    #Initialize P (memory array)
    max_size = 2*int(RanksumList1)
    P = [[[-1.0 for i in range(max_size + 1)] for j in range(Numsample2 + 1)] for k in range(Numsample1 + 1)]

    #Take index from r
    def idx(r):
	return int(r) + int(RanksumList1)

    #Recursive function
    def getPvalues(n, m, r):

	if n == 0.0:
            if m == 1.0:
		if r < 0.0:
                    return 0.0
		else:
                    return 1.0

	if n == 1.0:
            if m == 0.0:
		if r > 0.0:
                    return 1.0
		else:
                    return 0.0

	if n == 0.0:
            return getPvalues(n, m - 1, r)
	if m == 0.0:
            return getPvalues(n - 1, m, r - n - m)

	if P[n][m][idx(r)] != -1.0:
            return P[n][m][idx(r)]
	else:
            P[n][m][idx(r)] = (float(n)/float(n+m))*(getPvalues(n - 1, m, r - n - m)) +  (float(m)/float(n + m))*(getPvalues(n, m - 1, r))
            return P[n][m][idx(r)]

    Pvalue = getPvalues(Numsample1, Numsample2, RanksumList1)

    #Re-initialize P
    P = [[[-1.0 for i in range(max_size + 1)] for j in range(Numsample2 + 1)] for k in range(Numsample1 + 1)]
    PvalueC = 1 - (getPvalues(Numsample1, Numsample2, (RanksumList1 - 1)))

    return Pvalue, PvalueC


"""
Module creates the Matrix datatype
"""
import sys

class Matrix(object):
    def __init__(self, cols, rows):
        self.cols = cols
        self.rows = rows
        self.matrix = []

			
        for i in range(rows):
            each_row = []
            for j in range(cols):
                each_row.append(0)
            self.matrix.append(each_row)

    def setitem(self, col, row, v):
        self.matrix[col - 1][row - 1] = v
 
    def getitem(self, col, row):
        return self.matrix[col - 1][row - 1]

    def __repr__(self):
        outStr = ""
        for i in range(self.rows):
            outStr += 'Treat %s %s\n' % (i + 1, self.matrix[i])
        return outStr


"""
Module creates the Adjacency Matrix to represent graph.
It works in colaboration of other functions Adjacency Matrix, Contractible to find whether a Graph is Contractible or not.
Input: number of treatment, Pattern(string)
Output: True if Graph is contractible and vice-versa.
"""
def AdjacencyMatrix(numoftreat,Pattern):
	
    A = Matrix(numoftreat,numoftreat)
	
    CombList = []
    for comb in combinations(range(1, (numoftreat + 1)), 2): CombList.append(tuple(comb))
    count = 0
    for m,n in CombList:

        if Pattern[count] == '1':
            A.setitem(n, m, 0)
        elif Pattern[count] == '2':
            A.setitem(n, m, 1)
        elif Pattern[count] == '0':
            A.setitem(m, n, 1)
        count = count + 1

    return A


"""
Module checks the contractibility of the graph.
"""
def Contractibility(numoftreat, Pattern):
			
    CombList = []
	
    for comb in combinations(range(1, (numoftreat + 1)), 2): CombList.append(tuple(comb))
	
    A = AdjacencyMatrix(numoftreat, Pattern)

    # checking whether graph consists of all given vertices (Are there any vertices which are not part of graph )
    for h in range(1, (numoftreat + 1)):
	for o in range(1, (numoftreat + 1)):
            status = 0
            if A.getitem(h, o) != 0 or A.getitem(o, h) != 0:
		status = 1
		break

	if status == 0:
            return False

    # checking whether graph is a complete graph or not
    numberofedges = 0
    for m in range(1, (numoftreat + 1)):
	for n in range(1, (numoftreat + 1)):
            if A.getitem(m, n) == 1:
                numberofedges = numberofedges + 1

    maxedges = 0
    maxedges = (numoftreat*(numoftreat - 1))/2

    if numberofedges == maxedges:
	return True

    #Taking vertices that don't share any edge among them
    UnconnectedList = []
    for i, j in CombList:
	if A.getitem(i, j) == 0 and A.getitem(j, i) == 0:
            UnconnectedList.append(tuple([i, j]))	
	
    #Checking if all pair of unconnected vertices are pair of equivalent vertices, the graph is contractible iff the answer is Yes
    nonConFlag = 0
    for a, b in UnconnectedList:
	conFlag = 1
	for l in range(1,numoftreat + 1):
            if A.getitem(l,a) == A.getitem(l, b) and A.getitem(a, l) == A.getitem(b, l):
		pass
            else: 
		conFlag = 0
		break
	if conFlag == 0:
            nonConFlag = 1
            break

    if nonConFlag == 1:
	return False
    else:
	return True


"""
Module for topological sorting for contractive graph
"""
def ContractiveTopologicalSorting(numoftreat, nameoftreat, Pattern):
			
    # Checking whether graph connected or not (consists of all given vertices: Are there any vertices which are not part of graph?)
    A = AdjacencyMatrix(numoftreat, Pattern)
    for h in range(1,(numoftreat + 1)):
	status = 0
	for o in range(1,(numoftreat + 1)):
            if A.getitem(h,o) != 0 or A.getitem(o,h) != 0:
		status = 1
		break
	# If not connected
	if status == 0:
            return 'Not contractible'

    # If connected,
    # Constructing vertices list
    VList = []
    for i in range(len(nameoftreat)):
	VList.append(nameoftreat[i])

    # Constructing adjacency matrix for edges
    EMatrix = []
    for i in range(numoftreat):
	each_row = []
	for j in range(numoftreat):
            each_row.append(0)
	EMatrix.append(each_row)

    CombList = []
    for comb in combinations(range(1, (numoftreat + 1)), 2): CombList.append(tuple(comb))
    count = 0
    for m,n in CombList:

	if Pattern[count] == '1':
            EMatrix[n - 1][m - 1] = 0
	elif Pattern[count] == '2':
            EMatrix[n - 1][m - 1] = 1
	elif Pattern[count] == '0':
            EMatrix[m - 1][n - 1] = 1
	count = count + 1

    # Method for checking the source vertex
    def SourceVertex():
        for col in range(len(VList)):
            flag = True
            for row in range(len(VList)):
		if EMatrix[row][col] == 1:
                    flag = False
                    break
            if flag:
		return col
	return -1

    # Method for deleting a vertex from matrix
    def delVertex(vertex):
        if len(VList) > 1:
            #Move row up
            for row in range(vertex, len(VList) - 1):
		for col in range (len(VList)):
                    EMatrix[row][col] = EMatrix[row + 1][col]
            #Move col left			
            for col in range(vertex, len(VList) - 1):
		for row in range (len(VList) - 1):
                    EMatrix[row][col] = EMatrix[row][col + 1]

    # Checking whether graph is a complete graph or not
    A = AdjacencyMatrix(numoftreat,Pattern)
    numberofedges = 0
    for m in range(1, (numoftreat + 1)):
	for n in range(1, (numoftreat + 1)):
            if A.getitem(m,n) == 1:
                numberofedges = numberofedges + 1
	
    maxedges = 0
    maxedges = (numoftreat * (numoftreat - 1))/2

    #if complete graph
    if numberofedges == maxedges:

	#Topological Sorting
	TList = []
	while(len(VList) > 0):
            source=SourceVertex()
            if source == -1:
		return False
            else:
		TList.append(VList[source])
		delVertex(source)
		VList.pop(source)

	results = ""
	for i in range(len(TList)):
            results = "%s%s" % (results, TList[i])
            if i < len(TList) - 1:
                results = "%s%s" % (results, ' > ')
	return results

    #if not complete graph
    else:

	# Cheking whether graph is contractible or not
	if not Contractibility(numoftreat, Pattern):
            return "Not contractible"

	# Construct contractible vertices list
	A = AdjacencyMatrix(numoftreat,Pattern)

	CombList = []
	for comb in combinations(range(1, (numoftreat + 1)), 2): CombList.append(tuple(comb))

	# Storing equivalent vertices list
	EquivalentList = []
	for i,j in CombList:
            if A.getitem(i, j) == 0 and A.getitem(j, i) == 0:
		EquivalentList.append(tuple([i, j]))

	# Storing contractible vertices list
	ContractibleList = []
	for a,b in EquivalentList:
            condition = 1
            for l in range(1, numoftreat + 1):
		if A.getitem(l, a) == A.getitem(l, b) and A.getitem(a,l) == A.getitem(b, l):
                    pass
		else: 
                    condition = 0
                    break

            if condition == 1:
                ContractibleList.append([a, b])

	# Topological Sorting
	tempList = []
	for i in range(len(nameoftreat)):
            tempList.append(nameoftreat[i])

	TList = []
	while(len(VList) > 0):
            source=SourceVertex()
            if source == -1:
                return False
            else:
		CVList=[]
		CVList.append(VList[source])
		for a,b in ContractibleList:
                    if (tempList[a - 1] in CVList)|(tempList[b - 1] in CVList):
			CVList.append(tempList[a - 1])
			CVList.append(tempList[b - 1])
		CVList = list(set(CVList))
		CVList.sort()

		if (len(CVList) > 0):
                    TList.append(CVList)
                    for cv in CVList:
			i=VList.index(cv)
			delVertex(i)
			VList.remove(cv)
		else:
                    TList.append(VList[source])
                    delVertex(source)
                    VList.pop(source)

	results = ""
	for i in range(len(TList)):
            for j in range(len(TList[i])):
		results="%s%s" % (results, TList[i][j])
		if j < len(TList[i]) - 1:
                    results="%s%s" % (results, ' ~ ')
            if i < len(TList) - 1:
		results="%s%s" % (results, ' > ')
	return results


"""
Module creates .dot files and calls graphviz to create .svg files for overview cluster
"""

def Overviewclusters(SignificantGenesindex, Patterndict, numoftreat, numofrep, nameoftreat, sign_level, result_dir_path):

    file_path = os.path.join(result_dir_path, 'overview.dot')

    ovfile = open(file_path, 'w')
    ovfile.write("strict digraph TDAGS {\n")
    ovfile.write("\tlabelloc = \"t\";\n\tcompound = false;\n")
    ovfile.write("\tsubgraph cluster_0 {\n")

    """
    Sorting Dictionary containing pattern clusters containing probesets by number of genes in the cluster specified by pattern
    """
    Sortedlist = []
    Patterndictionarys = {}
    for keys in Patterndict.keys():
        Patterndictionarys[keys] = len(Patterndict[keys])
    Sortedlist = sorted(Patterndictionarys, key = Patterndictionarys.__getitem__, reverse = True)

    cluster = 1

    for pattern in Sortedlist:
	ovfile.write("\t\tsubgraph cluster_" + str(cluster) + " {\n")
	ovfile.write("\t\t\tlabel= \"" + str(pattern) + "\";\n");

	if Contractibility(numoftreat,str(pattern)):
            ovfile.write("\t\t\tstyle = filled;\n\t\t\t")
	else:
            pass

	for treat in nameoftreat:
            ovfile.write("node " + "[label= " + str(treat) + ", shape = plaintext, fontsize=12] " + str(treat) + str(cluster) + ";")
            ovfile.write("\n\t\t\t")

	for row in range(1, (numoftreat + 1)):
            for col in range(1, (numoftreat + 1)):
		if AdjacencyMatrix(numoftreat, str(pattern)).getitem(row,col):
                    ovfile.write(str(nameoftreat[row - 1]) + str(cluster) + "->" + str(nameoftreat[col - 1]) + str(cluster) + ";")
	
        ovfile.write("\t\t}\n")
        pattern = ""
        cluster = cluster + 1

    ovfile.write("}}")
    ovfile.close()

    img_path = os.path.join(result_dir_path, 'overview.svg')
    os.system("fdp -Tsvg  " + file_path + " -o " + img_path)


"""
Module creates .dot files and calls graphviz to create .png files for first clusters
"""
def FirstLevelClusters(pattern, number, numoftreat, nameoftreat, first_dir_path):

    file_path = os.path.join(first_dir_path, str(pattern) + 'graphviz.dot')
    clusterfile = open(file_path, 'w')
    clusterfile.write ("strict digraph DAGS {\n\tsize = \"4,4!\" ; ratio =\"fill\"; ")
    clusterfile.write("subgraph cluster_0{\n\t\t\t")
    clusterfile.write ("labeldoc = \"t\";\n\t\t\t")
    clusterfile.write ("label = \"" + str(pattern) + "\";")

    for treat in nameoftreat:
        clusterfile.write("node	" + "[label= " + str(treat) + ", shape = plaintext, fontsize=20] " + str(treat) + ";")
                
    for row in range(1, (numoftreat + 1)):
        for col in range(1, (numoftreat + 1)):
            if AdjacencyMatrix(numoftreat, str(pattern)).getitem(row, col):
                clusterfile.write("\n" + str(nameoftreat[row - 1]) + "->" + str(nameoftreat[col - 1]) + ";")

    clusterfile.write("\n\t}}")
    clusterfile.close()

    img_path = os.path.join(first_dir_path, str(pattern) + "graph.png")
    os.system("dot -Tpng " + file_path + " -o " + img_path)


"""
Module write gene information to file
"""
def WriteToFiles(Genetuple, Pvaluedict, FoldChangedict, Header, probeset, genesymbol, unigeneid, reppubid, otherinfo, numoftreat, nameoftreat, numofrep, Expressionmatrix, file_path):
    
    gene_file = open(file_path, 'w')

    for i in range(0, 5):
        gene_file.write(str(Header[i]) + "\t")

    gene_file.write("Pvalue\t")

    for i in range(5, len(Header)):
        gene_file.write(str(Header[i]) + "\t")

    for i in range(1, numoftreat - 1):
        gene_file.write(nameoftreat[i] + "(fold-change)\t")

    gene_file.write(nameoftreat[numoftreat-1] + "(fold-change)\n")

    for index, value in Genetuple:
        gene_file.write(str(probeset[int(index)]) + "\t" + str(genesymbol[int(index)]) + "\t" + str(unigeneid[int(index)]) + "\t"+ str(reppubid[int(index)]) + "\t" + str(otherinfo[int(index)]) + "\t")

        gene_file.write(str(Pvaluedict[int(index)]) + "\t")

	for i in range(0, len(Expressionmatrix[int(index)])):
            gene_file.write(str(Expressionmatrix[int(index)][i]) + "\t")

        fcNum = len(FoldChangedict[int(index)])
        for i in range(0, fcNum - 1):
            gene_file.write(str(FoldChangedict[int(index)][i]) + "\t")
        gene_file.write(str(FoldChangedict[int(index)][fcNum - 1]) + "\n")

    gene_file.close()


"""
Module creates gene list files containing all the information for all significant genes 
"""
def AllGeneListFiles(geneIndex, Pvaluedict, FoldChangedict, Header, probeset, genesymbol, unigeneid, reppubid, otherinfo, numoftreat, nameoftreat, numofrep, Expressionmatrix, result_dir_path):

    Meanvalues = {}
    for s in geneIndex:
        Meanvalues[s] = Mean(Expressionmatrix[s])
    Genetuple = {}
    Genetuple = sortbyvalue(Meanvalues, reverse = True)

    file_path = os.path.join(result_dir_path, "AllGeneList.txt")
    WriteToFiles(Genetuple, Pvaluedict, FoldChangedict, Header, probeset, genesymbol, unigeneid, reppubid, otherinfo, numoftreat, nameoftreat, numofrep, Expressionmatrix, file_path)


def SignificantGeneListFiles(SignificantGenesindex, Pvaluedict, FoldChangedict, Header, probeset, genesymbol, unigeneid, reppubid, otherinfo, numoftreat, nameoftreat, numofrep, Expressionmatrix, result_dir_path):

    Meanvalues = {}
    for s in SignificantGenesindex:
        Meanvalues[s] = Mean(Expressionmatrix[s])
    Genetuple = {}
    Genetuple = sortbyvalue(Meanvalues, reverse = True)

    file_path = os.path.join(result_dir_path, "SignificantGeneList.txt")
    WriteToFiles(Genetuple, Pvaluedict, FoldChangedict, Header, probeset, genesymbol, unigeneid, reppubid, otherinfo, numoftreat, nameoftreat, numofrep, Expressionmatrix, file_path)


"""
Module creates gene list files containing all the information for first clusters 
"""
def PatternGeneListFiles(pattern, IndexPatterndict, Pvaluedict, FoldChangedict, Header, probeset, genesymbol, unigeneid, reppubid, otherinfo, numoftreat, nameoftreat, numofrep, Expressionmatrix, first_dir_path):

    Meanvalues = {}
    for s in IndexPatterndict[pattern]:
        Meanvalues[s] = Mean(Expressionmatrix[s])

    Genetuple = {}
    Genetuple = sortbyvalue(Meanvalues, reverse = True)

    file_path = os.path.join(first_dir_path, str(pattern) + "GeneList.txt")
    WriteToFiles(Genetuple, Pvaluedict, FoldChangedict, Header, probeset, genesymbol, unigeneid, reppubid, otherinfo, numoftreat, nameoftreat, numofrep, Expressionmatrix, file_path)


"""
Module creates gene list files containing all the information for second clusters 
"""
def MetaPatternGeneListFiles(symbol, Symboldict, Resultdict, IndexPatterndict, Pvaluedict, FoldChangedict, Header, probeset, genesymbol, unigeneid, reppubid, otherinfo, numoftreat, nameoftreat, numofrep, Expressionmatrix, second_dir_path):

    Patterns = Resultdict[symbol]
    Meanvalues = {}
    for pattern in Patterns:
        for s in IndexPatterndict[pattern]:
            Meanvalues[s] = Mean(Expressionmatrix[s])

    Genetuple = {}
    Genetuple = sortbyvalue(Meanvalues, reverse = True)

    file_path = os.path.join(second_dir_path, str(Symboldict[symbol]) + "GeneList.txt")
    WriteToFiles(Genetuple, Pvaluedict, FoldChangedict, Header, probeset, genesymbol, unigeneid, reppubid, otherinfo, numoftreat, nameoftreat, numofrep, Expressionmatrix, file_path)



"""
Module produces result dictionary for second clusters
It categorizes patterns into clusters depending on different cases as discussed in paper by Phan et.al ,2009.
Input: name of treatment, number of treatment, Patterndictionary (containing patterns and probeset)
Output: dictionary containing second clusters.
"""
import copy

def SecondLevelClusterDictionary(nameoftreat,numoftreat,Patterndictionary):
    Patterndict = Patterndictionary

    Combinationtreat = []
    for x in combinations(nameoftreat, 2): Combinationtreat.append(list(x))

    SecondLevelClusterDictionary = []

    for treat in nameoftreat:
		SecondLevelClusterDictionary.append(treat)

    for m in range(2, numoftreat):
		for n in combinations(nameoftreat,m): SecondLevelClusterDictionary.append(list(n))

    Cluster_treat_matrixG = []	

    for i in SecondLevelClusterDictionary:
		Cluster_treat_matrixG.append(['']*len(Combinationtreat))

    count = 0	
    for x in SecondLevelClusterDictionary:
		for y in range(0,len(Combinationtreat)):
			if set(Combinationtreat[y]).issubset(set(x)):
				Cluster_treat_matrixG[count][y] = '1'
		count = count + 1
	

    Cluster_treat_matrixL = copy.deepcopy(Cluster_treat_matrixG[:])

    case = 0
    statusG = ['']*(numoftreat*(numoftreat-1)/2)
    statusL = statusG[:]

    for a in SecondLevelClusterDictionary:
	for b in range(0,len(Combinationtreat)):

            if type(a) != type([]):
                if a in Combinationtreat[b]:
                    if Cluster_treat_matrixG[case][b] == '':
			if statusG[b] == '':
                            Cluster_treat_matrixG[case][b] = '0'
                            statusG[b] = '0'
                            Cluster_treat_matrixL[case][b] = '2'
                            statusL[b] = '2'
			else:
                            Cluster_treat_matrixG[case][b] = '2'
                            statusG[b] = '2'
                            Cluster_treat_matrixL[case][b] = '0'
                            statusL[b] = '0'
					
            else:
		for c in a:
                    if c in Combinationtreat[b]:
			if Cluster_treat_matrixG[case][b] == '':
                            if Combinationtreat[b].index(c) == 0:
				Cluster_treat_matrixG[case][b] = '0'
				Cluster_treat_matrixL[case][b] = '2'
                            else:
				Cluster_treat_matrixG[case][b] = '2'
				Cluster_treat_matrixL[case][b] = '0'				
				
	case = case + 1

    Resultdict = {}
    Symboldict = {}

    for pattern in Patterndict:
	pcount = 0
	for cluster in Cluster_treat_matrixG:
            counts = 0
            for p in pattern:
		condition = 1
		if p == cluster[counts] or cluster[counts] == '':
                    pass
		else:
                    condition = 0
                    break
		counts = counts + 1

            if condition == 1:
		p = SecondLevelClusterDictionary[pcount]
		q = SecondLevelClusterDictionary[pcount]

		if type(p) == type([]):
                    p = str('~'.join(p))
                    q = str('_'.join(q))

		keyn = str(p) +'>'

		if keyn in Resultdict.keys():
                    Resultdict[keyn].append(pattern)
		else:
                    Resultdict[keyn] = []
                    Symboldict[keyn] = str(q) + "G"
                    Resultdict[keyn].append(pattern)

            pcount = pcount + 1

    for pattern in Patterndict:
	pcount = 0
	for cluster in Cluster_treat_matrixL:
            counts = 0
            for p in pattern:
		condition = 1
		if p == cluster[counts] or cluster[counts] == '':
                    pass
		else:
                    condition = 0
                    break
		counts = counts + 1

            if condition == 1:
		p = SecondLevelClusterDictionary[pcount]
		q = SecondLevelClusterDictionary[pcount]

		if type(p) == type([]):
                    p = str('~'.join(p))
                    q = str('_'.join(q))

		keyn = '>' + str(p)

		if keyn in Resultdict.keys():
                    Resultdict[keyn].append(pattern)
		else:
                    Resultdict[keyn] = []
                    Symboldict[keyn] = str(q) + "L"
                    Resultdict[keyn].append(pattern)

            pcount = pcount + 1	

    return Resultdict, Symboldict


def MetaCluster(nameoftreat, numoftreat, pattern):

    Pattern = pattern

    Combinationtreat = []
    for x in combinations(nameoftreat, 2): Combinationtreat.append(list(x))

    MetaCluster = []

    for treat in nameoftreat:
	MetaCluster.append(treat)

    for m in range(2, numoftreat):
	for n in combinations(nameoftreat, m): MetaCluster.append(list(n))

    Cluster_treat_matrixG = []

    for i in MetaCluster:
	Cluster_treat_matrixG.append([''] * len(Combinationtreat))

    count = 0
    for x in MetaCluster:
	for y in range(0, len(Combinationtreat)):
            if set(Combinationtreat[y]).issubset(set(x)):
		Cluster_treat_matrixG[count][y] = '1'
	count = count + 1


    Cluster_treat_matrixL = copy.deepcopy(Cluster_treat_matrixG[:])

    case = 0
    statusG = [''] * (numoftreat * (numoftreat-1)/2)
    statusL = statusG[:]

    for a in MetaCluster:
	for b in range(0, len(Combinationtreat)):

            if type(a) != type([]):
                if a in Combinationtreat[b]:
                    if Cluster_treat_matrixG[case][b] == '':
			if statusG[b] == '':
                            Cluster_treat_matrixG[case][b] = '0'
                            statusG[b] = '0'
                            Cluster_treat_matrixL[case][b] = '2'
                            statusL[b] = '2'
			else:
                            Cluster_treat_matrixG[case][b] = '2'
                            statusG[b] = '2'
                            Cluster_treat_matrixL[case][b] = '0'
                            statusL[b] = '0'
					
            else:
		for c in a:
                    if c in Combinationtreat[b]:
			if Cluster_treat_matrixG[case][b] == '':
                            if Combinationtreat[b].index(c) == 0:
				Cluster_treat_matrixG[case][b] = '0'
				Cluster_treat_matrixL[case][b] = '2'
                            else:
				Cluster_treat_matrixG[case][b] = '2'
				Cluster_treat_matrixL[case][b] = '0'				
				
	case = case + 1

    metacluster = ''

    pcount = 0
    for cluster in Cluster_treat_matrixG:
        counts = 0
        for p in Pattern:
            condition = 1
            if p == cluster[counts] or cluster[counts] == '':
                pass
            else:
                condition = 0
                break
            counts = counts + 1

        if condition == 1:
            p = MetaCluster[pcount]
            if type(p) == type([]):
                p = str('~'.join(p))
            metacluster = str(p) +'>'

        pcount = pcount + 1

    pcount = 0
    for cluster in Cluster_treat_matrixL:
        counts = 0
        for p in Pattern:
            condition = 1
            if p == cluster[counts] or cluster[counts] == '':
                pass
            else:
                condition = 0
                break
            counts = counts + 1

        if condition == 1:
            p = MetaCluster[pcount]
            if type(p) == type([]):
                p = str('~'.join(p))
            metacluster = '>' + str(p)

        pcount = pcount + 1	

    return metacluster


"""
Module creates .dot files and calls graphviz to create .png files for second clusters
"""

def SecondLevelClusters(second_cluster, Allpattern, Contractratio, Resultdict, Symboldict, Sumdict, Patterndict, numoftreat, nameoftreat, second_dir_path):

    file_path = os.path.join(second_dir_path, str(Symboldict[second_cluster]) + "scgraphviz.dot")
    slc_file = open(file_path, 'w')
    slc_file.write("strict digraph FirstlevelCluster {\n")
    slc_file.write("\tlabelloc = \"t\";\n\tcompound = false;\n")

    slc_file.write("\tsubgraph cluster_0 {\n")
    slc_file.write("\tlabel= \"" + str(second_cluster) + "\";" + "\n")
    
    cluster = 1
    clustid = 0

    sortvaluelist = []
    sortvaluelist = sorted(Allpattern, key = Allpattern.__getitem__, reverse= True)  

    for	pattern	in sortvaluelist:
	cluster	= cluster + 1 
        clustid = clustid + 1
	slc_file.write("\t\tsubgraph	cluster_" + str(cluster) + " {\n")
	slc_file.write("\t\tlabel = \"" + str(pattern) + "\";\n")

	if Contractibility(numoftreat, str(pattern)):
	    slc_file.write("\t\t\tstyle = filled;\n\t\t\t")
	else:
	    slc_file.write("\t\t\t")

	for treat in nameoftreat:
	    slc_file.write("node	" + "[label = "+str(treat) + ", shape = plaintext, fontsize = 16] " + str(treat) + str(cluster) + ";")
	    slc_file.write("\n\t\t\t")

	for row	in range(1, (numoftreat + 1)):
	    for col in range(1, (numoftreat + 1)):
		if AdjacencyMatrix(numoftreat, str(pattern)).getitem(row,col):
		    slc_file.write(str(nameoftreat[row - 1]) + str(cluster) + "->" + str(nameoftreat[col - 1]) + str(cluster) + ";")		
	
    	slc_file.write("\t\t\t}\n")

    slc_file.write("}}")
    slc_file.close()

    img_path = os.path.join(second_dir_path, str(Symboldict[second_cluster]) + "scgraph.svg")
    
    os.system("fdp -Tsvg " + file_path + " -o " + img_path)
