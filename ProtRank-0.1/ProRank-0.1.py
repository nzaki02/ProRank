#!/usr/bin/python0
# c 2011 Jos Berenger, Nazar Zaki and Dmitry Efimov, 2011.JUN.20
# UAE University, FIT, Bioinformatics Lab more on @bioAE Twitter
# LIST OF CHANGES: added POSIX

from optparse import OptionParser

from collections import defaultdict
from math import *
import string
from numpy import *
import numpy as num
import numpy.matlib as M
from numpy.matlib import rand,zeros,ones,empty,eye

def sortedDictValues1(adict):
    items = adict.items()
    items.sort()
    return [value for key, value in items]

def hermitian(A, **kwargs):
    return num.transpose(A,**kwargs).conj()
T = num.transpose
H = hermitian


# build dictionary of unique proteins and creates a relations dictlist
# parseProteins(myproteins, id2protein, usedproteins, myrelations, "data/interactions.txt")

def parseProteins(p2id, id2p, usp, rel, pairsfile,options,args):
    inp = open (pairsfile,"r");
    k= 0
    for line in inp.readlines():
        p = string.split(line,"\t")
	x = p[0].strip();
        y = p[1].strip();
        #print "cheking pair  " + x + " " + y 
        if not p2id.has_key(x):
            #print "k =" + str(k) + " found new protein myproteins[" + str(len(p2id)) + "]  = " + x ;
            p2id[x] = len(p2id);
            usp[len(id2p)] = 1
            id2p[len(id2p)] = x
        if not p2id.has_key(y):
            #print "k =" + str(k) + " found new protein myproteins[" + str(len(p2id)) + "]  = " + y ;
            p2id[y] = len(p2id);
            usp[len(id2p)] = 1
            id2p[len(id2p)] = y
        if not p2id[x] in rel[p2id[y]]:
            rel[p2id[y]].append(p2id[x]);
            rel[p2id[x]].append(p2id[y]);
            k = k +1
        
    print "Number of unique interactions = ", k

    clean_low_realtions_less_than = 0

    
    ignored = 0
    maxlen = 0
    for k, v in rel.items():
        #print "k = ", k , "len ", len(v)," values = ", v
        if (len(v) <= int(clean_low_realtions_less_than)):
            #print "to erase", len(v), "  ",clean_low_realtions_less_than
            del rel[k]
            ignored = ignored + 1
        if maxlen < len(v):
            maxlen = len(v)
            curprotein = k
        #print "protein", k, "has connections", v
    #print "max number of neighbors has protein", id2p[curprotein]     
    #print ignored, " pp interactions with less than or totaling ", clean_low_realtions_less_than, " interactions to a single protein have been ignored"
    print ""

def parseWeights(p2id, pairsfile, W):
    inp = open (pairsfile,"r");
    for line in inp.readlines():
        p = string.split(line,"\t")
        x = p2id[p[0].strip()];
        y = p2id[p[1].strip()];
        weight = p[2].strip();
        W[x,y] = weight
        W[y,x] = weight

def createMatrixOfConnectivity(rel,N,A):
    for key in rel.iterkeys():
        x = key
        y = rel[key]
        A [x,y] = 1
        A [y,x] = 1

def makeSMatrix(p2id, Sfile, S):
    # ------ Read Similarity Matrix ------------
   # print "Checking if 1st row of similaritymatrix.csv is matches interactions.txt ... "
    sim_proteins = {}
    id2sim_proteins = {}
    inp = open (Sfile,"U");
    k = -1
    for line in inp.readlines():
        k = k + 1
        p = string.split(line,",")
        x = p[0].strip();
        if p2id.has_key(x):
            if not sim_proteins.has_key(x):
                #print " found new protein sim_proteins[" + str(len(sim_proteins)) + "]  = " + x ;
                sim_proteins[x] = k;
                id2sim_proteins[k] = x;
            else:
                print x + " protein row is repeated in file (!)"
                exit(0)
        else:
            print "Warning found a protein " + x +" in data/similaritymatrix.csv not found in interactions.txt (!)"
            print "missmatch! corrrect the error. This row will be ignored."
    
    #print "Importing..."
    inp = open (Sfile,"U");
    sim_cols = -1
    sim_rows = 0
    for line in inp.readlines():
        sim_rows = sim_rows + 1 
        p = string.split(line,",")
        sim_cols= len(p)-2
        if (p2id.has_key(p[0])): 
            i = p2id[p[0]]
            for k in range(2,len(p)):
                if ( id2sim_proteins.has_key(k-2) and p2id.has_key(id2sim_proteins[k-2])):
                    j = p2id[id2sim_proteins[k-2]] # k-2 ---> proteinname --> id of interactions
                    S[i,j] = 1.0*float(p[k].replace(",",".")) #russian decimal notation just in case
   #print "S matrix cols = " + str(S.shape[1]) + " rows = "  + str(S.shape[0]) 
    print "Similarity matrix cols = " + str(sim_cols) + " rows = " + str(sim_rows)
    if (sim_cols !=  S.shape[1] or sim_rows != S.shape[0]):
        print "Mismatch warning. Mismatch positions might be filled with 1.0 values (!)"

    for i in range(0,sim_rows):
        for j in range(0,sim_rows):
            if (i == j):
                S[i,j] = 0
            else:
                S[i,j] = S[i,j]/100.0

    #return S
        

#make matrix A using similarities algorithm 1
def makeMatrixSimilarities1(A,rel,N,S,m):
    # Normalize columns to sum no more than 1 ------------
    sumbycolumn = A.sum(1);
    for i in range(0,N):
        for j in range(0,N): #according to PageRank algorithm
            if A[i,j] == 0:
                A[i,j] = S[i,j] #probability of connection equals probability from similarities table
            else:
                if sumbycolumn[j,0] <> 0:
                    A[i,j]=S[i,j] + (A[i,j]/sumbycolumn[j,0]) #probability of connection equals probability form similarities table plus 1/number of neighbors

    sumbycolumn = A.sum(1);
    for i in range(0,N):
        for j in range(0,N): #normalization of matrix A
            if sumbycolumn[j,0] <> 0:
                A[i,j]=A[i,j]/sumbycolumn[j,0]

    #savetxt("matrixA.txt",  A)
    return A


#make matrix A using similarities algorithm 2
def makeMatrixSimilarities2(A,rel,N,S,m):
    #change matrix S according to number of possible interactions
    for i in range(0,N):
        sumbyrow = 0;
        for j in range(0,N):
            if A[i,j] <> 0:
                sumbyrow = sumbyrow + S[i,j];
        for j in range(0,N):
            if A[i,j] == 0:
                S[i,j] = 0
            else:
                if (sumbyrow <> 0):
                    S[i,j] = S[i,j]/sumbyrow

    # Normalize columns using matrix S
    for i in range(0,N):
        for j in range(0,N): #according to PageRank algorithm
            if A[i,j] <> 0:
                A[i,j]=S[j,i]

    #savetxt("matrixA.txt",  A)
    return A

#algorithm using just similarities
def algorithmS(A,N,rel,p2c):
    useThres = raw_input('Choose the threshold for probability:')
    if (useThres):
        thre = useThres
    else:
        thre = 0.5
    pr = defaultdict(list) #array for complexes

    for i in range(0,N):
        newcomplex = 0;
        for j in range(0,N):
            #if (j<=i):
            #    continue
            if A[i,j] < thre:
                A[i,j] = 0;
                continue
            if A[i,j]<>0:
                if newcomplex == 0:
                    newcomplex = 1;
                    pr[i].append(i);
                    p2c[i] = i
                pr[i].append(j);
                p2c[j] = i;
                A[i,j] = 0;
                A[j,i] = 0
    return pr

def graphIsConnected(G,Nn):
    E = zeros((Nn,Nn))
    for i in range(0,Nn):
        E[i,i] = 1
        
    TG = E + G
    Gcur = G.copy()
    for i in range(0,Nn):
        Gcur = Gcur*G
        TG = TG + Gcur
    
    for i in range(0,Nn):
        for j in range(0,Nn):
            if TG[i,j] == 0:
                return 0
    
    return 1

def findBridgeProteins(id2p, rel, p_f):
    for protein1 in id2p:
        set_of_neigh = set(rel[protein1])
        Nn = len(set_of_neigh)
        G = zeros((Nn,Nn))
        list_of_neigh = {}
        k = 0
        for id in set_of_neigh:
            list_of_neigh[k] = id
            k = k + 1

        for i in range(0,Nn):
            for j in range(i,Nn):
                if list_of_neigh[i] in rel[list_of_neigh[j]]:
                    G[i,j] = 1
                    G[j,i] = 1

        flag_not_connected = 0
        sumbyrow = G.sum(0);
        for i in range(0,Nn):
            if sumbyrow[0,i]*1.0/Nn < 0.2:
                flag_not_connected = 1
        if flag_not_connected == 1:
            p_f[protein1] = 2
        else:
            p_f[protein1] = 2-graphIsConnected(G,Nn)


def createWeightMatrix(A,S,N,p_f):
    for id in p_f:
        if p_f[id]<2:
            continue
        A[N-1,id] = 1
        A[id,N-1] = 1
        S[N-1,id] = 1000000.0
        S[id,N-1] = 1000000.0
    
def eigenCalc(xinitial,MM,N,maxiter):
    x = xinitial;
    x_norma = T(x)*x;
    xn = math.sqrt(x_norma[0,0]);
    k=0
    prev_xn = 0
    err = 100000
    for i in range(0, maxiter):
        k = k +1
        x = MM *x;
        x_norma = T(x)*x;
        xn = math.sqrt(x_norma[0,0]);
        #if (xn < 0.001):
        #    break; 
        #x = x/xn;
        err = xn - prev_xn
        prev_xn = xn

    print "Converged after " + str(k) +  " iterations" + " with err = " + str(err);
    return x

def PageRank(A,N,m,id2p,top,p_f):
    E = ones((N,N));
    MM = (1.0-m)*A + (m/N)*E;
    #MM = A;

    x = MM * (1.0/N) * ones((N,1));
    x = eigenCalc(x, MM,N,200)

    for i in range(0,N):
        top[i] = x[i][0,0]

    import operator
    s_t = sorted(top.iteritems(), key=operator.itemgetter(1))
    s_t.reverse();
    #print "sorted_top = ", sorted_top;
    f = open('listwithranks.txt', 'w')
    k = 0
    for protein in s_t:
        s = str(id2p[protein[0]]) + " " + str(protein[1]) + " flag " + str(p_f[protein[0]])
        print >> f, s
    f.close()

    return s_t
    
def findcomplexes(A,N,up,rel,id2p,p2c,s_t,top,p_f):

   #print "\n -- Predicted Complexes by rank score  ---- "
    #find all complexes according to the ranks
    pr = defaultdict(list) #array for complexes
    complexranks = {}
    k=1
    #print p_f
    for candidate in s_t:
        rank = candidate[1]
        c_key = candidate[0]
        if p_f[c_key]>1:
            continue
        if (up[c_key] == 0):
            continue
        k = k+1
        pr[k].append(c_key)
        complexranks[k] = rank
        p2c[c_key] = k
        up[c_key] = 0

        for id in rel[c_key]:
            #if (top[id]>rank_of_threshold):
            #    continue
            if (up[id] == 0):
                continue
            s = id2p[id]
            p2c[id] = k
            pr[k].append(id)
            up[id] = 0

    return pr


def findcomplexes_variance(A,N,up,rel,id2p,p2c,s_t):
    print "\n -- Predicted Complexes by rank score  ---- "
#find all complexes according to the ranks
    pr = defaultdict(list) #array for complexes
    complexranks = {}
    proteinranks = {}
    proteinvariance = {}
    k=0
    for candidate in s_t:
        c_key = candidate[0]
        c_rank = candidate[1]
        proteinranks[c_key] = c_rank
    
    for candidate in s_t:
        c_key = candidate[0]
        c_rank = candidate[1]
        measure_variance = 0
        for id in rel[c_key]:
            n_rank = proteinranks[id]
            measure_variance = measure_variance + (c_rank-n_rank)*(c_rank-n_rank)
        measure_variance = sqrt(measure_variance)
        proteinvariance[c_key] = measure_variance

    import operator
    s_t_var = sorted(proteinvariance.iteritems(), key=operator.itemgetter(1))
    s_t_var.reverse();
    
    for candidate in s_t_var:
        k = k+1
        c_key = candidate[0]
        if (up[c_key] == 0):
            continue
        pr[k].append(c_key)
        complexranks[k] = candidate[1]
        p2c[c_key] = k
        up[c_key] = 0
        
        for id in rel[c_key]:
            
            if (up[id] == 0):
                continue
            s = id2p[id]
            p2c[id] = k
            pr[k].append(id)
            up[id] = 0
    return pr

def mergeComplexesWithTwoProteins(sorted_complex_len,myrelations,predicted,protein2complex):
    #choose all complexes with two proteins
    for complnumber in sorted_complex_len:
        if (complnumber[1] != 2):
            continue
        key_complex = complnumber[0]
        neigh1 = set(myrelations[predicted[key_complex][0]]);
        neigh2 = set(myrelations[predicted[key_complex][1]]);
        inter = neigh1 & neigh2
        if len(inter) == 1:
            for common_neigh in inter:
                if protein2complex.has_key(common_neigh):
                    predicted[protein2complex[common_neigh]].append(predicted[key_complex][0]);
                    predicted[protein2complex[common_neigh]].append(predicted[key_complex][1]);
                    del predicted[key_complex]

def mergeComplexesWithThreeProteins(sorted_complex_len,myrelations,predicted,protein2complex):
    #choose all complexes with three proteins
    for complnumber in sorted_complex_len:
        if (complnumber[1] != 3):
            continue
        key_complex = complnumber[0]
        neigh1 = set(myrelations[predicted[key_complex][0]]);
        neigh2 = set(myrelations[predicted[key_complex][1]]);
        neigh3 = set(myrelations[predicted[key_complex][2]]);
        inter = neigh1 & neigh2 & neigh3
        if len(inter) == 1:
            for common_neigh in inter:
                if protein2complex.has_key(common_neigh):
                    predicted[protein2complex[common_neigh]].append(predicted[key_complex][0]);
                    predicted[protein2complex[common_neigh]].append(predicted[key_complex][1]);
                    predicted[protein2complex[common_neigh]].append(predicted[key_complex][2]);
                    del predicted[key_complex]


def mergeComplexesWithOneProtein(sorted_complex_len,predicted,myrelations,sorted_top):
    #choose all complexes with one protein
    for complnumber in sorted_complex_len:
        if (complnumber[1] != 1):
            continue
        key_complex = complnumber[0]
        id1 = predicted[key_complex][0];
        #choose the maximum rank neighbor of the protein in the low rank complex
        maxrank = -1;
        idmaxrank = -1;
        for id2 in myrelations[id1]:
            if (sorted_top[id2] > maxrank):
                idmaxrank = id2;
        #add protein to other complex
        try:
            predicted[protein2complex[idmaxrank]].append(id1);
        except:
            idmaxrank = -1
        #delete the complex with this core protein from predicted
        del predicted[key_complex];
        
        
def deleteSmallComplexes(sorted_complex_len,usedproteins,predicted):
    for complnumber in sorted_complex_len:
        if (complnumber[1] > 2):
            continue
        key_complex = complnumber[0]
        #usedproteins[predicted[key_complex][0]] = 1;
        try:
            del predicted[key_complex]
        except:
            continue

def mergeComplexes(predicted,rel):
    a = set(predicted.iterkeys())
    b = set(predicted.iterkeys())
    for c1 in a:
        for c2 in b:
            if c1 == c2:
                continue
            if not predicted.has_key(c1):
                continue
            if not predicted.has_key(c2):
                continue
            complex1 = set(predicted[c1])
            complex2 = set(predicted[c2])
            merge = 1
            for id in complex1:
                neigh_of_id = set(rel[id])
                if (2*len(neigh_of_id & complex2)<len(complex2)):
                    merge = 0
            if merge == 1:
                for id in complex1:
                    predicted[c2].append(id)
                del predicted[c1]

def main():
 
    use = "Usage: %prog -i [interaction_file] -s [similarity_matrix_file] "
    parser = OptionParser(usage = use)
    parser.add_option("-i", "--interactions", dest="interactions", action="store", default="interactions.txt", help="name of interactions file.txt ")
    parser.add_option("-s", "--similaritymatrix", dest="similaritymatrix", action="store", default="similaritymatrix.csv", help="name of file that contains similarity matrix in csv format .cvs")

    
    options, args = parser.parse_args()        
            
    myproteins = {};
    id2protein = {};
    protein2complex = {}
    usedproteins = {}; #list for checking of used proteins
    sorted_top = {}; #list of ranks for proteins
    myrelations = defaultdict(list)
    predicted = defaultdict(list)
    protein_flags = {}; #for bridge proteins
    top = {};
    m = 0.1
    parseProteins(myproteins, id2protein, usedproteins, myrelations, options.interactions,options,args)
    N = len(myproteins);

    A = zeros((N,N))
    P = zeros((N,N))
    S = zeros((N,N))
    W = zeros((N,N))
    MM = zeros((N,N))
    
    #parseWeights(myproteins, "data/interactions.txt", W)

    # Apply formula  PageRank ---------------------------
    makeSMatrix(myproteins, options.similaritymatrix,S)
    findBridgeProteins(id2protein,myrelations,protein_flags)
    createMatrixOfConnectivity(myrelations,N,A)

    A = makeMatrixSimilarities2(A,myrelations,N,S,m)

    sorted_top = PageRank(A,N,m,id2protein,top,protein_flags)
    predicted = findcomplexes(A,N,usedproteins,myrelations,id2protein,protein2complex,sorted_top,top,protein_flags)
    
    #here we need to merge complexes with one or two proteins
    complex_len = {}
    sorted_complex_len = {}
    for key in predicted.iterkeys():
        complex_len[key] = len(predicted[key])
    import operator
    sorted_complex_len = sorted(complex_len.iteritems(), key=operator.itemgetter(1))
    
    mergeComplexes(predicted,myrelations)
    deleteSmallComplexes(sorted_complex_len,usedproteins,predicted)
            
    #print "predicted = ", predicted
    
    # ---------- add to the file result.txt ----------

    f = open('Predicted_Complexes.txt', 'w')
    k = 0
    for compl in predicted.iterkeys():
        k = k + 1
        s = "C" + str(k) + ": "
        for prot in predicted[compl]:
            s = s + " "+ id2protein[prot]
        print >> f, s
    #print s
    f.close()

main()
import os
