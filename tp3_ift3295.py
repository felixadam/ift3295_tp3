#!/usr/bin/env python
# coding: utf-8

# In[2]:


import numpy as np
import math

readsFile = open("tRNAs.fasta", "r")
readsLines = readsFile.readlines()
readsFile.close()

reads = []
for i in range(1, len(readsLines), 2):
    reads.append(readsLines[i][:-1])


# In[3]:


def enumerateWords(seq,length):
    numberOfWords = len(seq)-length+1
    words = []
    for i in range(numberOfWords):
        words.append(seq[i:i+length])
    return words


# In[4]:


def findMatches(kmers,word,kernel):
    matches = []
    for i in range(len(kmers)):
        for j in range(len(word)):
            if((kernel[j] == '1') and (word[j] != kmers[i][j])):
                break
            if(j == len(word)-1):
                matches.append(i)
    return matches


# In[5]:


kmers = enumerateWords("abcdefghi",2)
print(findMatches(kmers,"de","11"))


# In[6]:


def extendHsp(seq1,i,seq2,j,kernel,e):
    
    nbOnes = 0
    for k in range(len(kernel)):
        if(kernel[k] == "1"):
            nbOnes += 1
    
    #score = nbOnes*5 + (len(kernel)-nbOnes)*-4
    score = nbOnes*5
    scoreToReturn = score
    scoreMax = score
    
    start1 = i
    end1 = i+len(kernel)-1
    start2 = j
    end2 = j+len(kernel)-1
    
    canExtendLeft = True if (start1 > 0 and start2 > 0) else False
    canExtendRight = True if (end1 < len(seq1)-1 and end2 < len(seq2)-1) else False
    
    while((canExtendLeft or canExtendRight) and (score > scoreMax-e)):
        scoreLeft = -float("Inf")
        scoreRight = -float("Inf")
        scoreBoth = -float("Inf")
        
        #print(score,start1,end1,start2,end2)
        if(canExtendLeft):
            scoreLeft = score + (5 if (seq1[start1-1] == seq2[start2-1]) else -4)
            
        if(canExtendRight):
            scoreRight = score + (5 if (seq1[end1+1] == seq2[end2+1]) else -4)
            
        if(canExtendRight and canExtendLeft):
            scoreBoth = scoreLeft + scoreRight -score
            
        scores = [scoreLeft,scoreRight,scoreBoth]
        #print(scores)
        indexScoreMax = np.argmax(np.array(scores))
        score = max(scoreLeft,scoreRight,scoreBoth)
        #print(score)
        
        if(score > scoreMax):
            scoreMax = score
        
        if(score > scoreMax-e):
            scoreToReturn = score
            if(indexScoreMax == 0 or indexScoreMax == 2):
                start1 -= 1
                start2 -= 1
            if(indexScoreMax == 1 or indexScoreMax == 2):
                end1 += 1
                end2 += 1
            
        canExtendLeft = True if (start1 > 0 and start2 > 0) else False
        canExtendRight = True if (end1 < len(seq1)-1 and end2 < len(seq2)-1) else False
        
    return (start1,end1,start2,end2,scoreToReturn)


# In[7]:


#seq1 = "abcdefghi"
#seq2 = "aaaaafghijklmn"

seq1 = "AGTTGGTTAGAGTATACGCCTGATAAGCGTA"
seq2 = "AGTCGGCTAACGCATACGCTTGATAAGCGTA"

print(extendHsp(seq1,20,seq2,20,"11111111111",5))


# In[8]:


def computeBitscore(score):
    return round((0.192*score-math.log(0.176))/math.log(2))


# In[9]:


def computeEvalue(m,n,b):
    return m*n*(2**-b)


# In[10]:


def computeScore(seq1,start1,end1,seq2,start2,end2):
    length = end1 - start1 + 1
    score = 0
    for i in range(length):
        if(seq1[start1+i] == seq2[start2+i]):
            score += 5
        else:
            score -= 4
            
    return score


# In[11]:


def mergeHsps(hsps,seq1,seq2):
    
    merged_hsps_new = hsps
    merged_hsps = []
    while(len(merged_hsps) != len(merged_hsps_new)):
        merged_hsps = merged_hsps_new
        merged_hsps_new = []
        hasMerged = np.zeros((len(merged_hsps)))
        
        for i in range(len(merged_hsps)-1):
            for j in range(i+1,len(merged_hsps)):
                start1_1 = merged_hsps[i][0]
                end1_1 = merged_hsps[i][1]
                start2_1 = merged_hsps[i][2]
                end2_1 = merged_hsps[i][3]
                
                start1_2 = merged_hsps[j][0]
                end1_2 = merged_hsps[j][1]
                start2_2 = merged_hsps[j][2]
                end2_2 = merged_hsps[j][3]
                                
                if(start1_1 < start1_2 and end1_1 > start1_2):
                    overlap = end1_1 - start1_2
                    if(end2_1 - overlap  == start2_2):
                        hasMerged[i] = 1
                        hasMerged[j] = 1
                        score = computeScore(seq1,start1_1,end1_2,seq2,start2_1,end2_2)
                        newHsp = (start1_1, end1_2, start2_1, end2_2, score)
                        merged_hsps_new.append(newHsp)
                elif(start1_1 > start1_2 and end1_2 > start1_1):
                    
                    overlap = end1_2 - start1_1
                    if(end2_2 - overlap == start2_1):
                        hasMerged[i] = 1
                        hasMerged[j] = 1
                        score = computeScore(seq1,start1_2,end1_1,seq2,start2_2,end2_1)
                        newHsp = (start1_2, end1_1, start2_2, end2_1, score)
                        merged_hsps_new.append(newHsp)

        for k in range(len(merged_hsps)):
            if(hasMerged[k] == 0):
                merged_hsps_new.append(merged_hsps[k])
                
    return merged_hsps_new
                


# In[12]:


hsps = [(0, 3, 9, 12, 11),(2,5,11,14,4), (0, 11, 15, 26, 15), (0, 3, 24, 27, 11), (1, 10, 0, 9, 23), (9, 13, 14, 18, 16), (11, 19, 10, 18, 27)]

print(mergeHsps(hsps,"AAGTCGGCTAACGCATACGCTTGATAAGCGTA","AGTTGGTTAGAGTATACGCCTGATAAGCGTA"))


# In[13]:


def getBestHsp(seq1,seq2,e,ss,kernel):
    #sequence being searched in database
    kmersSearched = enumerateWords(seq1,len(kernel))
    #sequence in database
    kmersSearching = enumerateWords(seq2,len(kernel))
    
    #chaque tuple de hsps: (searching,searched)/(seq2,seq11)
    hsps = []
    for i in range(len(kmersSearched)):
        matches = findMatches(kmersSearching,kmersSearched[i],kernel)
        
        if(matches != []):
            for j in range(len(matches)):
                hsps.append((matches[j],i))
                
    if(len(hsps) == 0):
        return None
                
    #print(hsps)
         
    extended_hsps = []
    for k in range(len(hsps)):
        extended_hsps.append(extendHsp(seq2,hsps[k][0],seq1,hsps[k][1],kernel,e))
        
    #print(extended_hsps)
            
    extended_hsps_cleaned = list(dict.fromkeys(extended_hsps))
    
    #extended_hsps_cleaned_sorted = sort(extended_hsps_cleaned)
    
    merged_hsps = mergeHsps(extended_hsps_cleaned,seq2,seq1)

    
     
    #return extended_hsps_cleaned
    #return sort(extended_hsps_cleaned)
    
    best = 0
    scoreMax = 0
    for l in range(len(merged_hsps)):
        if(merged_hsps[l][4] > scoreMax):
            scoreMax = merged_hsps[l][4]
            best = l
    
    
    return merged_hsps[best]
    


# In[14]:


print(getBestHsp("AGTTGGTTAGAGTATACGCCTGATAAGCGTA","AAGTCGGCTAACGCATACGCTTGATAAGCGTA",5,5,"11111"))


# In[44]:


# Alignement Smith-Waterman

def traceback(startIndex, matrix, chemins, read1, read2):
    
    tracebackList1 = []
    tracebackList2 = []
    
    i = startIndex[0]
    j = startIndex[1]
    
    while(chemins[i,j] != -1):
        if(chemins[i,j] == 0):
            tracebackList1.append("-")
            tracebackList2.append(read2[i-1])
            i = i-1
        elif(chemins[i,j] == 1):
            tracebackList1.append(read1[j-1])
            tracebackList2.append("-")
            j = j-1
        else:
            tracebackList1.append(read1[j-1])
            tracebackList2.append(read2[i-1])
            i = i-1
            j = j-1
            
    tracebackList1.reverse()
    tracebackList2.reverse()
    
    return(tracebackList1, tracebackList2)

def smithWaterman(seq1, seq2):
    matrix = np.ones((len(seq2)+1, len(seq1)+1))
    matrix[:,0] = 0
    matrix[0] = 0

    chemins = np.zeros((len(seq2)+1, len(seq1)+1))
    chemins = chemins - 1

    for i in range(1, len(seq2)+1):
        for j in range(1, len(seq1)+1):
            diagoValue = -1
            if(seq1[j-1] == seq2[i-1]):
                diagoValue = 1

            matrix[i,j] = max(0, matrix[i-1,j] - 1, matrix[i,j-1] - 1, matrix[i-1,j-1] + diagoValue)
            if(matrix[i,j] == 0):
                chemins[i,j] = -1
            elif(matrix[i,j] == matrix[i-1,j] - 1):
                chemins[i,j] = 0
            elif(matrix[i,j] == matrix[i,j-1] - 1):
                chemins[i,j] = 1
            else:
                chemins[i,j] = 2
    
    scoreMax = matrix.max()
    indice = np.unravel_index(np.argmax(matrix, axis=None), matrix.shape)
                        
    response = traceback(indice, matrix, chemins, seq1, seq2)
    print("response",response)
    # Pourcentage d'identité
    #print(response)
    match = 0
    for i in range (len(response[0])):
        if (response[0][i] == response[1][i]):
            match += 1;
    identity = match / len(response[0]) * 100
    
    return(scoreMax, response, len(response[0]), identity)

print(smithWaterman("CGTAGTCGGCTAACGCATACGCTTGATAAGCGTAAGAGCCC","GCGGATATGATGGAATTGGCAGACATGCTAGACTTAGGATCTAGTGATTATTATCTTGGGGGTTCAAGTCCCTCTATCCGTA"))


# In[45]:


file = "tRNAs.fasta"
sequence = "CGTAGTCGGCTAACGCATACGCTTGATAAGCGTAAGAGCCC"
e = 5
ss = 10
seed = "11111111111"

readsFile = open(file, "r")
readsLines = readsFile.readlines()
readsFile.close()

reads = []
for i in range(0, len(readsLines)):
    reads.append(readsLines[i][:-1])
    
#print(reads)


total = 0

for i in range(0,len(reads),2):
    #retourne (start1,end1,start2,end2,score)
    #print(i)
    bestHsp = getBestHsp(sequence,reads[i+1],e,ss,seed)
    if(bestHsp == None):
        continue
    bitscore = computeBitscore(bestHsp[4])
    evalue = computeEvalue(len(reads[i+1]),len(sequence),bitscore)
    
    if(evalue > ss):
        continue
    
    total += 1
    #retourn (alignment1,alignment2,score,identite)
    print(bestHsp, i+1)
    alignment = smithWaterman(bestHsp,reads[i+1])
    alignment = [1,2,4,5]
    
    print("%s\tscore: %.2f\tident: %.2f" % (reads[i],alignment[0],alignment[3]))
    #print(alignment[1][0])
    #print(alignment[1][1])
    #print(int(bestHsp[2]),int(bestHsp[3]+1))
    print("\n# Best HSP score: %.2f, bitscore: %.2f, evalue: %.2e" % (bestHsp[4],bitscore,evalue))
    print("%d  %s  %d" % (bestHsp[2], sequence[int(bestHsp[2]):int(bestHsp[3])+1], bestHsp[3]))
    print("%d  %s  %d\n\n" % (bestHsp[0], reads[i+1][int(bestHsp[0]):int(bestHsp[1])+1], bestHsp[1]))
    print("----------------------------------------\n")
    
print("Total: "+ str(total))
    


# In[96]:


import argparse
args = []

parser = argparse.ArgumentParser(description='PLAST : Primitive Local Alignment Search Tool')
parser.add_argument('-db', action='store', type=argparse.FileType('r'))
parser.add_argument('-E', action='store', type=int)
parser.add_argument('-ss', action='store', type=int)
parser.add_argument('-seed', action='store')
args = parser.parse_args()
print(args)

