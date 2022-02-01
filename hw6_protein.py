"""
Protein Sequencing Project
Name:
Roll Number:
"""

from itertools import count
from random import randrange
from re import A
import re
from typing import final
import hw6_protein_tests as test

project = "Protein" # don't edit this

### WEEK 1 ###

'''
readFile(filename)
#1 [Check6-1]
Parameters: str
Returns: str
'''
def readFile(filename):
    f = open(filename,'r')
    Lines = f.readlines()
    text=''
    for i in Lines:
        text=text+i.strip("\n")
    return text

'''
dnaToRna(dna, startIndex)
#2 [Check6-1]
Parameters: str ; int
Returns: list of strs
'''
def dnaToRna(dna, startIndex):
    condonlist = []
    var = ["UGA","UAG","UAA"]
    for word in range(startIndex,len(dna),3):
        dna = dna.replace("T","U")
        condon = dna[word:word+3]
        if condon not in var:
            condonlist.append(condon)
        else:
            condonlist.append(condon)
            break
    return condonlist


'''
makeCodonDictionary(filename)
#3 [Check6-1]
Parameters: str
Returns: dict mapping strs to strs
'''
def makeCodonDictionary(filename):
    import json

    f = open (filename, "r")
    data = json.loads(f.read())
    keysList = list(data.keys())
    valuesList=list(data.values())
    cdict={}
    for i in range(len(valuesList)):
        for j in range(len(valuesList[i])):

            cdict[valuesList[i][j].replace("T", "U")]=keysList[i]
    return cdict

    

'''
generateProtein(codons, codonD)
#4 [Check6-1]
Parameters: list of strs ; dict mapping strs to strs
Returns: list of strs
'''
def generateProtein(codons, codonD):

    rna=codons
    list=[]
    for keys in rna:
        list.append(codonD[keys])
    if rna[0] == "AUG":
        list[0]="Start"
    return list
    


'''
synthesizeProteins(dnaFilename, codonFilename)
#5 [Check6-1]
Parameters: str ; str
Returns: 2D list of strs
'''
def synthesizeProteins(dnaFilename, codonFilename):
    
    dna=readFile(dnaFilename)
    cd=makeCodonDictionary(codonFilename)
    count=0 ## count unused bases
    final=[]
    k=0
    while k < len(dna):
        sent=dna[k:k+3]
        if sent=="ATG":
            A=dnaToRna(dna,k)
            B=generateProtein(A,cd)
            final.append(B)
            k = k+3*len(A)
        else:
            count=count+1
            k=k+1
    return final


def runWeek1():
    print("Human DNA")
    humanProteins = synthesizeProteins("data/human_p53.txt", "data/codon_table.json")
    print("Elephant DNA")
    elephantProteins = synthesizeProteins("data/elephant_p53.txt", "data/codon_table.json")


### WEEK 2 ###

'''
commonProteins(proteinList1, proteinList2)
#1 [Check6-2]
Parameters: 2D list of strs ; 2D list of strs
Returns: 2D list of strs
'''
def commonProteins(proteinList1, proteinList2):
    common=[]
    for i in range(len(proteinList1)):
        A=proteinList1[i]
        if A in proteinList2:
            if A not in common:
                common.append(A)
    return common


'''
combineProteins(proteinList)
#2 [Check6-2]
Parameters: 2D list of strs
Returns: list of strs
'''
def combineProteins(proteinList):
    proteinList1=sum(proteinList,[])
    return proteinList1


'''
aminoAcidDictionary(aaList)
#3 [Check6-2]
Parameters: list of strs
Returns: dict mapping strs to ints
'''
def aminoAcidDictionary(aaList):
    dict={}
    for word in aaList:
        dict[word]= aaList.count(word)
    return dict


'''
findAminoAcidDifferences(proteinList1, proteinList2, cutoff)
#4 [Check6-2]
Parameters: 2D list of strs ; 2D list of strs ; float
Returns: 2D list of values
'''
def findAminoAcidDifferences(proteinList1, proteinList2, cutoff):
    cmb1 = combineProteins(proteinList1) 
    cmb2 = combineProteins(proteinList2)
    total1 = len(cmb1)
    total2 = len(cmb2)
    prodict1=aminoAcidDictionary(cmb1)
    prodict2=aminoAcidDictionary(cmb2)
    cmb = list(set(cmb1 + cmb2))
    final=[]
    for i in cmb:
        if i != "Start" and i != "Stop":
            if i not in final:
                if i in cmb1 and i in cmb2:
                    freq1 = prodict1[i]/total1
                    freq2 = prodict2[i]/total2
                elif i in cmb1 and i not in cmb2:
                    freq1 = prodict1[i]/total1
                    freq2 = 0
                elif i not in cmb1 and i in cmb2:
                    freq1 = 0 
                    freq2 = prodict2[i]/total2
                freqdiff = abs(freq1-freq2)
                if freqdiff > cutoff:
                    list1=[i,freq1,freq2]
                    final.append(list1)
    return final
    

'''
displayTextResults(commonalities, differences)
#5 [Check6-2]
Parameters: 2D list of strs ; 2D list of values
Returns: None
'''
def displayTextResults(commonalities, differences):
    list = []
    list1 = []
    str = ''
    print(commonalities)
    print("The following proteins occurred in both DNA Sequences:")
    
    for j in commonalities:
        j.remove("Start")
        j.remove("Stop")
        list.append(j)

    for k in list:
        if len(k) > 0:
            k = '-'.join(k)
            list1.append(k)       
    list1.sort()
    for l in list1:
        str += ' '+ l +"\n"
    print(str)
    print("The following amino acids occurred at very different rates in the two DNA sequences:")
    for b in differences:
        wrd = b[0]
        seq1 = round(b[1]*100,2)
        seq2 = round(b[2]*100,2)
        print(f"{wrd}:{seq1}% in Seq1, {seq2}% in seq2")
    return 




def runWeek2():
    humanProteins = synthesizeProteins("data/human_p53.txt", "data/codon_table.json")
    elephantProteins = synthesizeProteins("data/elephant_p53.txt", "data/codon_table.json")

    commonalities = commonProteins(humanProteins, elephantProteins)
    differences = findAminoAcidDifferences(humanProteins, elephantProteins, 0.005)
    displayTextResults(commonalities, differences)


### WEEK 3 ###

'''
makeAminoAcidLabels(proteinList1, proteinList2)
#2 [Hw6]
Parameters: 2D list of strs ; 2D list of strs
Returns: list of strs
'''
def makeAminoAcidLabels(proteinList1, proteinList2):
    return


'''
setupChartData(labels, proteinList)
#3 [Hw6]
Parameters: list of strs ; 2D list of strs
Returns: list of floats
'''
def setupChartData(labels, proteinList):
    return


'''
createChart(xLabels, freqList1, label1, freqList2, label2, edgeList=None)
#4 [Hw6] & #5 [Hw6]
Parameters: list of strs ; list of floats ; str ; list of floats ; str ; [optional] list of strs
Returns: None
'''
def createChart(xLabels, freqList1, label1, freqList2, label2, edgeList=None):
    import matplotlib.pyplot as plt
    return


'''
makeEdgeList(labels, biggestDiffs)
#5 [Hw6]
Parameters: list of strs ; 2D list of values
Returns: list of strs
'''
def makeEdgeList(labels, biggestDiffs):
    return


'''
runFullProgram()
#6 [Hw6]
Parameters: no parameters
Returns: None
'''
def runFullProgram():
    return


### RUN CODE ###

# This code runs the test cases to check your work
if __name__ == "__main__":
    print("\n" + "#"*15 + " WEEK 1 TESTS " +  "#" * 16 + "\n")
    test.week1Tests()
    print("\n" + "#"*15 + " WEEK 1 OUTPUT " + "#" * 15 + "\n")
    runWeek1()

    ## Uncomment these for Week 2 ##
    
    print("\n" + "#"*15 + " WEEK 2 TESTS " +  "#" * 16 + "\n")
    test.week2Tests()
    print("\n" + "#"*15 + " WEEK 2 OUTPUT " + "#" * 15 + "\n")
    runWeek2()
    

    ## Uncomment these for Week 3 ##
    
    print("\n" + "#"*15 + " WEEK 3 TESTS " +  "#" * 16 + "\n")
    test.week3Tests()
    print("\n" + "#"*15 + " WEEK 3 OUTPUT " + "#" * 15 + "\n")
    runFullProgram()
    
