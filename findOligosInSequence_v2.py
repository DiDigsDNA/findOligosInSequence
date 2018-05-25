#!/usr/bin/env python

'''Searches nucleotide fasta sequences for exact matches to oligonucleotides (i.e. primers, probes) or their
reverse complements. Accepts degenerate bases in oligonucleotides. Outputs sequences requiring analytial
follow-up (i.e. in DNA analysis software) to a fasta file. Outputs genbank accession numbers, along with
details of primers not found to an information .txt file.'''

'''Author: Diane Eisler, Molecular Microbiology & Genomics, BCCDC Public Health Laboratory,Sept 2017'''
import sys,string,os, time, Bio
from Bio import Seq, SeqIO, SeqUtils, Alphabet, SeqRecord
from Bio.Alphabet import IUPAC
from Bio.Seq import Seq

fastaToParse = sys.argv[1] #fasta file of nucleotide sequences to parse
oligosToFind = sys.argv[2] #fasta file of primer/probe sequences to search for
gb_accessionsHandle = sys.argv[3] + ".csv" #user-specified output file name
gb_accessions= open(gb_accessionsHandle,'w') #output file of GB access #'s and details re: which oligos NOT found
outputFastaHandle = sys.argv[3] + ".fasta" #from user=specified output name
outputFasta = open(outputFastaHandle, 'w') #output fasta containing sequences requiring further analysis
console = sys.stderr #user info printed to the screen

def printUserOutput():
    localtime = time.asctime(time.localtime(time.time())) #date and time of analysis
    console.write("---------------------------------------------------------------------------\n")
    console.write("RESULTS OF SEARCH FOR EXACT MATCHES TO OLIGONUCLEOTIDES IN QUERY SEQUENCES:\n")
    console.write("---------------------------------------------------------------------------\n")
    console.write("Sequences Searched: %s \n" % (fastaToParse))
    console.write("Oligonucleotides Searched for: %s \n" % (oligosToFind))
    console.write("Analysis as of: " + localtime + "\n")
    console.write("IDs of sequences with primer/probe mismatches: %s" % (gb_accessionsHandle))
    console.write("\nFasta sequences written to: " + outputFastaHandle + "\n\n")
    #console.write("\n\n---------------------------------------------------------------------------\n")
    return

def initializeCSVOutputFile(oligoList):
 	gb_accessions.write("Sequence_ID")
 	#write the names of all oligos as column headers
 	for i, oligo in enumerate(oligoList):
 		gb_accessions.write(',' + oligo.id)
 	return 

def parseSeqRecordForOligo(record,oligo):
    '''Parse SeqRecord for oligo and return True if found and False if not.'''
    results = SeqUtils.nt_search(str(record.seq),oligo) #search in SeqRecord sequence for oligo
    if (len(results) > 1):
        return True #if list > 1 item, a match position was found
    else: #print "Did NOT find %s in %s" % (ol.id, record.id)
        return False

def printSearchResults(oligoFoundDict, record):
    '''Check oligoFoundDict for oligos NOT found in record and, in case of unfound oligos,
    print the SeqRecord id to consol and information file, followed by names of unfound oligos.'''
    #print search results to console
    console.write("______________________________________________________\n")
    console.write ("%s:\n" % (record.id))
    #print search results to information file
    gb_accessions.write("\n%s" % (record.id)) #write sequence name in csv file
    for oligo_name in oligoFoundDict:
        if oligoFoundDict[oligo_name] == False:
            #print ',' followed by oligo name
            #gb_accessions.write(",%s" % (oligo_name))
            gb_accessions.write(",NOT_FOUND")
            #print info to console
            console.write(oligo_name + " NOT found\n")
        else: #insert a comma into the csv file
        	gb_accessions.write(",found")
    return

def searchSequenceForOligoSet(record, oligoList, rc_oligoList):
    '''Searches a SeqRecord for each oligo in the oligo list (and for its reverse complement). '''
    oligoFoundDict = {} #dictionaries to store found/not found state of oligos and revcomp's
    oligoSetResults = [] #empty list to hold search result tuples
    for i, oligo in enumerate(oligoList): #oligos in oligoList and rc_oligoList already in same order
        oligoFound = parseSeqRecordForOligo(record, oligo) #search for oligo in sequence
        rc_oligoFound = parseSeqRecordForOligo(record, rc_oligoList[i]) #search for oligo rev comp
    	#print GB accession # to file only if there is an unfound primer/probe
        if (oligoFound == False) and (rc_oligoFound == False):
        	oligoFoundDict[oligo.id]=False
        else:
            oligoFoundDict[oligo.id]=True
    if False in oligoFoundDict.values(): 
        printSearchResults(oligoFoundDict,record) #print SeqRecord id (GB acc #) and unfound oligos to file
        SeqIO.write(record, outputFasta, "fasta") #output SeqRecord to fasta file
    return


#main code starts here...
printUserOutput() #write headers to console

with open(oligosToFind,'r') as oligoFile:
    #read in primers from fasta file as Sequence objects, uppercase, add to oligoList
    oligoList = [oligo.upper() for oligo in SeqIO.parse(oligoFile, "fasta", alphabet=IUPAC.ambiguous_dna)]
    initializeCSVOutputFile(oligoList) #print headers to csv output file
    rc_oligoList = [] # reverse complements of oligos in oligoList
    for oligo in oligoList: #reverse complement oligo and add to rc_oligoList
        rev_comp_oligo = oligo.reverse_complement()
        rev_comp_oligo.id = oligo.id + "_RC"
        rc_oligoList.append(rev_comp_oligo)
    for i in range(len(oligoList)):
        print("oligo %s: %s\treverse complement: %s (%s)" % (oligoList[i].id, oligoList[i].seq, rc_oligoList[i].seq,
            rc_oligoList[i].id))

with open(fastaToParse,'r') as inFile:
    #read nucleotide sequences from fasta file into SeqRecords, uppercases and adds to a list
    seqList = [rec.upper() for rec in list(SeqIO.parse(inFile, "fasta", alphabet=IUPAC.ambiguous_dna))]
    print("\nSearching for oligonucleotides in %i sequences..." % len(seqList))
    for record in seqList: #parse SeqRecord for exact match to each oligo and reverse complement
        searchSequenceForOligoSet(record, oligoList, rc_oligoList)

console.write("\n------------------------------------------------------\n")
console.write("END OF RESULTS\n")
console.write("------------------------------------------------------\n")


print("Search Results located at: " + gb_accessionsHandle) #print names of info and fasta output files
print ("Fasta sequences for further analysis: " + outputFastaHandle + "\n")
inFile.close()
gb_accessions.close()
