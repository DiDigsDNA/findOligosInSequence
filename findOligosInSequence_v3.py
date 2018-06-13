#!/usr/bin/env python

'''Searches nucleotide fasta sequences for exact matches to oligonucleotides (i.e. primers, probes) or their
reverse complements. Accepts degenerate bases in oligonucleotides. Outputs: (i)sequences requiring analytical
follow-up (i.e. in DNA analysis software) to a fasta file,(ii) sequence id's (usually genbank accession numbers)
and  details of which oligos were not found to an information .txt file, and (iii) sequence id's and found/
found information on oligos to a csv file which can be viewed and sorted in spreadsheet software'''

'''Author: Diane Eisler, Molecular Microbiology & Genomics, BCCDC Public Health Laboratory,Sept 2017'''
import sys,string,os, time, Bio
from Bio import Seq, SeqIO, SeqUtils, Alphabet, SeqRecord
from Bio.Alphabet import IUPAC
from Bio.Seq import Seq

fastaToParse = sys.argv[1] #fasta file of nucleotide sequences to parse
oligosToFind = sys.argv[2] #fasta file of primer/probe sequences to search for
gb_accessionsHandle = sys.argv[3]  #user-specified output file name
gb_access_csv = open(gb_accessionsHandle + ".csv",'w') # csv list of GB access #'s and found/unfound state of oligos
gb_accessions_txt= open(gb_accessionsHandle + ".txt",'w') #output file of GB access #'s and details re: which oligos NOT found
outputFastaHandle = sys.argv[3] + ".fasta" #from user=specified output name
outputFasta = open(outputFastaHandle, 'w') #output fasta containing sequences requiring further analysis
console = sys.stderr #user info printed to the screen

def configureOutputFile():
    '''Writes run information to the text output file.'''
    localtime = time.asctime(time.localtime(time.time())) #date and time of analysis
    gb_accessions_txt.write("---------------------------------------------------------------------------\n")
    gb_accessions_txt.write("RESULTS OF SEARCH FOR EXACT MATCHES TO OLIGONUCLEOTIDES IN QUERY SEQUENCES:\n")
    gb_accessions_txt.write("---------------------------------------------------------------------------\n\n")
    gb_accessions_txt.write("Sequences Searched: %s \n\n" % (fastaToParse))
    gb_accessions_txt.write("Oligonucleotides Searched for: %s \n\n" % (oligosToFind))
    gb_accessions_txt.write("Analysis as of: " + localtime + "\n\n")
    gb_accessions_txt.write("The following Genbank accession numbers require further analysis in Geneious:\n")
    gb_accessions_txt.write("(Note: only sequences with primer/probe mismatches are listed.)")
    gb_accessions_txt.write("\n\n---------------------------------------------------------------------------")
    #if none of the query genomes require further analysis, print a message to outfile
    return
    
def printConsoleOutput():
    '''Prints user output to console.'''
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
    '''Writes headers (i.e. names of oligos) to csv output file.'''
 	
    gb_access_csv.write("Sequence_ID")
 	#write the names of all oligos as column headers
    for i, oligo in enumerate(oligoList):
        gb_access_csv.write(',' + oligo.id)
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
    print the SeqRecord id to console, csv and text output files, followed by info on unfound oligos.'''
    console.write("______________________________________________________\n")
    console.write ("%s:\n" % (record.id))
    #print search results to information file
    gb_accessions_txt.write("\n\n%s:" % (record.id)) #write sequence name to txt file
    gb_access_csv.write("\n%s" % (record.id)) #write sequence name to csv file
    for oligo_name in oligoFoundDict: #write FOUND/NOT_FOUND to csv file for each oligo
        if oligoFoundDict[oligo_name] == False:
            gb_accessions_txt.write("\n'%s' NOT found" % (oligo_name)) #write unfound oligo name to text file
            gb_access_csv.write(",NOT_FOUND")
            console.write(oligo_name + " NOT found\n") #write "Not found" to console
        else: #insert a comma into the csv file
        	gb_access_csv.write(",found")
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
printConsoleOutput() #print run data to console

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


print("Search Results text file : %s.txt" % gb_accessionsHandle) #print names of info and fasta output files
print("Search Results csv file: %s.csv" % gb_accessionsHandle) #print names of info and fasta output files
print ("Fasta sequences for further analysis: " + outputFastaHandle + "\n")
inFile.close()
gb_access_csv.close()
gb_access_csv.close()
