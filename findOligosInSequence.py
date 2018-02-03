'''Searches sequences in fasta format for exact matches to oligonucleotides (i.e. primers, probes) or their
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
gb_accessionsHandle = sys.argv[3] + ".txt" #user-specified output file name
gb_accessions= open(gb_accessionsHandle,'w') #output file of GB access #'s and details re: which oligos NOT found
outputFastaHandle = sys.argv[3] + ".fasta" #from user=specified output name
outputFasta = open(outputFastaHandle, 'w') #output fasta containing sequences requiring further analysis

def configureOutputFile():
    localtime = time.asctime(time.localtime(time.time())) #date and time of analysis
    gb_accessions.write("---------------------------------------------------------------------------\n")
    gb_accessions.write("RESULTS OF SEARCH FOR EXACT MATCHES TO OLIGONUCLEOTIDES IN QUERY SEQUENCES:\n")
    gb_accessions.write("---------------------------------------------------------------------------\n\n")
    gb_accessions.write("Sequences Searched: %s \n\n" % (fastaToParse))
    gb_accessions.write("Oligonucleotides Searched for: %s \n\n" % (oligosToFind))
    gb_accessions.write("Analysis as of: " + localtime + "\n\n")
    gb_accessions.write("The following Genbank accession numbers require further analysis in Geneious:\n")
    gb_accessions.write("(Note: only sequences with primer/probe mismatches are listed.)")
    gb_accessions.write("\n\n---------------------------------------------------------------------------")
    #if none of the query genomes require further analysis, print a message to outfile
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
    print "______________________________________________________"
    print "%s:" % (record.id)
    #print search results to information file
    #gb_accessions.write("\n______________________________________________________\n%s:" % (record.id))
    gb_accessions.write("\n\n%s:" % (record.id))
    #gb_accessions.write("\n%s:" % (record.id))
    for oligo in oligoFoundDict:
        if oligoFoundDict[oligo] == False:
            #print SeqRecord id of oligo
            gb_accessions.write("\n'%s' NOT found" % (oligo.id))
            print oligo.id + " NOT found"
    return

def searchSequenceForOligoSet(record):
    '''Searches a SeqRecord for each oligo in the oligo list (and for its reverse complement). '''
    oligoFoundDict = {} #dictionaries to store found/not found state of oligos and revcomp's
    #rc_oligoDict = {}
    oligoSetResults = [] #empty list to hold search result tuples
    for i, oligo in enumerate(oligoList): #oligos in oligoList and rc_oligoList already in same order
        oligoFound = parseSeqRecordForOligo(record, oligo) #search for oligo in sequence
        rc_oligoFound = parseSeqRecordForOligo(record, rc_oligoList[i]) #search for oligo rev comp
        oligoResults = (oligoFound,rc_oligoFound) #store both results in a tuple
        oligoSetResults.append(oligoResults) #DO I REALLY NEED TO STORE THIS?
    #I ONLY WNAT TO PRINT THE GB ACCESS TO TEXT FILE IF THERE IS AN UNFOUND PRIMER
        if (oligoFound == False) and (rc_oligoFound == False):
            oligoFoundDict[oligo]=False
        else:
            oligoFoundDict[oligo]=True
    if False in oligoFoundDict.values(): 
        printSearchResults(oligoFoundDict,record) #print SeqRecord id (GB acc #) and unfound oligos to file
        SeqIO.write(record, outputFasta, "fasta") #output SeqRecord to fasta file
    return

#main code starts here...
configureOutputFile() #write headers to output text file

with open(oligosToFind,'r') as oligoFile:
    #read in primers from fasta file as Sequence objects, uppercase, add to oligoList
    oligoList = [oligo.upper() for oligo in SeqIO.parse(oligoFile, "fasta", alphabet=IUPAC.ambiguous_dna)]
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
    print "\nSearching for oligonucleotides in %i sequences..." % len(seqList)
    for record in seqList: #parse SeqRecord for exact match to each oligo and reverse complement
        searchSequenceForOligoSet(record)
print "\nEND OF RESULTS\n"
gb_accessions.write("\n------------------------------------------------------")
gb_accessions.write("\n\nEND OF RESULTS")
gb_accessions.write("\nOutput file: " + gb_accessionsHandle)
gb_accessions.write("\nFasta sequences for analysis: " + outputFastaHandle)
print("Search Results Output file: " + gb_accessionsHandle) #print names of info and fasta output files
print ("Fasta sequences for further analysis: " + outputFastaHandle + "\n")
inFile.close()
gb_accessions.close()
