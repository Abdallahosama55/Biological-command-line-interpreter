"""""*******************************************************************************************************************
Name:Abdallah Osama Mohamed
ID:20198053
Name:Amr Hosney Eid
ID:20198059
Name:Amr Alaa Ali
ID:20198060
Name : Hadeel Ali Salim
ID:20198097
******************************** TA: Dina Amr *************** GROUP :B1 ********************BIOPYTHON*****************
"""""

""" Hallo in My Bio python Project this is a Docmentaion to Run this CODE BE Carefull with this Tips  

**************************************************Tips to run this code***********************************************
1-firstly open Cmd 
2-secand : change dir by using  (cd) to project dirctorey  for example  cd E/USER/pycharm_project/biopython
3-thired :write python  main.py  function_Name attbutes  for example  E/USER/pycharm_project/biopython>python main.py gc AGCT

***********************************************************************************************************************
simple example to run and test functions 

1-python main.py gc AGCT
3-python main.py seq_alignment ACTGCC GTCAAG -o out.txt
4-python main.py seq_alignment_files s1.fasta s2.fasta -o out.txt
5-python main.py convert_to_fasta ls-orchid.gbk
6-python main.py transcribe AGCTTACTAG
7*-python main.py filter_nbases AGNNCTTACTAGNN
"""

from Bio.Seq import Seq
import sys, getopt
from Bio import SeqIO
from Bio.Blast import NCBIWWW
from Bio.Blast import NCBIXML
from Bio.Blast import NCBIXML
from Bio import pairwise2
from Bio.pairwise2 import format_alignment

########################################################################################################################
# function to calculate gc precentage in DNA
# #___gc___ Calculate G+C content, return percentage

"""This command takes a seq as a string and returns the gc percentage of it."""
def gc (seq):
   n = seq.count('n') + seq.count('N')
   gc_percentage = float(seq.count('c') + seq.count('C') + seq.count('g')
                  + seq.count('G')) *100/ (len(seq) - n)
   return gc_percentage

######################################################################################################################
# function to find complement of DNA
def complement_DNA(dna):
    N_bases = {'A':'T', 'C':'G','G':'C', 'T':'A','N':'N','a':'t','c':'g','g':'c', 't':'a', 'n':'n'}
    DNA = list(dna)
    DNA = [N_bases [i] for i in DNA]
    return ''.join(DNA)
################################################################################################################
# #__transcribe___The Seq object that defined in the Bio.Seq module combines a transcribe function
def Make_Transcribe(my_dna):
    """This command takes a seq as a string and returns its transcription."""
    my_dna=Seq(my_dna)
    return my_dna.transcribe()


##############################################################################################################

def reverse_string(sequance):
    return sequance[::-1]

def reversecomplement(sequance):
    sequance = reverse_string(sequance)
    sequance = complement_DNA(sequance)
    return sequance
###############################################################################################################
def clc_nbases(DAN_STRING_SEQ):
    DAN_STRING_SEQ=DAN_STRING_SEQ.upper()
    sum_n=DAN_STRING_SEQ.count("N")
    return sum_n

###############################################################################################################
def valid_seq(sequance, type):
    x=False
    """This command takes a seq and a type (protein, dna, rna) and returns a Boolean value of whether it’s a valid
    type or not """
    type = type.lower()
    if type == 'protein':
        sequance = sequance.upper()
        for i in sequance:
            if i not in 'ABCDEFGHIKLMNPQRSTVWXYZ':
                x= False
            else:
                x=True


        if x==True:
            return "it's a valid protein"
        else:
            return "it's not valid protein"
######################################################################################################################
    elif type == 'dna':
        flag_dna = False
        sequance = sequance.upper()
        for i in sequance:
            if i not in 'AGCT':
                flag_dna= False
            else:
                flag_dna = True

        if flag_dna==True:
            return "it's a valid Dna"
        else:
            return "it's not valid Dna"
############@################################@#############################@###############$#################@######
    elif type == 'rna':
        flag_RNA = False
        sequance=sequance.upper()
        for i in sequance:
            if i not in 'AGCU':
                flag_RNA= False
            else:
                flag_RNA = True
        if flag_RNA==True:
            return "it's a valid Rna"
        else:
            return "it's not valid Rna"
    else:
        return 'Invalid Type or invalid Sequance'
#########################################################################################################################

def filter_nbases(seqaunce):
    """This command takes a seq and returns the Seq after removing n bases."""

    for i in seqaunce:
        if i not in 'ATGCN' or i not in 'atgcn':
            return 'THere is wrong seq'
    seqaunce = seqaunce.replace("N", "")
    return seqaunce

#########################################################################################################################
def seq_alignment(seqaunce_one,sequance_2,out=""):
    seq1 = seqaunce_one.upper()
    seq2 = sequance_2.upper()

    alignmet = []
    for i in seq1:
        if i not in 'AGCT':
            print('this is wrong sequace')
            return
    for j in seq2:
        if j not in 'AGCT':
            print('this is wrong sequace')
            return
    alignments = pairwise2.align.globalxx(seq1, seq2)  # global alignment
    for alignment in alignments:
        if out =="":
            print(alignment)
            print(format_alignment(*alignment))
            return 0
        else:
            file_output = open(out, 'w')
            for alignment in alignments:
                nonFormattedAlignment = str(alignment)
                file_output.write(nonFormattedAlignment)
                file_output.write('\n')
                formattedAlignment = str(format_alignment(*alignment))
                file_output.write(formattedAlignment)
                file_output.write('\n')
            file_output.close()
            print("output File done seccessfuly ")

'''''''''if we want  seq qlignment soultion from scratch
    score = 0
    gap = -1
    match = 3
    mismatch = -2
    if(len(seq1)<= len(seq2)):
        for i in range(len(seq1)):
            if (seq1[i]==seq2[i]):
                score=score+match
                alignmet.append(seq1[i])
            elif(seq1[i]==" "):
                score=score+gap
            elif(seq1[i]=="N"):
                score=score+gap
            elif(seq1[i]=="n"):
                score=score+gap
            else:
                score=score+mismatch
        if (out!=""):
            file=open(out,"w")
            file.write("".join(alignmet))
            file.close()
            print("file done")
        return " the score is : " + str(score) +  "  the alignment is : " + "".join(alignmet)
    else: print("please enter seq1 and seq2 in same range ")
'''''''''''
#########################################################################################################################
def seq_alignment_files(f_align1, f2_align, outputfile=""):
    """This command takes 2 fasta files as input, each file contains a single sequence. It reads the 2 sequences from
    files and get all their alignments along with the score. The -o is an optional parameter if we need the output to
    be written on a file instead of the screen. """
    try:
        sequcance1 = SeqIO.read(f_align1, 'fasta')
        sequcance2 = SeqIO.read(f2_align, 'fasta')
    except OSError as Error:
        print(Error)
        return 'this wrong file name ,please enter anothe'
    alignments = pairwise2.align.globalxx(sequcance1, sequcance2)  # global alignment
    if outputfile == '':
        for align in alignments:
            print(align)
            print(format_alignment(*align))
    else:

        f = open(outputfile, 'w')
        for alignment in alignments:
            nonFormattedAlignment = str(alignment)
            f.write(nonFormattedAlignment)
            f.write('\n')
            formattedAlignment = str(format_alignment(*alignment))
            f.write(formattedAlignment)
            f.write('\n')
        f.close()

    print('Alignmnet Done to File ', outputfile)

########################################################################################################################
########################################################################################################################
def online_alignment(seq,filepatth=""):
    seq=seq.upper()
    fastaSequence = ">sequence_unknown\n"+seq
    blast_record = NCBIWWW.qblast("blastn", "nr", sequence=fastaSequence)
    blast_record = NCBIXML.read(blast_record)
    if(len(str(filepatth))>0):
        file=open(filepatth,"wt")
        for alignment in blast_record.alignments:
            for hsp in alignment.hsps:
                file.write('*Alignment* \n')
                file.write('sequence: '+ str(alignment.title)+"\n")
                file.write('length:'+str(alignment.length)+"\n")
                file.write('e value:'+ str(hsp.expect)+"\n")
                file.write(str(hsp.query)+"\n")
                file.write(str(hsp.match)+"\n")
                file.write(str(hsp.sbjct)+"\n")
            file.close()
    else:
        for alignment in blast_record.alignments:
            for hsp in alignment.hsps:
                print('*Alignment*')
                print('sequence:', alignment.title)
                print('length:', alignment.length)
                print('e value:', hsp.expect)
                print(hsp.query)
                print(hsp.match)
                print(hsp.sbjct)

#print(online_alignment("CATGCTACGGTGCTAAAAGCATTACGCCCTATAGTGATTTTCGAGACATACTGTGTTTTTAAATATAGTATTGCC",filepatth="out.txt"))

##########################################################################################################################

def merge_fasta(file1, file2,output='',*files):
    """This command takes any number of fasta files (at least two) and merge their contents into one fasta output
    file. There is an option to write the merge result in a file using -o option, otherwise the merge result will be
    displayed on the console. """
    try:
        file1 = SeqIO.parse(file1, 'fasta')
        file2 = SeqIO.parse(file2, 'fasta')
    except OSError as Error:
        print(Error)
        return
    FilesList = []
    FilesList.append(file1)
    FilesList.append(file2)
    for file in files:
        try:
            file = SeqIO.parse(file, 'fasta')
        except OSError as Error:
            print(Error)
            return
        FilesList.append(file)

    if output == '':
        for file in FilesList:
            for record in file:
                print(record.id)
                print(record.description)
                print(record.seq)
    else:

        Merged_file=open(str(output),'wt')
        for file in FilesList:
            for record in file:
                Merged_file.write(str(record.id+"\n"))
                Merged_file.write(str(record.description+"\n"))
                Merged_file.write(str(record.seq+"\n"))



#print(merge_fasta("s1.fasta","s2.fasta","output.fasta","s3.fasta",))
################################################ # convert_to_fasta file########################################################################
def tranform_to_fasta(file_genebank):
    # # convert_to_fasta file
    """This command converts the input genbank file with multiple records onto a fasta formatted file. The output is
    to be written in a different output fasta file. """
    if file_genebank[-3:] == 'gbk':
        output = file_genebank[:-3] + 'fasta'
    try:
        with open(file_genebank, "r") as input:
                with open(output, "w") as output:
                    my_sequences = SeqIO.parse(input, "genbank")  # printing the records
                    # count Record iterator to the write function that converte to fasta file
                    c_iterator = SeqIO.write(my_sequences, output, "fasta")
        print("There are %i records converted into a fasta file" % c_iterator)
    except OSError as Error:
            print(Error)
            return
    else:
        print('File must be genbank\n',tranform_to_fasta.__doc__)

##########################################################################################################################
"""""EROOR handling passing function with getopt attrbutes """
def get_opt_function():
    print('gc: \n', gc.__doc__, '\n')
    print('transcribe:\n', Make_Transcribe.__doc__, '\n')
    print('reverse_complement:\n', reversecomplement.__doc__, '\n')
    print('calc_nbases:\n', clc_nbases.__doc__, '\n')
    print('is_valid:\n', valid_seq.__doc__, '\n')
    print('filter_nbases:\n', filter_nbases.__doc__, '\n')
    print('seq_alignment:\n', seq_alignment.__doc__, '\n')
    print('seq_alignment_files:\n', seq_alignment_files.__doc__, '\n')
    print('online_alignment:\n', online_alignment.__doc__, '\n')
    print('merge_fasta:\n', merge_fasta.__doc__, '\n')
    print('convert_to_fasta:\n', tranform_to_fasta.__doc__, '\n')


try:
    Loption, Arg_value = getopt.gnu_getopt(sys.argv[1:], 'o:h', 'help')
except getopt.GetoptError as error:
    print(error)
    sys.exit()
if Arg_value != []:
    if Arg_value[0] == 'gc':
        if not Loption:
            if len(Arg_value) == 2:
                print(gc(Arg_value[1]))
            else:
                print('WRONG IN NUMBER OF PAPRMETER :\n', gc.__doc__)
        else:
            print('this is invalid option in command %s\n' % Loption[0][0], gc.__doc__)
    elif Arg_value[0] == 'transcribe':
        if not Loption:
            if len(Arg_value) == 2:
                print(Make_Transcribe(Arg_value[1]))
            else:
                print('WRONG IN NUMBER OF PAPRMETER :\n', Make_Transcribe.__doc__)
        else:
            print('this is invalid option in command %s\n' % Loption[0][0], Make_Transcribe.__doc__)
    elif Arg_value[0] == 'reverse_complement':
        if not Loption:
            if len(Arg_value) == 2:
                print(reversecomplement(Arg_value[1]))
            else:
                print('WRONG IN NUMBER OF PAPRMETER :\n', reversecomplement.__doc__)
        else:
            print('this is invalid option in command %s\n' % Loption[0][0], reversecomplement.__doc__)
    elif Arg_value[0] == 'calc_nbases':
        if not Loption:
            if len(Arg_value) == 2:
                print(clc_nbases(Arg_value[1]))
            else:
                print('WRONG IN NUMBER OF PAPRMETER :\n', clc_nbases.__doc__)
        else:
            print('this is invalid option in command %s\n' % Loption[0][0], clc_nbases.__doc__)
    elif Arg_value[0] == 'is_valid':
        if not Loption:
            if len(Arg_value) == 3:
                print(valid_seq(Arg_value[1], Arg_value[2]))
            else:
                print('WRONG IN NUMBER OF PAPRMETER :\n', valid_seq.__doc__)
        else:
            print('this is invalid option in command %s\n' % Loption[0][0], valid_seq.__doc__)

    elif Arg_value[0] == 'filter_nbases':
        if not Loption:
            if len(Arg_value) == 2:
                print(filter_nbases(Arg_value[1]))
            else:
                print('WRONG IN NUMBER OF PAPRMETER :\n', filter_nbases.__doc__)
        print('this is invalid option in command %s\n' % Loption[0][0], filter_nbases.__doc__)

    elif Arg_value[0] == 'seq_alignment':
        if not Loption:
            if len(Arg_value) == 3:
                seq_alignment(Arg_value[1], Arg_value[2])
            else:
                print('Less or more number of prammeters:\n', seq_alignment.__doc__)
        else:
            if '-o' in Loption[0]:
                if len(Arg_value) == 3:
                    seq_alignment(Arg_value[1], Arg_value[2], Loption[0][1])
                else:
                    print('WRONG IN NUMBER OF PAPRMETER :\n', seq_alignment.__doc__)
            else:
                print('Wrong option %s not recognized' % Loption[0][0])
    elif Arg_value[0] == 'seq_alignment_files':
        if not Loption:
            if len(Arg_value) == 3:
                seq_alignment_files(Arg_value[1], Arg_value[2])
            else:
                print('WRONG IN NUMBER OF PAPRMETER :\n', seq_alignment_files.__doc__)
        else:
            if '-o' in Loption[0]:
                if len(Arg_value) == 3:
                    seq_alignment_files(Arg_value[1], Arg_value[2], Loption[0][1])
                else:
                    print('WRONG IN NUMBER OF PAPRMETER :\n', seq_alignment_files.__doc__)
            else:
                print('Wrong option %s not recognized' % Loption[0])
    elif Arg_value[0] == 'online_alignment':
        if not Loption:
            if len(Arg_value) == 2:
                online_alignment(Arg_value[1])
            else:
                print('WRONG IN NUMBER OF PAPRMETER :\n', online_alignment.__doc__)
        else:
            if '-o' in Loption[0]:
                if len(Arg_value) == 2:
                    online_alignment(Arg_value[1], Loption[0][1])
                    print('Alignmnet Done to File ', Loption[0][1])
                else:
                    print('WRONG IN NUMBER OF PAPRMETER :\n', online_alignment.__doc__)
            else:
                print('Wrong option %s not recognized' % Loption[0])
    elif Arg_value[0] == 'merge_fasta':
        if not Loption:
            if len(Arg_value) == 3:
                merge_fasta(Arg_value[1], Arg_value[2])

            elif len(Arg_value) > 3 and not Loption:
                files = []
                for i in range(3, len(Arg_value)):
                 files.append(Arg_value[i])
                 print(files)

                merge_fasta(Arg_value[1], Arg_value[2],files)
            else:
                print('less argumnents for merge_fasta\n', merge_fasta.__doc__)
        else:

            if len(Arg_value) == 3 and '-o' in Loption[0]:
                merge_fasta(Arg_value[1], Arg_value[2], Loption[0][1])
            elif len(Arg_value) > 3 and '-o' in Loption[0]:
                files = []
                for i in range(3, len(Arg_value)):
                    files.append(Arg_value[i])
                files = tuple(files)
                merge_fasta(Arg_value[1], Arg_value[2], Loption[0][1],files)
            else:
                print('less argumnents for merge_fasta or Wrong option\n', merge_fasta.__doc__)
    elif Arg_value[0] == 'convert_to_fasta':
        if not Loption:
            if len(Arg_value) == 2:
                tranform_to_fasta(Arg_value[1])
            else:
                print('Less or more number of prammeters:\n',tranform_to_fasta.__doc__)
        else:
            print('option is not valid in this command %s\n' % Loption[0][0], tranform_to_fasta.__doc__)

    else:
        print('This command is wrong')
        get_opt_function()

elif Loption:
    if '-h' in Loption[0] or '--help' in Loption[0]:
        get_opt_function()

##############################################################################################################################

