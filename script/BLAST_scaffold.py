import argparse
from Bio.Seq import Seq
from Bio.Alphabet import IUPAC

parser = argparse.ArgumentParser()
parser.add_argument("-b",dest="blastresult",action="store",help="blastresult",required=True)
parser.add_argument("-q",dest="blastquery",action="store",help="blastquery",required=True)
parser.add_argument("-s",dest="blastsubject",action="store",help="blastsubject",required=True)
parser.add_argument("-o",dest="ofile",action="store",help="output scaffold",required=True)
args = parser.parse_args()

def readfasta(fastafile):
    fasta_dict = {}
    for line in open(fastafile):
        if line[0:1] == ">":
            seqname = line.rstrip().split(" ")[0][1:]
            subject_dict[seqname] = ""
        else:
            subject_dict[seqname] += line.rstrip()
    return fasta_dict

subject_dict = readfasta(args.blastsubject)
query_dict = readfasta(args.blastquery)

blast_dict = {}
for line in open(args.blastresult):
    data = line.rstrip().split("\t")
    if data[0] not in blast_dict.keys(): blast_dict[data[0]] = []
    blast_dict[data[0]].append(data)

for key in sorted(blast_dict.keys()):
    if len(blast_dict[key]) == 1 or len(list(set([data[2] for data in blast_dict[key]]))) == 1: continue
    tmpwrite = ""
    for i in range(len(blast_dict[key])-1):
        idata = blast_dict[key][i]
        jdata = blast_dict[key][i+1]
        if idata[2] == jdata[2]: continue
        if int(idata[8]) < int(idata[9]) and int(idata[3])-int(idata[9]) > 1000: flag = True
        elif int(idata[8]) > int(idata[9]) and int(idata[9]) > 1000: flag = True
        elif int(jdata[8]) < int(jdata[9]) and int(jdata[8]) > 1000: flag = True
        elif int(jdata[8]) > int(jdata[9]) and int(jdata[3])-int(jdata[8]) > 1000: flag = True
        elif int(idata[3]) > 10000 and int(idata[5]) < 10000: flag = True
        elif int(jdata[3]) > 10000 and int(jdata[5]) < 10000: flag = True
        else:
            tmpconnect = ""
            if int(idata[8]) < int(idata[9]): tmpconnect += subject_dict[idata[2]][:int(idata[9])]
            else: tmpconnect += str(Seq(subject_dict[idata[2]][int(idata[9])-1:],
                                         IUPAC.ambiguous_dna).reverse_complement())
            if int(idata[7]) <= int(jdata[6]):
                tmpconnect += query_dict[idata[0]][int(idata[7])-1:int(jdata[6])]
                if int(jdata[8]) < int(jdata[9]):
                    tmpconnect += subject_dict[jdata[2]][int(jdata[8])-1:]
                else:
                    tmpconnect += str(Seq(subject_dict[jdata[2]][:int(jdata[8])],
                                          IUPAC.ambiguous_dna).reverse_complement())
            else:
                if int(jdata[8]) < int(jdata[9]):
                    tmpconnect += subject_dict[jdata[2]][int(idata[7])-int(jdata[6])+int(jdata[8]):]
                else:
                    tmpconnect += str(Seq(subject_dict[jdata[2]][:int(jdata[8])-(int(idata[7])-int(jdata[6])+1)],
                                          IUPAC.ambiguous_dna).reverse_complement())
            with open(args.ofile, "w") as o:
                o.write(">"+idata[2]+"\n"+tmpconnect+"\n")
            del(subject_dict[idata[2]])
            del(subject_dict[jdata[2]])
            flag = True
            break
    else:
        flag = False
    if flag == True:
        break

with open(args.ofasta, "w") as o:
    for key in subject_dict.keys():
        o.write(">" + key + "\n" + subject_dict[key] + "\n")
