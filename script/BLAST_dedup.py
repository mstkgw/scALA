import argparse

parser = argparse.ArgumentParser()
parser.add_argument("-i",dest="ifasta",action="store",help="input contig",required=True)
parser.add_argument("-b",dest="blastresult",action="store",help="blastresult",required=True)
parser.add_argument("-o",dest="ofasta",action="store",help="dedup contig",required=True)
args = parser.parse_args()

# read input contig
query_dict = {}
for line in open(args.ifasta):
    if line[0:1] == ">":
        seqname = line.rstrip().split(" ")[0][1:]
        query_dict[seqname] = ""
    else:
        query_dict[seqname] += line.rstrip()

# dedup
otxt = []
for line in open(args.blastresult):
    data = line.rstrip().split("\t")
    qseqid, qlen, sseqid, slen, pident, length, qstart, qend, sstart, send = data
    qlen = int(qlen)
    slen = int(slen)
    length = int(length)
    qstart = int(qstart)
    qend = int(qend)
    sstart = int(sstart)
    send = int(send)
    if qseqid not in query_dict: continue
    if qseqid != sseqid: continue
    if length < 6000: continue
    if min(qstart, qend) == min(sstart, send) and max(qstart, qend) == max(sstart, send): continue
    if 1 < min(qstart, qend, sstart, send) < max(qstart, qend, sstart, send) < qlen: continue
    if 1 in [qstart, qend]:
        otxt.append(">" + qseqid + "\n" + query_dict[qseqid][max(qstart, qend):] + "\n")
    elif 1 in [sstart, send]:
        otxt.append(">" + qseqid + "\n" + query_dict[qseqid][max(sstart, send):] + "\n")
    elif qlen in [qstart, qend]:
        otxt.append(">" + qseqid + "\n" + query_dict[qseqid][:min(qstart, qend)] + "\n")
    elif qlen in [sstart, send]:
        otxt.append(">" + qseqid + "\n" + query_dict[qseqid][:min(sstart, send)] + "\n")
    del(query_dict[qseqid])

# output
with open(args.ofasta, "w") as o:
    o.write("".join(otxt))
    for key in sorted(query_dict.keys()):
        o.write(">" + key + "\n" + query_dict[key] + "\n")
