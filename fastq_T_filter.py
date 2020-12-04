import sys
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

def main(fastq_in):

    f =  open(fastq_in)

    output =  open(".".join([fastq_in, "Ubias"]) , "w")

    for record in SeqIO.parse(f, "fastq"):

        if record.seq[0]=="T":

            fastq_out = SeqRecord(record.seq, id=record.id, description="")
            fastq_out.letter_annotations["phred_quality"] = record.letter_annotations["phred_quality"]
            output.write(fastq_out.format("fastq"))


if __name__ == "__main__":
    main(sys.argv[1])
