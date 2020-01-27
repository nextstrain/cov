"""
Mask initial bases from alignment FASTA
"""
import argparse
import Bio
import Bio.SeqIO
from Bio.Seq import Seq


if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description="Mask initial bases from alignment FASTA",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    parser.add_argument("--alignment", required=True, help="FASTA file of alignment")
    parser.add_argument("--mask-from-beginning", type = int, help="number of bases to mask from start")
    parser.add_argument("--mask-from-end", type = int, help="number of bases to mask from end")
    parser.add_argument("--mask-sites", nargs='+', type = int,  help="list of sites to mask")
    parser.add_argument("--mask-from-X-to-X", nargs="+", type=int, help="start and end bp of masking")
    parser.add_argument("--output", required=True, help="FASTA file of output alignment")
    args = parser.parse_args()

    begin_length = 0
    if args.mask_from_beginning:
        begin_length = args.mask_from_beginning
    end_length = 0
    if args.mask_from_end:
        end_length = args.mask_from_end
    mask_from_to = [0,0]
    if args.mask_from_X_to_X:
        mask_from_to = args.mask_from_X_to_X

    with open(args.output, 'w') as outfile:
        for record in Bio.SeqIO.parse(args.alignment, 'fasta'):
            seq = str(record.seq)
            start = "N" * begin_length
            middle = seq[begin_length:-end_length] if end_length != 0 else seq[begin_length:len(seq)]
            end = "N" * end_length
            seq_list = list(start + middle + end)
            if args.mask_sites:
                for site in args.mask_sites:
                    seq_list[site-1] = "N"
            if args.mask_from_X_to_X:
                Ns = list("N" * (mask_from_to[1]-mask_from_to[0]))
                seq_list[mask_from_to[0]:mask_from_to[1]] = Ns
            record.seq = Seq("".join(seq_list))
            Bio.SeqIO.write(record, outfile, 'fasta')
