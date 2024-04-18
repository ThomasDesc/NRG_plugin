import argparse
import pymol

#python process_result_flexaid.py -f RESULT_0.pdb -o RESULT_0_processed.pdb

parser = argparse.ArgumentParser(description="Your script description here")
parser.add_argument("-f", type=str, required=True, help="Path to the PDB result file")
parser.add_argument("-o", type=str, required=True, help="output_filename")
args = parser.parse_args()

file = args.f
output = args.o

with open(file, 'r') as t1:
    with open(output,'w') as t2:
        text=t1.readlines()
        for line in text:
            if 'REMARK' not in line:
                if 'LIG  9999' in line:
                    a_name=line[12:17].split()[0]+line[9:11]+' '*(5-len(line[12:17].split()[0]+line[9:11]))
                    new_line=line[:12]+a_name+line[17:21]+'L'+line[22:]
                    t2.write(new_line)
                else:
                    t2.write(line)
