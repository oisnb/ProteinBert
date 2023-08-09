import pandas as pd
from Bio.Seq import Seq
from Bio.SeqUtils.IsoelectricPoint import IsoelectricPoint as IP

input_path = "data/Gene_11.txt"
output_path = "data/"
mark = "DNA"

def loadData():
    datalist = []
    ProteinSet = []
    with open(input_path, 'r') as f:
        for line in f:
            if mark in line:
                position = line.find(":")
                gene_seq = Seq(line[position+2:])
                # 将DNA序列转化为蛋白质序列；蛋白质序列中最后的*代表终止，可以删除
                protein_seq = gene_seq.translate()
                print(protein_seq)
                position = protein_seq.find("*")
                protein_seq = protein_seq[:position]
                ProteinSet.append(protein_seq)
                datalist.append([gene_seq, protein_seq])

    return datalist


if __name__ == '__main__':
    datalist = loadData()
    df = pd.DataFrame(data=datalist, columns=['DNASet', 'ProteinSet'])
    df.to_csv(output_path + 'afp11_positive.csv', index=False)
