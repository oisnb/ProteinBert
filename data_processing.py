###把序列全部抽出来
import pandas as pd
import os

def positive_protein_csv(input_csv_path, output_csv_path):
    # 读取CSV文件
    df_raw = pd.read_csv(input_csv_path)
    datalist = []

    # 假设蛋白质序列所在列名为"protein_sequence"，标签所在列名为"label"
    protein_sequences = df_raw["ProteinSet"]
    label = 1
    for seq in protein_sequences:
        datalist.append([label, seq])


    df = pd.DataFrame(data=datalist, columns=['label', 'seq'])
    df.to_csv(output_csv_path, index=False)

def negative_protein_csv(input_csv_path, output_csv_path):
    # 读取CSV文件
    df = pd.read_csv(input_csv_path)
    df_out = pd.read_csv(output_csv_path)
    datalist = []

    # 假设蛋白质序列所在列名为"protein_sequence"，标签所在列名为"label"
    protein_sequences = df["ProteinSet"]
    label = 0
    for seq in protein_sequences:
        datalist.append([label, seq])

    new_data_df = pd.DataFrame(datalist, columns=['label', 'seq'])
    df_out = pd.concat([df_out, new_data_df], ignore_index=True)
    df_out.to_csv(output_csv_path, index=False)

if __name__ == "__main__":
    # 输入CSV文件路径和输出txt文件路径
    input_csv_path = "data/afp11_positive.csv"
    input2_csv_path = "data/negative.csv" #最后afp11的负样本要重新抽样
    output_csv_path = "data/afp11_Final.csv"
    if os.path.exists(output_csv_path):
        os.remove(output_csv_path)

    positive_protein_csv(input_csv_path, output_csv_path)

    negative_protein_csv(input2_csv_path, output_csv_path)


