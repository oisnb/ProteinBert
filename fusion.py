###把序列全部抽出来
import pandas as pd
import os

def positive_protein_csv(input_csv_path_1 , input_csv_path_2, output_csv_path):
    # 读取CSV文件
    df_raw = pd.read_csv(input_csv_path_1)
    df_2 = pd.read_csv(input_csv_path_2)

    # 假设蛋白质序列所在列名为"protein_sequence"，标签所在列名为"label"
    protein_sequences = df_2["ProteinSet"]

    new_data_df = pd.DataFrame(protein_sequences, columns=["ProteinSet"])
    df_out = pd.concat([df_raw, new_data_df], ignore_index=True)
    df_out.to_csv(output_csv_path, index=False)

if __name__ == "__main__":
    # 输入CSV文件路径和输出txt文件路径 将两个afp1_5positive合为一个文件
    input_csv_path = "data/add_ 1_5positive.csv"
    input2_csv_path = "data/positive.csv"
    output_csv_path = "data/afp1_5_positive.csv"
    if os.path.exists(output_csv_path):
        os.remove(output_csv_path)

    positive_protein_csv(input_csv_path, input2_csv_path, output_csv_path)



