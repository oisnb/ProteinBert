import pandas as pd
import random

def split_data(input_csv_path, train_csv_path, test_csv_path, test_ratio=0.33):
    # 读取CSV文件
    df = pd.read_csv(input_csv_path)

    # 计算测试集的大小
    test_size = int(len(df) * test_ratio)

    # 随机打乱数据
    df = df.sample(frac=1, random_state=42)

    # 划分数据为测试集和训练集
    test_data = df[:test_size]
    train_data = df[test_size:]

    # 将测试集和训练集保存为CSV文件
    test_data.to_csv(test_csv_path, index=False)
    train_data.to_csv(train_csv_path, index=False)

if __name__ == "__main__":
    # 输入CSV文件路径和输出训练集、测试集的CSV文件路径
    input_csv_path = "data/afp1_5_Final_10000.csv"
    train_csv_path = "protein_benchmarks/igem_afp1_5_train_data_update_10000.csv"
    test_csv_path = "protein_benchmarks/igem_afp1_5_test_data_update_10000.csv"

    # 设置测试集占总数据的比例，默认为1/3
    test_ratio = 0.33

    split_data(input_csv_path, train_csv_path, test_csv_path, test_ratio)
