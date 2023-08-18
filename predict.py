import pickle
import time
import os
import pandas as pd
from proteinbert import OutputType, OutputSpec, load_pretrained_model
from proteinbert.finetuning import predict_our
from Bio import SeqIO


# 完整的预测函数
def prediction(seqences_data, output_path, model, threshold):
    # seqences_data 为要预测的序列  output_path 为暂时保存的路径，可以不管   model 为afp1_5或afp11模型所在地   threshold为阈值
    if os.path.exists(output_path + '.csv'):
        os.remove(output_path + '.csv')
    df = pd.DataFrame(data=seqences_data, columns=['seq'])
    df.to_csv(output_path + '.csv', index=False)

    test_set_file_path = os.path.join(output_path + '.csv')
    test_set = pd.read_csv(test_set_file_path).dropna().drop_duplicates()

    OUTPUT_TYPE = OutputType(False, 'binary')
    UNIQUE_LABELS = [0, 1]
    OUTPUT_SPEC = OutputSpec(OUTPUT_TYPE, UNIQUE_LABELS)

    pretrained_model_generator, input_encoder = load_pretrained_model()

    with open("proteinbert_models/" + model, "rb") as f:  # 加载model_generator模型文件
        model_generator = pickle.load(f)

    # 预测结果
    result = predict_our(model_generator, input_encoder, OUTPUT_SPEC, test_set['seq'],
                         start_seq_len=512, start_batch_size=32, threshold=threshold)
    if os.path.exists(output_path + '.csv'):
        os.remove(output_path + '.csv')
    return result


# afp1_5相关参数 path和name是要预测的数据集的路径
path = 'data/predict_data/'  # 改这里
name = os.listdir(path)
output_path_1_5 = "predict_ready_1_5"
seq_1_5 = []
info_dict_1_5 = {}
model_1_5 = "model_afp1_5.pkl"
threshold_1_5 = 0.7


start_time = time.time()
# 遍历每个文件 构建afp1_5的测试集以及相应的字典
for each in name:
    species = str(each)
    species = species[:species.find('.faa')]  # 取出文件名
    position = 0
    for record in SeqIO.parse(path + each, "fasta"):
        sequence = str(record.seq)
        if '*' in sequence:
            sequence = sequence[:sequence.find('*')]  # 有些蛋白质最后以*结尾 去掉*即可
        seq_1_5.append(sequence)
        info_dict_1_5[sequence] = [species, position, None]
        position += 1

predict_1_5 = prediction(seq_1_5, output_path_1_5, model_1_5, threshold_1_5)

# 遍历每一个蛋白质 看看预测结果是真还是假 预测结果为真则输出相关信息 预测结果为假则删除该蛋白质 意味着1_5字典中只存了预测为真的蛋白质
for protein in predict_1_5:
    if predict_1_5[protein]:
        print('Sequence: ', protein, 'Predict: ', predict_1_5[protein], 'Species: ', info_dict_1_5[protein][0])
    else:
        info_dict_1_5.pop(protein, None)

# 再循环一遍文件 用一个dict表示1_5管辖的上下10个11 以及新的seq_11 然后用seq_11来输入afp_11模型预测结果
# 根据1_5和11都对来判断

# 设置afp11相关参数
seq_11_or_afp13 = []
dict_1_5to11_and_13 = {}  # {afp1_5: [afp11, afp11, afp11] }
surrounding = 10  # 上下10个
output_path_11 = "predict_ready_11"
model_11 = "model_afp11.pkl"
threshold_11 = 0.74

# 设置afp13相关参数
output_path_13 = "predict_ready_13"
model_13 = "model_afp13_threshold_0.7.pkl"
threshold_13 = 0.7

for sequence_afp1_5 in info_dict_1_5.keys():
    species = info_dict_1_5[sequence_afp1_5][0]
    sequences_all = [record for record in SeqIO.parse(path + species + '.faa', "fasta")]
    position = info_dict_1_5[sequence_afp1_5][1]
    info_dict_1_5[sequence_afp1_5][1] = sequences_all[position].description
    start = max(0, position - surrounding)  # 起点
    end = min(len(sequences_all), position + surrounding + 1)  # 终点
    for j in range(start, end):  # 此循环中遍历到的蛋白质就是上下各10个蛋白质
        if position != j:
            sequence_surrounding = str(sequences_all[j].seq)
            if '*' in sequence_surrounding:
                sequence_surrounding = sequence_surrounding[:sequence_surrounding.find('*')]
            seq_11_or_afp13.append(sequence_surrounding)
            if sequence_afp1_5 not in dict_1_5to11_and_13:
                dict_1_5to11_and_13[sequence_afp1_5] = [sequence_surrounding]
            else:
                dict_1_5to11_and_13[sequence_afp1_5].append(sequence_surrounding)

    # position = 0
    # for record in SeqIO.parse(path + species + '.faa', "fasta"):
    #     sequence = str(record.seq)
    #     if '*' in sequence:
    #         sequence = sequence[:sequence.find('*')]
    #     if sequence_afp1_5 == sequence:
    #         info_dict_1_5[sequence_afp1_5][1] = record.description
    #         start = max(0, position - surrounding)  # 起点
    #         end = min(len(sequences_all), position + surrounding + 1)  # 终点
    #         for j in range(start, end):  # 此循环中遍历到的蛋白质就是上下各10个蛋白质
    #             if position != j:
    #                 sequence_surrounding = str(sequences_all[j].seq)
    #                 if '*' in sequence_surrounding:
    #                     sequence_surrounding = sequence_surrounding[:sequence_surrounding.find('*')]
    #                 seq_11_or_afp13.append(sequence_surrounding)
    #                 if sequence_afp1_5 not in dict_1_5to11_and_13:
    #                     dict_1_5to11_and_13[sequence_afp1_5] = [sequence_surrounding]
    #                 else:
    #                     dict_1_5to11_and_13[sequence_afp1_5].append(sequence_surrounding)
    #     else:
    #         position += 1

predict_1_5_and_11 = prediction(seq_11_or_afp13, output_path_11, model_11, threshold_11)
predict_1_5_and_13 = prediction(seq_11_or_afp13, output_path_13, model_13, threshold_13)

# 修改版
positive_afp_11 = []
# 遍历1_5字典中的所有蛋白质 如果此蛋白质的上下各10个蛋白质组合中有一个蛋白质为真则输出最后结果 并且输出这个1_5蛋白质的相关信息
for afp_1_5 in info_dict_1_5.keys():
    for protein_afp11 in dict_1_5to11_and_13[afp_1_5]:
        if predict_1_5_and_11[protein_afp11]:
            positive_afp_11.append(protein_afp11)
    for protein_afp13 in dict_1_5to11_and_13[afp_1_5]:
        if protein_afp13 not in positive_afp_11 and predict_1_5_and_13[protein_afp13]:
            print('This is a True eCIS!')
            print('Protein:', afp_1_5, ' Species: ', info_dict_1_5[afp_1_5][0], 'Description: ',
                  info_dict_1_5[afp_1_5][1])
            break

end_time = time.time()
print('Time: ', end_time - start_time)

# 原版
# for afp_1_5 in info_dict_1_5.keys():
#     for protein_afp11 in dict_1_5to11_and_13[afp_1_5]:
#         if predict_1_5_and_11[protein_afp11]:
#             for protein_afp13 in dict_1_5to11_and_13[afp_1_5]:
#                 if protein_afp13 != protein_afp11 and predict_1_5_and_13[protein_afp13]:
#                     print('This is a True eCIS!')
#                     print('Protein:', afp_1_5, ' Species: ', info_dict_1_5[afp_1_5][0], 'Description: ',
#                           info_dict_1_5[afp_1_5][1])
#                     break
#             break
