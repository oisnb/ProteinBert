import pickle
import pandas as pd
import os
from proteinbert import OutputType, OutputSpec, FinetuningModelGenerator, load_pretrained_model, finetune, evaluate_by_len
from proteinbert.finetuning import predict_our
from Bio import SeqIO

# 完整的预测函数
def prediction(seqences_data, output_path, model ,threshold):
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

    with open("proteinbert_models/" + model, "rb") as f:  #加载model_generator模型文件
        model_generator = pickle.load(f)

    # 预测结果
    result = predict_our(model_generator, input_encoder, OUTPUT_SPEC, test_set['seq'], \
                               start_seq_len=512, start_batch_size=32, threshold=threshold)
    if os.path.exists(output_path + '.csv'):
        os.remove(output_path + '.csv')
    return result

# afp1_5相关参数 path和name是要预测的数据集的路径
path = 'data/predict_data/' #改这里
name = os.listdir(path)
output_path_1_5 = "predict_ready_1_5"
seq_1_5 = []
info_dict_1_5 = {}
model_1_5 = "model_afp1_5.pkl"
threshold_1_5 = 0.7

# 遍历每个文件 构建afp1_5的测试集以及相应的字典
for each in name:
    species = str(each)
    species = species[:species.find('.faa')] # 取出文件名
    for record in SeqIO.parse(path + each, "fasta"):
        sequence = str(record.seq)
        if '*' in sequence:
            sequence = sequence[:sequence.find('*')] # 有些蛋白质最后以*结尾 去掉*即可
        seq_1_5.append(sequence)
        info_dict_1_5[sequence] = [species, record.description]


predict_1_5 = prediction(seq_1_5, output_path_1_5, model_1_5, threshold_1_5)

#遍历每一个蛋白质 看看预测结果是真还是假 预测结果为真则输出相关信息 预测结果为假则删除该蛋白质 意味着1_5字典中只存了预测为真的蛋白质
for protein in predict_1_5:
    if predict_1_5[protein] == True:
        print('Sequence: ', protein, 'Predict: ', predict_1_5[protein], 'Species: ', info_dict_1_5[protein][0], 'Description: ', info_dict_1_5[protein][1])
    else:
        info_dict_1_5.pop(protein, None)


# 再循环一遍文件 用一个dict表示1_5管辖的上下10个11 以及新的seq_11 然后用seq_11来输入afp_11模型预测结果
# 根据1_5和11都对来判断

# 设置afp11相关参数
seq_11 = []
# info_dict_11 = {}
dict_1_5to11 = {} # {afp1_5: [afp11, afp11, afp11] }
surrounding = 10 # 上下10个
output_path_11 = "predict_ready_11"
model_11 = "model_afp11.pkl"
threshold_11 = 0.7

# 遍历取上下各10个蛋白质
for each in name:
    sequences_all = [record for record in SeqIO.parse(path + each, "fasta")]
    for record in SeqIO.parse(path + each, "fasta"):
        sequence_afp1_5 = str(record.seq)
        if '*' in sequence_afp1_5:
            sequence_afp1_5 = sequence_afp1_5[:sequence_afp1_5.find('*')]
        if sequence_afp1_5 in info_dict_1_5: # 如果当前查找的蛋白质在1_5字典中则查找其上下各10个蛋白质
            for i, seq in enumerate(sequences_all):
                if seq.id == record.id: # 先找到这个蛋白质的位置
                    position = i
                    break
            start = max(0, position - surrounding) # 起点
            end = min(len(sequences_all), position + surrounding + 1) # 终点
            for j in range(start, end): # 此循环中遍历到的蛋白质就是上下各10个蛋白质
                if position != j:
                    sequence_afp11 = str(sequences_all[j].seq)
                    if '*' in sequence_afp11:
                        sequence_afp11 = sequence_afp11[:sequence_afp11.find('*')]
                    seq_11.append(sequence_afp11)
                    if sequence_afp1_5 not in dict_1_5to11:
                        dict_1_5to11[sequence_afp1_5] = [sequence_afp11]
                    else:
                        dict_1_5to11[sequence_afp1_5].append(sequence_afp11)

predict_1_5_and_11 = prediction(seq_11, output_path_11, model_11, threshold_11)

# 遍历1_5字典中的所有蛋白质 如果此蛋白质的上下各10个蛋白质组合中有一个蛋白质为真则输出最后结果 并且输出这个1_5蛋白质的相关信息
for afp_1_5 in info_dict_1_5.keys():
    for protein in dict_1_5to11[afp_1_5]:
        if predict_1_5_and_11[protein] == True:
            print('This is a True eCIS!')
            print('Species: ', info_dict_1_5[afp_1_5][0], 'Description: ', info_dict_1_5[afp_1_5][1])
            break
