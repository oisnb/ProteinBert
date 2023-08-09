import pickle
import pandas as pd
import os
from proteinbert import OutputType, OutputSpec, FinetuningModelGenerator, load_pretrained_model, finetune, evaluate_by_len
from proteinbert.finetuning import predict_our
from Bio import SeqIO

with open('data/afp1_5.pkl', 'rb') as file:
    info_dict_1_5 = pickle.load(file)
#再循环一遍文件 用一个dict表示1_5管辖的上下10个11 以及新的seq_11 然后用seq_11来输入afp_11模型预测结果
#根据1_5和11都对来判断

path = 'data/predict_data/' #改这里
name = os.listdir(path)
output_path = "predict_ready"
seq_11 = []
# info_dict_11 = {}
dict_1_5to11 = {} #{afp1_5: [afp11, afp11, afp11] }
surrounding = 10 #上下10个基因

for each in name:
    sequences_all = [record for record in SeqIO.parse(path + each, "fasta")]
    for record in SeqIO.parse(path + each, "fasta"):
        sequence_afp1_5 = str(record.seq)
        if '*' in sequence_afp1_5:
            sequence_afp1_5 = sequence_afp1_5[:sequence_afp1_5.find('*')]
        if sequence_afp1_5 in info_dict_1_5:
            for i, seq in enumerate(sequences_all):
                if seq.id == record.id:
                    position = i
                    break
            start = max(0, position - surrounding)
            end = min(len(sequences_all), position + surrounding + 1)
            for j in range(start, end):
                if position != j:
                    sequence_afp11 = str(sequences_all[j].seq)
                    if '*' in sequence_afp11:
                        sequence_afp11 = sequence_afp11[:sequence_afp11.find('*')]
                    seq_11.append(sequence_afp11)
                    if sequence_afp1_5 not in dict_1_5to11:
                        dict_1_5to11[sequence_afp1_5] = [sequence_afp11]
                    else:
                        dict_1_5to11[sequence_afp1_5].append(sequence_afp11)


df = pd.DataFrame(data=seq_11, columns=['seq'])
df.to_csv(output_path + '.csv', index=False)

test_set_file_path = os.path.join(output_path + '.csv')
test_set = pd.read_csv(test_set_file_path).dropna().drop_duplicates()

OUTPUT_TYPE = OutputType(False, 'binary')
UNIQUE_LABELS = [0, 1]
OUTPUT_SPEC = OutputSpec(OUTPUT_TYPE, UNIQUE_LABELS)

pretrained_model_generator, input_encoder = load_pretrained_model()

with open("proteinbert_models/model_afp11.pkl", "rb") as f: #记得改 改成按afp11的generator
    model_generator = pickle.load(f)

y_pred_afp11 = predict_our(model_generator, input_encoder, OUTPUT_SPEC, test_set['seq'], \
        start_seq_len = 512, start_batch_size = 32)

for afp_1_5 in info_dict_1_5.keys():
    for protein in dict_1_5to11[afp_1_5]:
        if y_pred_afp11[protein] == True:
            print('This is a True eCIS!')
            print('Species: ', info_dict_1_5[afp_1_5][0], 'Description: ', info_dict_1_5[afp_1_5][1])
            break