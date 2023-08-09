import pickle
import pandas as pd
import os
from proteinbert import OutputType, OutputSpec, FinetuningModelGenerator, load_pretrained_model, finetune, evaluate_by_len
from proteinbert.finetuning import predict_our
from Bio import SeqIO


path_afp1_5 = 'data/predict_data/' #改这里
name = os.listdir(path_afp1_5)
output_path = "predict_ready"
seq = []
info_dict = {}

for each in name:
    species = str(each)
    species = species[:species.find('.faa')]
    for record in SeqIO.parse(path_afp1_5 + each, "fasta"):
        sequence = str(record.seq)
        if '*' in sequence:
            sequence = sequence[:sequence.find('*')]
        seq.append(sequence)
        info_dict[sequence] = [species, record.description]



# fa_seq = SeqIO.parse(path_afp1_5, "fasta") #迭代器读取
# for record in fa_seq:
#     print(1)
#     for feature in record.features:
#         if feature.type == "CDS" and "translation" in feature.qualifiers:
#             protein_sequence = feature.qualifiers["translation"][0]
#             seq.append(protein_sequence)
#             info_dict[protein_sequence] = [record.description, record.name, feature.location]

df = pd.DataFrame(data=seq, columns=['seq'])
df.to_csv(output_path + '.csv', index=False)



test_set_file_path = os.path.join(output_path + '.csv')
test_set = pd.read_csv(test_set_file_path).dropna().drop_duplicates()

OUTPUT_TYPE = OutputType(False, 'binary')
UNIQUE_LABELS = [0, 1]
OUTPUT_SPEC = OutputSpec(OUTPUT_TYPE, UNIQUE_LABELS)

pretrained_model_generator, input_encoder = load_pretrained_model()

with open("proteinbert_models/model_afp1_5.pkl", "rb") as f:
    model_generator = pickle.load(f)

y_pred = predict_our(model_generator, input_encoder, OUTPUT_SPEC, test_set['seq'], \
        start_seq_len = 512, start_batch_size = 32)

for protein in y_pred:
    if y_pred[protein] == True:
        print('Sequence: ', protein, 'Predict: ', y_pred[protein], 'Species: ', info_dict[protein][0], 'Description: ', info_dict[protein][1])
    else:
        info_dict.pop(protein, None)

with open('data/afp1_5.pkl', 'wb') as file:
    pickle.dump(info_dict, file)

if os.path.exists(output_path + '.csv'):
    os.remove(output_path + '.csv')

