import os

from Bio.Cluster import kcluster
import numpy as np
import pandas as pd
import random

from matplotlib import pyplot as plt
from sklearn import metrics
from sklearn.cluster import KMeans
from sklearn.decomposition import PCA
from sklearn.feature_extraction.text import TfidfVectorizer
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.SeqFeature import SeqFeature, FeatureLocation
from Bio import SeqIO
import random

source = []
locus_tag = []
DNA = []
Protein = []
path = 'D://mypython//data//all_DNA_genomics//'

sum = 0  # 总个数
sample_numbers = 10  # 期望个数
counter = 0

temp_gene = None
temp_sequence = None
for file in os.listdir(path):
    position = file.find("_cds")
    replicon = file[0:position]
    for record in SeqIO.parse(path + file, "fasta"):
        number = random.random()
        if number > 0.99 and sum < sample_numbers:
            source.append(str(replicon))
            description = record.description
            start = description.find('tag=')
            end = description.find('] [')
            if start != -1:
                component = description[start + 4:end]
            else:
                component = ''
            locus_tag.append(str(component))
            seq = record.seq
            protein_seq = seq.translate()
            DNA.append(str(seq))
            Protein.append(str(protein_seq))
            sum += 1
        elif sum >= sample_numbers:
            counter = 2
            break
    if counter == 2:
        break

#  聚类
vectorizer = TfidfVectorizer()
tfidf_matrix = vectorizer.fit_transform(Protein)
kmeans = KMeans(n_clusters=2, random_state=0, n_init='auto')
kmeans.fit(tfidf_matrix)
labels = kmeans.labels_
result = np.vstack((source, locus_tag, DNA, Protein, labels)).T
data = pd.DataFrame(data=result, columns=['Replicon', 'Component', 'DNASet', 'ProteinSet', 'clusterId'])
records = pd.DataFrame(data=None, columns=['Replicon', 'Component', 'DNASet', 'ProteinSet', 'clusterId'])
#  按照聚类结果分组
groups = data.groupby('clusterId')
for name, group in groups:
    random_records = group.sample(1)
    records = pd.concat([records, random_records])

records.to_csv('Negative_clustered_10000.csv', index=False)

# 评价效果
print(metrics.silhouette_score(tfidf_matrix, labels, metric='cosine'))
