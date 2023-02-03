from os.path import join, isfile
from os import listdir
from threading import Thread
from queue import Queue
from Bio import SearchIO
from Bio.Blast import NCBIXML
from Bio.Blast import NCBIWWW
from Bio.Seq import Seq
from Bio import SeqIO
from Bio import motifs
from Bio.SeqUtils import GC
import threading
import pandas as pd 
import numpy as np
import itertools
import pickle
import os


import warnings
warnings.filterwarnings("ignore")


li=['A','C','G','T']

li_2 = []
for i in itertools.product(['A','C','G','T'], repeat=2):
    li_2.append(''.join(map(str, i)))


li_3 = []
for i in itertools.product(['A','C','G','T'], repeat=3):
    li_3.append(''.join(map(str, i)))


li_4 = []
for i in itertools.product(['A','C','G','T'], repeat=4):
    li_4.append(''.join(map(str, i)))


def lncrna_feat(file_name):
    print('LncRNA Feature Generation')
    lncrna_df=pd.DataFrame(columns=['lncrna','RNA_sequence'])
    i=0
    for seq_record in SeqIO.parse(file_name, "fasta"):
        lncrna_df.loc[i,'lncrna']=seq_record.id
        lncrna_df.loc[i,'RNA_sequence']=str(seq_record.seq)
        i=i+1
    lncrna_df['RNA_sequence']=lncrna_df['RNA_sequence'].str.upper()
    lncrna_df['RNA_sequence']=lncrna_df['RNA_sequence'].str.replace('U','T')

    for i in range(0,lncrna_df.shape[0]):
        lncrna_df.loc[i,'lncrna_seq_length']=len(lncrna_df.loc[i,'RNA_sequence'])
        lncrna_df.loc[i,'lncrna_GC_perc']=GC(lncrna_df.loc[i,'RNA_sequence']) 
        
        
    for i in range(0,lncrna_df.shape[0]):
        for j in li:
            lncrna_df.loc[i,j+'_lncrna_count']=lncrna_df.loc[i,'RNA_sequence'].count(j)


        for j in li_2:
            lncrna_df.loc[i,j+'_lncrna_count']=lncrna_df.loc[i,'RNA_sequence'].count(j)


        for j in li_3:
            lncrna_df.loc[i,j+'_lncrna_count']=lncrna_df.loc[i,'RNA_sequence'].count(j)

        for j in li_4:
            lncrna_df.loc[i,j+'_lncrna_count']=lncrna_df.loc[i,'RNA_sequence'].count(j)
            
    return lncrna_df
    

def mrna_feat(file_name):
    print('mRNA Feature Generation')
    
    mrna_df=pd.DataFrame(columns=['mrna','RNA_sequence'])
    i=0
    for seq_record in SeqIO.parse(file_name, "fasta"):
        mrna_df.loc[i,'mrna']=seq_record.id
        mrna_df.loc[i,'RNA_sequence']=str(seq_record.seq)
        i=i+1
    mrna_df['RNA_sequence']=mrna_df['RNA_sequence'].str.upper()
    mrna_df['RNA_sequence']=mrna_df['RNA_sequence'].str.replace('U','T')
    
    for i in range(0,mrna_df.shape[0]):
        mrna_df.loc[i,'mrna_seq_length']=len(mrna_df.loc[i,'RNA_sequence'])
        mrna_df.loc[i,'mrna_GC_perc']=GC(mrna_df.loc[i,'RNA_sequence']) 
        
        
    for i in range(0,mrna_df.shape[0]):
        for j in li:
            mrna_df.loc[i,j+'_mrna_count']=mrna_df.loc[i,'RNA_sequence'].count(j)


        for j in li_2:
            mrna_df.loc[i,j+'_mrna_count']=mrna_df.loc[i,'RNA_sequence'].count(j)


        for j in li_3:
            mrna_df.loc[i,j+'_mrna_count']=mrna_df.loc[i,'RNA_sequence'].count(j)

        for j in li_4:
            mrna_df.loc[i,j+'_mrna_count']=mrna_df.loc[i,'RNA_sequence'].count(j)
    
    return mrna_df


def merge_data(lncrna_df,mrna_df):
    print('Merging Data')
    
    lncrna_df.drop(['RNA_sequence'],inplace=True,axis=1)
    mrna_df.drop(['RNA_sequence'],inplace=True,axis=1)
    lnc_mrna_df=lncrna_df.assign(key=1).merge(mrna_df.assign(key=1)).drop('key', 1)
    lnc_m_ind=lnc_mrna_df[['lncrna','mrna']]
    lnc_mrna_df.drop(['lncrna','mrna'],inplace=True,axis=1)
    return lnc_m_ind,lnc_mrna_df


def make_pred(lnc_m_ind,lnc_mrna_df):
    print('Predicting')
    
    
    filename_1 = 'Models/LNCRNA_MRNA_NORMALIZER_MOUSE.sav'
    filename_2 = 'Models/LNCRNA_MRNA_PCA_MOUSE.sav'
    normaling=pickle.load(open(filename_1, 'rb'))
    pcaing=pickle.load(open(filename_2, 'rb'))
    
    X_Test=lnc_mrna_df.values
    X_Test_Norm=normaling.transform(X_Test)
    X_Test_PCA=pcaing.transform(X_Test_Norm)[:,:10]  
    
    model_TR = pickle.load(open('Models/Decision_Tree_MOUSE.sav', 'rb'))
    model_KNN = pickle.load(open('Models/KNN_MOUSE.sav', 'rb'))
    model_RF = pickle.load(open('Models/RF_MOUSE.sav', 'rb'))
    model_LGBM = pickle.load(open('Models/LGBM_MOUSE.sav', 'rb'))
    
    lnc_m_ind['Decision_Tree(%)']=np.round(model_TR.predict_proba(X_Test_PCA)[:,1]*100,2)
    lnc_m_ind['K-Nearest-Neighbours(%)']=np.round(model_KNN.predict_proba(X_Test_PCA)[:,1]*100,2)
    lnc_m_ind['Random_Forest(%)']=np.round(model_RF.predict_proba(X_Test_PCA)[:,1]*100,2)
    lnc_m_ind['LightGBM(%)']=np.round(model_LGBM.predict_proba(X_Test_PCA)[:,1]*100,2)
    
    lnc_m_ind['Decision_Tree']=model_TR.predict(X_Test_PCA)
    lnc_m_ind['K-Nearest-Neighbours']=model_KNN.predict(X_Test_PCA)
    lnc_m_ind['Random_Forest']=model_RF.predict(X_Test_PCA)
    lnc_m_ind['LightGBM']=model_LGBM.predict(X_Test_PCA)

    lnc_m_ind['Final_Target_Score']=np.round((lnc_m_ind['Decision_Tree']+lnc_m_ind['K-Nearest-Neighbours']+lnc_m_ind['Random_Forest']+lnc_m_ind['LightGBM'])*100/4,2)

    lnc_m_ind.drop(['Decision_Tree','K-Nearest-Neighbours','Random_Forest','LightGBM'],inplace=True,axis=1)


    lnc_m_ind.to_csv('Output_Data/Predicted_Val.csv',index=False)

if __name__ == "__main__":
    lnc_file='Input_Data/query_lncrna.fa'
    m_file='Input_Data/target_mrna.fa'
    lncrna_df=lncrna_feat(lnc_file)
    mrna_df=mrna_feat(m_file)
    lnc_m_ind,lnc_mrna_df=merge_data(lncrna_df,mrna_df)
    make_pred(lnc_m_ind,lnc_mrna_df)

