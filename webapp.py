import re
from collections import Counter
import pandas as pd
import pickle
import numpy as np
import pickle
import webbrowser
from threading import Timer
import joblib
import base64
import pep as p
from sklearn import preprocessing
import sys
import time
from tempfile import NamedTemporaryFile
import os
import streamlit as st

img1 ="2-1.png"


st.markdown(
    """
    <style>
    .logo-img {
    display: inline-block;
    margin-left: auto;
    margin-right: auto;
    height: 120px; 
  }



  #images{
    text-align:center;
  }
  </style>
    """,
    unsafe_allow_html=True
    )

st.markdown(
    f"""
    <div class="container">
      <div class="column">
        <img class="logo-img" src="data:image/png;base64,{base64.b64encode(open(img1, "rb").read()).decode()}">
        </div>
    </div>
    """,
    unsafe_allow_html=True
)



# Set page title
st.title('Anti-Cancer Peptide Prediction Tool')


# Preprocess function
alphabets = ['A', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'K', 'L', 'M', 'N', 'P', 'Q', 'R', 'S', 'T', 'V', 'W', 'Y']
punct = '!?,.@#'
maxlen = 280

st.subheader('Single Peptide Predictor')
str_features_input = st.text_input('')
if str_features_input!= '':
    str_features=[]
    str_features += [str_features_input]
    print(type(str_features))
    print(str_features)
    f= open("./checkfasta.txt","w+")
    count=1
    for i in range(len(str_features)):
        f.write(">" + str(count) + "\n" + str_features[i])
        count = count + 1
    f.close()

    fastas=p.readFasta('./checkfasta.txt')

    GDPC=p.GDPC(fastas)
    GAAC=p.GAAC(fastas)
    AAC=p.AAC(fastas)
    GTPC=p.GTPC(fastas)
    DPC=p.DPC(fastas)
    TPC=p.TPC(fastas)
    CTDC=p.CTDC(fastas)
    CTDD=p.CTDD(fastas) 
    CTDT=p.CTDT(fastas)

    GDPC=pd.DataFrame(data=GDPC)
    GDPC.columns = GDPC.iloc[0]
    GDPC=GDPC.drop(GDPC.index[0])
    GDPC.drop(['#'],inplace=True,axis=1)
    # GDPC.to_csv('./GDPC_check.csv')

    GAAC=pd.DataFrame(data=GAAC)
    GAAC.columns = GAAC.iloc[0]
    GAAC=GAAC.drop(GAAC.index[0])
    GAAC.drop(['#'],inplace=True,axis=1)

    AAC=pd.DataFrame(data=AAC)
    AAC.columns = AAC.iloc[0]
    AAC=AAC.drop(AAC.index[0])
    AAC.drop(['#'],inplace=True,axis=1)

    GTPC=pd.DataFrame(data=GTPC)
    GTPC.columns = GTPC.iloc[0]
    GTPC=GTPC.drop(GTPC.index[0])
    GTPC.drop(['#'],inplace=True,axis=1)

    DPC=pd.DataFrame(data=DPC)
    DPC.columns = DPC.iloc[0]
    DPC=DPC.drop(DPC.index[0])
    DPC.drop(['#'],inplace=True,axis=1)

    TPC=pd.DataFrame(data=TPC)
    TPC.columns = TPC.iloc[0]
    TPC=TPC.drop(TPC.index[0])
    TPC.drop(['#'],inplace=True,axis=1)

    CTDC=pd.DataFrame(data=CTDC)
    CTDC.columns = CTDC.iloc[0]
    CTDC=CTDC.drop(CTDC.index[0])
    CTDC.drop(['#'],inplace=True,axis=1)

    CTDD=pd.DataFrame(data=CTDD)
    CTDD.columns = CTDD.iloc[0]
    CTDD=CTDD.drop(CTDD.index[0])
    CTDD.drop(['#'],inplace=True,axis=1)

    CTDT=pd.DataFrame(data=CTDT)
    CTDT.columns = CTDT.iloc[0]
    CTDT=CTDT.drop(CTDT.index[0])
    CTDT.drop(['#'],inplace=True,axis=1)

    df1=pd.concat([TPC,DPC,GAAC,GDPC,AAC,GTPC,CTDC,CTDD,CTDT],axis=1)
    df1.to_csv('df1_check.csv',index=False)

    with open('peptide_rf_new.pkl', 'rb') as model_file:
        model = pickle.load(model_file)

    prediction = model.predict(df1.loc[0:1:2])
    df=pd.DataFrame(data=zip(fastas,prediction),columns=['peptide','Predictions'])
    df.to_csv('peptide_db.csv',mode='a',header=False,index=False)
    st.success(prediction)
    

else:
    st.write('Please enter peptide sequence')

###File Upload###
temp = st.file_uploader("",type=['txt'])
#temp = temp.decode()

buffer = temp
temp_file = NamedTemporaryFile(delete=False)
if buffer:
    temp_file.write(buffer.getvalue())
    st.write((temp.name))


if buffer is None:
    st.text("Please upload file (.txt)")

else:
    f= open("./checkfastafile.txt","w+")
    for i in temp:
        clean_line = i.strip().decode( "utf-8" )
        f.write(str(clean_line)+'\n')
    print(f)
    f.close()

    fastas=p.readFasta('./checkfastafile.txt')
    GDPC=p.GDPC(fastas)
    GAAC=p.GAAC(fastas)
    AAC=p.AAC(fastas)
    GTPC=p.GTPC(fastas)
    DPC=p.DPC(fastas)
    TPC=p.TPC(fastas)
    CTDC=p.CTDC(fastas)
    CTDD=p.CTDD(fastas) 
    CTDT=p.CTDT(fastas)

    GDPC=pd.DataFrame(data=GDPC)
    GDPC.columns = GDPC.iloc[0]
    GDPC=GDPC.drop(GDPC.index[0])
    GDPC.drop(['#'],inplace=True,axis=1)

    GAAC=pd.DataFrame(data=GAAC)
    GAAC.columns = GAAC.iloc[0]
    GAAC=GAAC.drop(GAAC.index[0])
    GAAC.drop(['#'],inplace=True,axis=1)

    AAC=pd.DataFrame(data=AAC)
    AAC.columns = AAC.iloc[0]
    AAC=AAC.drop(AAC.index[0])
    AAC.drop(['#'],inplace=True,axis=1)

    GTPC=pd.DataFrame(data=GTPC)
    GTPC.columns = GTPC.iloc[0]
    GTPC=GTPC.drop(GTPC.index[0])
    GTPC.drop(['#'],inplace=True,axis=1)

    DPC=pd.DataFrame(data=DPC)
    DPC.columns = DPC.iloc[0]
    DPC=DPC.drop(DPC.index[0])
    DPC.drop(['#'],inplace=True,axis=1)

    TPC=pd.DataFrame(data=TPC)
    TPC.columns = TPC.iloc[0]
    TPC=TPC.drop(TPC.index[0])
    TPC.drop(['#'],inplace=True,axis=1)

    CTDC=pd.DataFrame(data=CTDC)
    CTDC.columns = CTDC.iloc[0]
    CTDC=CTDC.drop(CTDC.index[0])
    CTDC.drop(['#'],inplace=True,axis=1)

    CTDD=pd.DataFrame(data=CTDD)
    CTDD.columns = CTDD.iloc[0]
    CTDD=CTDD.drop(CTDD.index[0])
    CTDD.drop(['#'],inplace=True,axis=1)

    CTDT=pd.DataFrame(data=CTDT)
    CTDT.columns = CTDT.iloc[0]
    CTDT=CTDT.drop(CTDT.index[0])
    CTDT.drop(['#'],inplace=True,axis=1)

    df1=pd.concat([TPC,DPC,GAAC,GDPC,AAC,GTPC,CTDC,CTDD,CTDT],axis=1)
    # df1.to_csv('C:/Users/ihamm/peptide(2)/peptide/file_upload/features_file.csv',index=False)

    with open('peptide_rf_new.pkl', 'rb') as model_file:
        model = pickle.load(model_file)

    prediction = model.predict(df1)
    df=pd.DataFrame(data=zip(fastas,prediction),columns=['Peptide','Predictions'])
    # if os.path.isfile('C:/Users/ihamm/peptide(2)/peptide/file_upload/peptide_db.csv'):
    #   os.remove('C:/Users/ihamm/peptide(2)/peptide/file_upload/peptide_db.csv')
    st.success("Prediction file generated")
    csv = df.to_csv(index=False)
    st.download_button('Download!', csv, file_name='.%d'%time.time()+'.csv')
    del csv
    # b64 = base64.b64encode(csv.encode()).decode()  # some strings <-> bytes conversions necessary here
    # href = f'<a href="data:file/csv;base64,{b64}">Download CSV File</a> (right-click and save as &lt;some_name&gt;.csv)'
    # st.markdown(href, unsafe_allow_html=True)
