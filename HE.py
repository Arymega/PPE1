import pandas as pd
import numpy as np

RAWSS=pd.DataFrame(pd.read_excel('RAW/RAWHE.xlsx'))

shape=RAWSS.shape
noi=shape[1]-1
nop=shape[0]
parraw=np.zeros((noi,nop))
ind=RAWSS.columns
for i in range(0,noi):
    parraw[i,:]=RAWSS[ind[i+1]]

def MASTER(par):
    LMTD=((par[0]-par[3])-(par[2]-par[1]))/np.log((par[0]-par[3])/(par[2]-par[1]))
    FT=par[12]/LMTD
    if FT>1:
        FT=1
    MTDreal=LMTD*FT    
    UD=par[13]/par[14]/MTDreal
    UC=1/(1/UD-par[16])
    return [par[13],par[14],FT,UD,UC,par[16]]

parres=np.zeros((noi,6))
parres[0,:]=MASTER(parraw[0,:])
RES=pd.DataFrame({ind[1]:parres[0,:]})
for i in range(0,noi):
    parres[i,:]=MASTER(parraw[i,:])
    RES[ind[i+1]]=pd.Series(parres[i,:])
RES=RES.rename(index={0:"Q (kcal/h)",1:"A (m2)",2:"FT (MTD/LMTD)",\
                      3:"UD (kcal/m2/h/K)",4:"UC (kcal/m2/h/K)",\
                      5:"Rd (m2 h K/kcal)"})
RES.to_excel('RES/HE.xlsx')