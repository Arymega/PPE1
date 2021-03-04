import pandas as pd
import numpy as np
from scipy.integrate import quad
import HE

RAWCPV=pd.DataFrame(pd.read_excel('RAW/CP.xlsx',sheet_name='Gas'))
RAWCPL=pd.DataFrame(pd.read_excel('RAW/CP.xlsx',sheet_name='Liquid'))
RAWANT=pd.DataFrame(pd.read_excel('RAW/CP.xlsx',sheet_name='Antoine'))
RAWDHV=pd.DataFrame(pd.read_excel('RAW/CP.xlsx',sheet_name='Vaporization Heat'))
RAWLD=pd.DataFrame(pd.read_excel('RAW/CP.xlsx',sheet_name='Liquid Density'))
RAWOP=pd.DataFrame(pd.read_excel('RAW/HE101JC1.xlsx'))
RAWSS=pd.DataFrame(pd.read_excel('RES/HE.xlsx'))

SS=RAWSS['101-JC1']
q=SS[0]
a=SS[1]
ft=SS[2]
ud=SS[3]
uc=SS[4]
rd=SS[5]
DESIGN=[q,ud,rd,1,1]

def LD(T,PAR):
    T=T+273.15
    return PAR[0]+PAR[1]*T+PAR[2]*T**2+PAR[3]*T**3
def CPL(T,PAR):
    T=T+273.15
    return PAR[0]+PAR[1]*T+PAR[2]*T**2+PAR[3]*T**3+PAR[4]*T**4
def HL(T,PAR):
    return quad(CPL,0,T,args=PAR)
def CPV(T,PAR):
    T=T+273.15
    return 8314*(PAR[0]+PAR[1]*T+PAR[2]*T**2+PAR[3]*T**-2)
def HV(T,PAR):
    return quad(CPV,0,T,args=PAR)

OP=RAWOP['101-JC1']
MTD=ft*((OP[0]-OP[3])-(OP[2]-OP[1]))/np.log((OP[0]-OP[3])/(OP[2]-OP[1]))

CPLWater=np.array(RAWCPL.loc[RAWCPL['J/kmol/K']=='Water']).flatten()
CPLWater=CPLWater[1:6]
CPVAir=np.array(RAWCPV.loc[RAWCPV['Heat Capacity/R']=='Air']).flatten()
CPVAir=CPVAir[1:5]

def UDRDL(Tin,Tout,F): #T in C, F in kmol/h
    Hin=F*HL(Tin,CPLWater)[0]
    Hout=F*HL(Tout,CPLWater)[0]
    Q=(Hout-Hin)/1e3/4.1868 #kcal/h
    UD=abs(Q/a/MTD)
    RD=1/UD-1/uc
    EFF=UD/ud
    HEATLOAD=abs(Q/q)
    return [Q,UD,RD,EFF,HEATLOAD]

def UDRDV(Tin,Tout,F): #T in C, F in kmol/h
    Hin=F*HV(Tin,CPVAir)[0]
    Hout=F*HV(Tout,CPVAir)[0]
    Q=(Hout-Hin)/1e3/4.1868 #kcal/h
    UD=abs(Q/a/MTD)
    RD=1/UD-1/uc
    EFF=UD/ud
    HEATLOAD=abs(Q/q)
    return [Q,UD,RD,EFF,HEATLOAD]

#SHELL
TSIN=OP[0]
TSOUT=OP[2]
FS=OP[4]*101325/8.314/273.15/1e3
SHELL=UDRDV(TSIN,TSOUT,FS)

#TUBE
TTIN=OP[1]
TTOUT=OP[3]
FT=OP[5]/18
TUBE=UDRDL(TTIN,TTOUT,FT)

RES=pd.DataFrame({'DESIGN':DESIGN,'SHELL':SHELL,'TUBE':TUBE})
RES=RES.rename(index={0:"Q (kcal/h)",1:"Ud (kcal/m2/h/K)",\
                      2:"Rd (m2 h K/kcal)",3:"Efficiency (Ud/Ud design)",\
                      4:"Heat load (Q/Q design)"})
print(RES)
RES.to_excel('RES/101-JC1.xlsx')