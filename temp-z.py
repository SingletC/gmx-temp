import groio
import pandas as pd
import matplotlib.pyplot as plt
import os

path = "./100nm-0/"
NA=6.02e23
rNA=1/NA*0.001   # g to kg
mass_table={'OW':16*rNA,
            'HW':1*rNA,
            'H0':1*rNA,
            'C0':12*rNA,
            'N0':14*rNA,
            'MW':0*rNA,
            'CL':35.5*rNA,
            'NA':23*rNA}
deg_freedom={'OW':2,
            'HW':2,
            'H0':3,
            'C0':3,
            'N0':3,
            'MW':3,
            'CL':3,
            'NA':3}
kb=1.38064852e-23
rkb=1/kb
rolling_windows=2000
def get_temp_vs_z(file):
    title,atom,box= groio.parse_file(file)
    ###inital Dataframe
    df=pd.DataFrame(atom)
    del atom
    df=df.drop(columns=['resid','resname','atomid'])
    df=df[df['atom_name']!='MW' ]  # delete lone pair
    df['m']=(df['atom_name'].map(lambda x:mass_table[x[0:2]]))
    df['T']=(df.apply(lambda x: (x['vx']*x['vx']+x['vy']*x['vy']+x['vz']*x['vz'])*0.5*x['m']*1e6*2*rkb,1))
    df['T']=(df.apply(lambda x: (x['T']/deg_freedom[x['atom_name'][0:2]]),1))
    df=df[['z','T']]
    df.sort_values(by='z',inplace=True)
    df.index = pd.date_range('1/1/2000', periods=len(df), freq='s')
    df= (pd.DataFrame (df[['T','z']].resample('200s').mean() ))
    df.reset_index(drop=True,inplace=False)
    return df[['z','T']]


df = pd.DataFrame(columns=['z', 'T'])
files = os.listdir(path)
plt.figure(figsize=(5, 4))
for file in range(100,401,1):
    filepath=path+(str(file))+".gro"
    print(filepath)
    if len(df) == 0:
        df = get_temp_vs_z((filepath))
        dftemp=df

        continue
    data = get_temp_vs_z((filepath))
    df.index = data.index


    df['T'] = (df['T'] + data['T'])
    df['Ttemp']=df['T']

    df['Tavg']=df['Ttemp'].map(lambda x: x/(file-99))

    try:
        data.to_csv(path+str(file)+"Tempvsz.csv", columns=['z', 'T'])
        df.to_csv(path+"Tempvsz.csv",columns= ['z','Tavg'],index =False)
    except OSError:
        pass



