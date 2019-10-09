import groio
import pandas as pd

title,atom,box= groio.parse_file("./data/0.gro")
NA=6.02*pow(10,23)
rNA=1/NA
mass_table={'OW':16*rNA,
            'HW':1*rNA,
            'C':12*rNA,
            'N':14*rNA,
            'MW':0*rNA,
            'CL':35.5*rNA,
            'NA':23*rNA}
kb=1.38064852*pow(10,(-24))
rkb=1/kb
dz=0.1  #nm
###inital Dataframe
df=pd.DataFrame(atom)
df['m']=df.apply(lambda x:mass_table[df['atom_name'][1][0:2]], axis=1)
df['T']=df.apply(lambda x:df['vx']*df['vx']+df['vy']*df['vy']+df['vz']*df['vz']*0.5*df['m']*pow(1,-6)*(2/3)*rkb,1)
###   T=(vx^2+vy^2+vz^2)*0.5*m*1e-6*(2/3)*rkb
print (df)