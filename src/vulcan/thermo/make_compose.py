import numpy as np
com_file = 'nomass_all_compose.txt'


compo = np.genfromtxt(com_file,names=True,dtype=None)
compo_row = list(compo['species'])

ost = 'species			H	O	C	He	N	S	P	Na	K	Si 	mass      \n'
with open(com_file) as f:
    for line in f.readlines():                

        if not line.startswith("species") and line.strip(): 

            col = line.split()
            indx = compo_row.index(col[0])
            
            # check sulfer
            mass = compo[indx][1]*1.008 + compo[indx][2]*16. + compo[indx][3]*12.011 + compo[indx][4]*4.003\
            + compo[indx][5]*14.007 + compo[indx][6]*32.065 + compo[indx][7]*30.974 +  compo[indx][8]*22.990\
            + compo[indx][9]*39.10 + compo[indx][10]*28.085

            line = line[:-1]
            line = line.rstrip() + '\t' + str(mass) + '\n' 
            ost += line
            

with open('all_compose.txt', 'w+') as f: f.write(ost)