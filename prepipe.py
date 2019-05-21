#!/usr/bin/python

import sys, os
import vulcan_cfg
# for constructing the symbolic Jacobian matrix
from sympy import sin, cos, Matrix
from sympy import Symbol
from sympy import lambdify

ofname = 'chem_funs.py'
gibbs_text = vulcan_cfg.gibbs_text


# read the network and produce the .txt table for chemdf
# Re-arrange the numerbers in the network
def read_network():
    
    
    Rf, Rindx = {}, {}
    i = 1
    special_re = False

    ofstr = '# H/O/C chemistry \n\n'
    re_label = '#R'
    new_network = ''
       
    with open(vulcan_cfg.network, 'r') as f:
        for line in f.readlines():

            # switch to 3-body and dissociation reations 
            if line.startswith("# 3-body"): re_label = '#M'
            
            if line.startswith("# special"): 
                special_re = True # switch to reactions with special forms (hard coded)                   
                re_label = '#S'
            
            # skip common lines and blank lines
            # ========================================================================================
            if not line.startswith("#") and line.strip() and special_re == False: # if not starts
            
                Rf[i] = line.partition('[')[-1].rpartition(']')[0].strip()
                li = line.partition(']')[-1].strip()
                columns = li.split()
                #
                # if line.partition('[')[0].isspace():
                #     line = str(i) + line
                # else:
                #     Rindx[i] = int(line.partition('[')[0].strip())
                #     # updating the numerical index in the network (1, 3, ...)
                line = '{:4d}\t{:s}'.format(i, "".join(line.partition('[')[1:]))
                
                ofstr += re_label + str(i) + '\n'
                ofstr +=  Rf[i] + '\n'
            
                i += 2
            # ========================================================================================    
            elif special_re == True and line.strip() and not line.startswith("#"):

                #Rindx[i] = int(line.partition('[')[0].strip())
                Rf[i] = line.partition('[')[-1].rpartition(']')[0].strip()
                line = '{:4d}\t{:s}'.format(i, "".join(line.partition('[')[1:]))
                ofstr += re_label + str(i) + '\n'
                ofstr +=  Rf[i] + '\n'
        
                i += 2
            new_network += line
    
    with open(vulcan_cfg.network, 'w+') as f: f.write(new_network)
    
    return ofstr


def make_chemdf(re_table, ofname):
    '''
    make the file HOC.py cpnstaining the function chemdf for calculating the chemical production/loss term
    '''

    chem_dict = {}
    reac_dict = {}      # reaction dict. packed with v_i
    exp_reac_dict = {}  # explicit reaction dict. 
    i = -1
    j = -1 #index of forward reaction(odd number:1,3,5,...)
    count = 0 # count for even/odd term in v_exp

    reac_list = []
    rate_dict = {}

    for line in re_table.splitlines():
        '''
        reac : e.g. [[1, 'H'], [1, 'H'], [1, 'M']]
        reac_args: e.g. ['H2', 'H', 'M']
        '''
        #skip the space and #
        if line == '':
            continue
        elif line[0] == '#':
            continue
        else:
            reac = []
            mol_reac = []
            mol_prod = []
            prod = []
            rate_exp = ''
            rate_str = ""
            v_str = ""
            v_exp = ''  # to store the explicit expression of v_i (for jacobian)
            # R = True:reactants  R = False:products
            R = True    # if reads '->'
            N = False
        
            #Need to take car of M!!!
            for term in line.split():
                if term == '+':
                    continue
                elif term == "->":
                    # R=True:reactants R=False:products
                    R = False
                    continue
            
                # mol_list = [stoi-number, species name]    
                mol_list = term.split("*")
                if len(mol_list) == 1:
                    mol_list = [1] + mol_list
                else:
                    mol_list = [int(mol_list[0])] + mol_list[1:]
                mol = mol_list[1]
                stoi = int(mol_list[0])

                # check if the species is already included
                if not chem_dict.has_key(mol) and not term=='M':
                    i += 1
                    chem_dict.update({mol : i})
                    reac_dict.update({chem_dict[mol] : ""})
                    exp_reac_dict.update({chem_dict[mol] : ""})
                # if R is true, it's the reactants
                if R:
                    reac.append(mol_list)
                    mol_reac.append(mol)
                else:
                    prod.append(mol_list)
                    mol_prod.append(mol)

            j += 2
        
            reac_args = list(set(mol_reac + mol_prod)) #Remove repeating elements
            # because set exclude duplicates

            # v_i() is the rate equation function for i
            rate_str = "v_" + str(j) + "(k, M, "

            reac_noM = [ele for ele in reac if not ele[1]=='M' ]
            prod_noM = [ele for ele in prod if not ele[1]=='M' ]
            reac_args_noM = [ele for ele in reac_args if not ele=='M' ]
        
            for term in [ele for ele in reac_args if not ele=='M' ] :
                rate_str += "y[" + str(chem_dict[term]) + "], "  
            rate_str = rate_str[0:-2] + ")"
        
            v_exp += 'k[' + str(j) + ']*'
            for term in reac:
                if term[0]!=1:
                    if term[1]== 'M':
                        v_exp += term[1] + "**" + str(term[0]) + '*'
                    else:
                        v_exp += 'y[' + str(chem_dict[term[1]]) + ']' + "**" + str(term[0]) + '*'
                else:
                    if term[1]== 'M':
                        v_exp += term[1] + '*'
                    else:
                        v_exp += 'y[' + str(chem_dict[term[1]]) + ']*'              
            v_exp = v_exp[0:-1]  # Delete the last '*'

            v_exp += ' - k[' + str(j+1) + ']*'
            for term in prod:
                if term[0]!=1:
                    if term[1]== 'M':
                        v_exp += term[1] + "**" + str(term[0]) + '*'
                    else:
                        v_exp += 'y[' + str(chem_dict[term[1]]) + ']' + "**" + str(term[0]) + '*'
                else:
                    if term[1]== 'M':
                        v_exp += term[1] + '*'
                    else:
                        v_exp += 'y[' + str(chem_dict[term[1]]) + ']*'
            v_exp = v_exp[0:-1]  # Delete the last '*'
                
            for term in reac_noM:
                # term[0] os the stoi-number of the species
                reac_dict[chem_dict[term[1]]] += " -" + str(term[0]) + "*" + rate_str     
                if term[0]==1:
                    exp_reac_dict[chem_dict[term[1]]] += " -" + "(" + v_exp + ')'
                    count += 1
                else:
                    exp_reac_dict[chem_dict[term[1]]] += " -" + str(term[0]) + "*(" + v_exp + ')'
                    count += 1
            
            for term in prod_noM:
                reac_dict[chem_dict[term[1]]] += " +" + str(term[0]) + "*" + rate_str           
                if term[0]==1:
                    exp_reac_dict[chem_dict[term[1]]] += " +" + "(" + v_exp + ')'
                    count += 1
                else:
                    exp_reac_dict[chem_dict[term[1]]] += " +" + str(term[0]) + "*(" + v_exp + ')'
                    count += 1

            v_str = "#" + line + "\n"
            v_str += "v_" + str(j) + " = lambda k, M, "

            for term in reac_args_noM:
                v_str += term + ", "
            v_str = v_str[0:-2] + " : "
            v_str += 'k[' + str(j) + ']*'
            for term in reac:
                if term[0]!=1:
                    v_str += term[1] + "**" + str(term[0]) + '*'
                else:
                    v_str += term[1] + '*'
            v_str = v_str[0:-1] # Delete the last '*'

            v_str += ' - k[' + str(j+1) + ']*'
            for term in prod:
                if term[0]!=1:
                    v_str += term[1] + "**" + str(term[0]) + '*'
                else:
                    v_str += term[1] + '*'
            v_str = v_str[0:-1]

            reac_list.append(v_str)
        
            #ouput of each single rate from k1...
            rate_exp += 'k[' + str(j) + ']*'
            for term in reac:
                if term[1] == 'M':
                    rate_exp += 'M*'
                else:
                    if term[0]==1:
                        rate_exp += ('y['+str(chem_dict[term[1]])+']*' )
                    else:
                        rate_exp += ('y['+str(chem_dict[term[1]])+']**'+str(term[0])+'*' )
        
            rate_exp = rate_exp[0:-1]  # Delete the last '*'
            rate_dict[j] = rate_exp
            rate_exp = ''
            rate_exp += 'k[' + str(j+1) + ']*' #j+1 even number for reverse index
            for term in prod:
                if term[1] == 'M':
                    rate_exp += 'M*'
                else:
                    if term[0]==1:
                        rate_exp += ('y['+str(chem_dict[term[1]])+']*' )
                    else:
                        rate_exp += ('y['+str(chem_dict[term[1]])+']**'+str(term[0])+'*' )
            rate_exp = rate_exp[0:-1]  # Delete the last '*'
            rate_dict[j+1] = rate_exp 


    # save output
    chem_dict_r = {}
    spec_list = []
        
    ofstr = "#!/usr/bin/python\n\nfrom scipy import *\nimport numpy as np\nfrom phy_const import kb, Navo\nimport vulcan_cfg\n\n"
    ofstr += "'''\n## Reaction ##\n\n"
    ofstr += re_table + "\n\n"

    ofstr += "## Mapping ##\n\n"
    for term in chem_dict:
        chem_dict_r.update({chem_dict[term] : term})
    for term in reac_dict:
        ofstr += chem_dict_r[term] + ': y[' + str(term) + '], '
    ofstr+='\n\n'
    for term in reac_dict:
        ofstr += chem_dict_r[term] + "\t" + str(term) + "\t" + reac_dict[term] + "\n"
    for i in chem_dict_r:
        spec_list.append(chem_dict_r[i])
  
    ofstr += "'''\n\n"
    
    ofstr += '#species list\n'
    ofstr += 'spec_list = ' + str(spec_list)
    ofstr += '\n# the total number of species'
    ofstr += '\nni = ' + str(len(chem_dict.keys()))
    ofstr += '\n# the total number of reactions (forward and reverse)'
    ofstr += '\nnr = ' + str(len(rate_dict.keys()))
    
    ost = '\n\ndef chemdf(y, M, k): \n' # Note: making M as input!!!
    ost += '\t y = np.transpose(y) \n'.expandtabs(3)
    ost += '\t dydt = np.zeros(shape=y.shape) \n'.expandtabs(3)

    for num in reac_dict:
        ost += '\t dydt['.expandtabs(3) + str(num) + '] = ' + reac_dict[num] + '\n'

    ost += '\t dydt = np.transpose(dydt) \n'.expandtabs(3)
    ost += '\t return dydt \n\n'.expandtabs(3)

    ost += 'def df(y, M, k):\n'
    ost += '\t df_list = [] \n'.expandtabs(3)
    for num in exp_reac_dict:
        ost += '\t df_list.append( '.expandtabs(3) +exp_reac_dict[num] + ' )\n'    
    ost += '\t return df_list \n\n'.expandtabs(3)
    ofstr += ost

    for term in reac_list:
        ofstr += term + "\n\n"

    with open(ofname, "w") as of:
        of.write(ofstr)
        
    # return (ni, nr)
    return (len(chem_dict.keys()), len(rate_dict.keys()))
        

               
def make_Gibbs(re_table, gibbs_text, ofname): 
    '''
    Calculating the equilibrium constants (K_eq) from the Gibbs free energy to reverse the reaction rates.
    To DO: combine the repetitive parts of make_chemdf and make_Gibbs into one finction
    '''
    chem_dict = {}
    reac_dict = {}      # reaction dict. packed with v_i
    exp_reac_dict = {}  # explicit reaction dict. 
    i = -1
    j = -1 #index of forward reaction(odd number:1,3,5,...)
    count = 0 # count for even/odd term in v_exp
    reac_list = []
    rate_dict = {}
    
    reac_list, rate_dict, gibbs_dict = [], {}, {}
    gstr = ''
    
    for line in re_table.splitlines():
        if line == '':
            continue
        elif line[0] == '#':
            continue
        else:
            reac = []
            mol_reac = []
            mol_prod = []
            prod = []
            rate_exp = ''
            rate_str = ""
            v_str = ""
            v_exp = ''  # to store the explicit expression of v_i (for jacobian)
            gibbs = 'np.exp( -('
            R = True
            N = False
            reac_num, prod_num = 0, 0
        
            #Need to take car of M!!!
            for term in line.split():
                #print 'term' + term
                if term == '+':
                    continue
                elif term == "->":
                    # R=True:reactants R=False:products
                    R = False
                    continue
            
                # mol_list = [stoi-number, species name]    
                mol_list = term.split("*")
                if len(mol_list) == 1:
                    mol_list = [1] + mol_list
                else:
                    mol_list = [int(mol_list[0])] + mol_list[1:]
                mol = mol_list[1]
                stoi = int(mol_list[0])
            

                # check if the species is already included
                if not chem_dict.has_key(mol) and not term=='M':
                    i += 1
                    chem_dict.update({mol : i})
                    reac_dict.update({chem_dict[mol] : ""})
                    exp_reac_dict.update({chem_dict[mol] : ""})
                # if R is true, it's the reactants
                if R:
                    reac.append(mol_list)
                    mol_reac.append(mol)
                else:
                    prod.append(mol_list)
                    mol_prod.append(mol)


            j += 2
        
            reac_args = list(set(mol_reac + mol_prod)) #Remove repeating elements
            # because set exclude duplicates

            # v_i() is the rate equation function for i
            rate_str = "v_" + str(j) + "("

            reac_noM = [ele for ele in reac if not ele[1]=='M' ]
            prod_noM = [ele for ele in prod if not ele[1]=='M' ]
            reac_args_noM = [ele for ele in reac_args if not ele=='M' ]
 
            for term in [ele for ele in reac_args if not ele=='M' ] :
                rate_str += "y[" + str(chem_dict[term]) + "], "  
            rate_str = rate_str[0:-2] + ")"
        
            v_exp += 'k[' + str(j) + ']*'
            for term in reac:
                if term[0]!=1:
                    if term[1]== 'M':
                        v_exp += term[1] + "**" + str(term[0]) + '*'
                    else:
                        v_exp += 'y[' + str(chem_dict[term[1]]) + ']' + "**" + str(term[0]) + '*'
                else:
                    if term[1]== 'M':
                        v_exp += term[1] + '*'
                    else:
                        v_exp += 'y[' + str(chem_dict[term[1]]) + ']*'              
            v_exp = v_exp[0:-1]  # Delete the last '*'

            v_exp += ' - k[' + str(j+1) + ']*'
            for term in prod:
                if term[0]!=1:
                    if term[1]== 'M':
                        v_exp += term[1] + "**" + str(term[0]) + '*'
                    else:
                        v_exp += 'y[' + str(chem_dict[term[1]]) + ']' + "**" + str(term[0]) + '*'
                else:
                    if term[1]== 'M':
                        v_exp += term[1] + '*'
                    else:
                        v_exp += 'y[' + str(chem_dict[term[1]]) + ']*'
            v_exp = v_exp[0:-1]  # Delete the last '*'
                
            for term in reac_noM:
                # term[0] is the stoi-number of the species
                reac_dict[chem_dict[term[1]]] += " -" + str(term[0]) + "*" + rate_str     
                if term[0]==1:
                    exp_reac_dict[chem_dict[term[1]]] += " -" + "(" + v_exp + ')'
                    count += 1
                else:
                    exp_reac_dict[chem_dict[term[1]]] += " -" + str(term[0]) + "*(" + v_exp + ')'
                    count += 1
            
            for term in prod_noM:
                reac_dict[chem_dict[term[1]]] += " +" + str(term[0]) + "*" + rate_str           
                if term[0]==1:
                    exp_reac_dict[chem_dict[term[1]]] += " +" + "(" + v_exp + ')'
                    count += 1
                else:
                    exp_reac_dict[chem_dict[term[1]]] += " +" + str(term[0]) + "*(" + v_exp + ')'
                    count += 1
                
            ######################## for constructing Gibbs free energy ########################  
            for term in reac_noM:
                gibbs += '-' + str(term[0]) + '*' + "gibbs_sp('" + str(term[1]) +"',T)"
                reac_num += term[0]
            
            for term in prod_noM:
                gibbs += '+' + str(term[0]) + '*' + "gibbs_sp('" + str(term[1]) +"',T)"
                prod_num += term[0]
        
            gibbs += ' ) )'
            if prod_num - reac_num != 0:
                gibbs += '*(corr*T)**' + str(reac_num - prod_num)   
            gibbs_dict[j] = gibbs 
            ######################## for constructing Gibbs free energy ######################## 

            v_str = "#" + line + "\n"
            v_str += "v_" + str(j) + " = lambda "

            for term in reac_args_noM:
                v_str += term + ", "
            v_str = v_str[0:-2] + " : "
            # j: ->
            # j+1: <-
            v_str += 'k[' + str(j) + ']*'
            for term in reac:
                if term[0]!=1:
                    v_str += term[1] + "**" + str(term[0]) + '*'
                else:
                    v_str += term[1] + '*'
            v_str = v_str[0:-1] # Delete the last '*'

            v_str += ' - k[' + str(j+1) + ']*'
            for term in prod:
                if term[0]!=1:
                    v_str += term[1] + "**" + str(term[0]) + '*'
                else:
                    v_str += term[1] + '*'
            v_str = v_str[0:-1]

            reac_list.append(v_str)
        
            #ouput of each single rate from k1...
            rate_exp += 'k[' + str(j) + ']*'
            for term in reac:
                if term[1] == 'M':
                    rate_exp += 'M*'
                else:
                    if term[0]==1:
                        rate_exp += ('y['+str(chem_dict[term[1]])+']*' )
                    else:
                        rate_exp += ('y['+str(chem_dict[term[1]])+']**'+str(term[0])+'*' )
        
            rate_exp = rate_exp[0:-1]  # Delete the last '*'
            rate_dict[j] = rate_exp
            rate_exp = ''
            rate_exp += 'k[' + str(j+1) + ']*' #j+1 even number for reverse index
            for term in prod:
                if term[1] == 'M':
                    rate_exp += 'M*'
                else:
                    if term[0]==1:
                        rate_exp += ('y['+str(chem_dict[term[1]])+']*' )
                    else:
                        rate_exp += ('y['+str(chem_dict[term[1]])+']**'+str(term[0])+'*' )
            rate_exp = rate_exp[0:-1]  # Delete the last '*'
            rate_dict[j+1] = rate_exp 

    
    with open(gibbs_text) as g:
        for line in g:
            gstr += line

    gstr += '\n\n'
    gstr += '# Gibbs free energy:\n'
    gstr += 'def Gibbs(i,T):\n'
    gstr += '\t G={}\n'.expandtabs(3)
    for _ in gibbs_dict:
        gstr += '\t G['.expandtabs(3) +str(_)+'] = lambda T: ' + str(gibbs_dict[_]) + '\n'

    gstr += '\t return G[i](T) \n\n'.expandtabs(3)

    with open(ofname, 'a+') as f:
        f.write(gstr)
        

def make_jac(ni, nr, ofname):
    '''
    to make the analytical Jocobian matrix of chemdf
    '''
    M = Symbol('M')
    y, k = [], []

    for i in range(ni):
        y.append( Symbol('y[:,'+str(i)+']') )
    for i in range(nr+1):
        k.append( Symbol('k['+str(i)+']') )
    
    # chemistry is the "ofname" module
    dy = Matrix(chemistry.df(y,M,k))
    x = Matrix(y)
    jac = dy.jacobian(x)

    jstr = '\ndef symjac(y, M, k): \n'
    jstr += '\t nz = vulcan_cfg.nz\n'.expandtabs(3)
    jstr += '\t dfdy = np.zeros(shape=[ni*nz, ni*nz])   \n'.expandtabs(3)
    jstr += '\t indx = [] \n'.expandtabs(3)
    jstr += '\t for j in range(ni): \n'.expandtabs(3)
    jstr += '\t indx.append( np.arange(j,j+ni*nz,ni) ) \n'.expandtabs(7)

    for i in range(ni):
        for j in range(ni):
            jstr += '\t dfdy[indx['.expandtabs(3) + str(i) + '], indx[' + str(j) +']] = ' + str(jac[i,j]) + '\n'

    jstr += '\t return dfdy \n\n'.expandtabs(3)

    # save the output function
    with open (ofname, 'a+') as f: f.write(jstr)
    

if __name__ == "__main__":   
    re_table = read_network()
    (ni, nr) = make_chemdf(re_table, ofname)
    make_Gibbs(re_table, gibbs_text, ofname)
    # import the "ofname" module as chemistry for make_jac to read df
    chemistry = __import__(ofname[:-3])
    make_jac(ni, nr, ofname)