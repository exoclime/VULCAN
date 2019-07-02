import numpy as np 
from collections import defaultdict
from scipy import interpolate
import itertools
import operator
import pydot
import os, sys
from PIL import Image 
import pickle
import math

# import VULCAN modules
import store, op
import vulcan_cfg
from phy_const import kb, Navo
import chem_funs

# The chemical network in vulcan_cfg.py needs to be the same as the one used in the vulcan ouput file (vul_data)
# no arguments or not setting '-n' (no re-making chem_funs.py) option
if len(sys.argv) < 2 or sys.argv[1] != '-n': 
    # running prepipe to construch chem_funs.py
    print ('Making chem_funs.py ...')
    python_executable = sys.executable
    os.system(python_executable + ' make_chem_funs.py')
else: pass

# the pressure level to analise (in cgs)
p_ana_list = [1e0]

vul_data = 'output/HD189.vul'
with open(vul_data, 'rb') as handle:
  data = pickle.load(handle)
vulcan_spec = data['variable']['species']

### read in the basic chemistry data
with open(vulcan_cfg.com_file, 'r') as f:
    columns = f.readline() # reading in the first line
    num_ele = len(columns.split())-2 # number of elements (-2 for removing "species" and "mass") 
type_list = ['int' for i in range(num_ele)]
type_list.insert(0,'U20'); type_list.append('float')
compo = np.genfromtxt(vulcan_cfg.com_file,names=True,dtype=type_list)
# dtype=None in python 2.X but Sx -> Ux in python3
compo_row = list(compo['species'])

species = chem_funs.spec_list
re_dict = chem_funs.re_dict
re_wM_dict = chem_funs.re_wM_dict
np.set_printoptions(threshold=np.inf)

# the pathway (conversion of the 2 species) to analise
conv_sp = ('N2', 'HCN')
# Chooseing the common element of the two species to track, e.g. for CH4 -> CO, one should be interested in the flow of C
atom = 'N'
# species to analize
tmp_list, ana_sp = list(species), []
for sp in tmp_list:
    if not compo[compo_row.index(sp)][atom] == 0: 
        ana_sp.append(sp)
print ('Including species with ' + str(atom) + ':\n' + str(ana_sp))

# all the combination from ana_sp
ana_comb = list(itertools.combinations(ana_sp, 2))

remove_fw_re = []
remove_re = list(remove_fw_re)
[remove_re.append(re+1) for re in remove_fw_re]
remove_pair = []
for pair in remove_pair:  
    if pair in ana_comb: ana_comb.remove(pair)
    if pair[::-1] in ana_comb: ana_comb.remove(pair[::-1])
        
for p_ana in p_ana_list:
    # Find the index of pco closest to p_ana
    p_indx = min( range(len(data['atm']['pco'])), key=lambda i: abs(data['atm']['pco'][i]-p_ana))
    # pressure and T from the vulcan data
    pp = data['atm']['pco'][p_indx]
    T_ana = data['atm']['Tco'][p_indx]
    # the third body
    mm = p_ana*1.e6/kb/T_ana

    # store the combination of prod and reac
    re_comb = {}

    # store all the connections
    edges = defaultdict(list)

    # store the net reaction rates of a pair of species, and the fastest reaction channel
    netRate_pair, maxRate = {}, {}
    # store the contribution (percentage) of the fastest reaction between the two species
    max_contri = {}
    # initializing zero
    for pair in ana_comb: 
        netRate_pair[pair] = 0.
        # reverse
        netRate_pair[pair[::-1] ] = 0.
        # first element stores the rate and the 2md element stores the reaction number 
        maxRate[pair] = [0. , np.nan]
        maxRate[pair[::-1]] = [0., np.nan]

        edges[pair[0]].append(pair[1]) 
        edges[pair[1]].append(pair[0])
        edges[pair[0]] = list(set(edges[pair[0]]))
        edges[pair[1]] = list(set(edges[pair[1]]))

    # Setting the current working directory to the script location
    abspath = os.path.abspath(sys.argv[0])
    dname = os.path.dirname(abspath)
    os.chdir(dname)

    from chem_funs import ni, nr  # number of species and reactions in the network
    np.set_printoptions(threshold='nan')  # print all for debuging


    # Reading rates
    re_prod, re_reat = {}, {}
    for re in range(1,nr+1): 
        if not re in remove_re:
            re_prod[re] = [prod for prod in re_dict[re][0] if prod in ana_sp]
            re_reat[re] = [reat for reat in re_dict[re][1] if reat in ana_sp]

            comb = []
            for prod in re_prod[re]:
                for reat in re_reat[re]:
                    comb.append((prod, reat)) #tuple
            re_comb[re] = comb

            # collecting rates (not rate constants) for each combination in ana_comb
            for pair in re_comb[re]:
                if pair in ana_comb or pair[::-1] in ana_comb:
                    # for the reverse of the photodissociation
                    #print (re)
                    if np.all(data['variable']['k'][re]== 0): 
                        #print (str(re) + '  passs')
                        pass
                    
                    else:
                        rate_const = data['variable']['k'][re][p_indx]
                        this_rate = rate_const
                        for sp in re_wM_dict[re][0]: # Looping all the "reactants" including M

                            if not sp == 'M': 
                                this_rate *= float(data['variable']['y'][p_indx,vulcan_spec.index(sp)]) 
                    
                            elif sp == 'M': this_rate *= mm
                            else: raise IOError ('\nUnknow species name in the network.')

                        if np.abs(this_rate) > maxRate[pair][0]: 
                            maxRate[pair][0] = np.abs(this_rate)
                            maxRate[pair][1] = re
                        netRate_pair[pair] +=  this_rate
                        
    # Calculating the contribution of the main reactions in maxRate
    for pair in ana_comb:
        if netRate_pair[pair] == 0:
            max_contri[pair] = np.nan
        elif netRate_pair[pair[::-1] ] == 0:
            max_contri[pair[::-1] ] = np.nan
        else:
            max_contri[pair] = maxRate[pair][0]/netRate_pair[pair]
            max_contri[pair[::-1] ] = maxRate[pair[::-1] ][0]/netRate_pair[pair[::-1] ]



    # Graph Data
    '''
    ana_sp = nodes  
    edges_dict = edges
    maxRate = distance
    '''

    class Graph:
      def __init__(self):
        self.nodes = ana_sp
        self.edges = edges
        self.distance = maxRate
        
    g = Graph()

    def shortest_path(graph, ini, end, skip):
        edges = graph.edges.copy()
        distance = graph.distance.copy()

        # Initialization
        distance_indx = defaultdict(lambda:np.inf)
        distance_indx[ini] = 0
        unvisited_distance = {node: np.inf for node in list(graph.nodes)}
        visited, path, path_re = [], [], []
        prev = {}
        nodes = set(graph.nodes)
        # the current node
        now = ini

        while unvisited_distance: # while unvisited_distance is not empty
            for neighbor in edges[now]:
                if neighbor in unvisited_distance.keys(): # never check the visited nodes again!
                    rate = distance[(now, neighbor)][0]
                    if not rate == 0 and end not in re_prod[ distance[(now, neighbor)][1] ] and distance[(now, neighbor)][1] not in skip\
                    and (ini not in re_prod[ distance[(now, neighbor)][1] ] and ini not in re_reat[ distance[(now, neighbor)][1] ] or now == ini): # meaning it is connected (otherewise rate would be 0)
                    # and the end product is not allowed to participate the reaction

                        temp_distance = distance_indx[now] + 1./rate 

                        if temp_distance < distance_indx[neighbor]: 
                            # distance_indx is only for recording footprints; unvisited_distance is for the algorithm 
                            prev[neighbor] = now # recording the pathway
                            distance_indx[neighbor] = temp_distance  
                            unvisited_distance[neighbor] = temp_distance 
    
            # marking the node "JUST VISITED"
            visited.append(now) 
            del unvisited_distance[now]

            # select the node that is marked with the smallest distance_indx from ALL unvisited nodes
            next_now = min(unvisited_distance,key=unvisited_distance.get) # key = xx.get gives the access to the values
            now = next_now
            
            # check if stops 
            if now == end:
                # visited sotres all the "footprints" including all attemps
                visited.append(now)
                path.append(now)
                goback = now
                while not goback==ini:
                    path.append(prev[goback])
                    path_re.append(maxRate[(goback,prev[goback])][1])
                    goback = prev[goback]
                return (path[::-1], path_re[::-1], distance_indx)
            elif unvisited_distance[next_now] == np.inf:
                print ('Not connected!')
                break
        

    # Plotting     
    # Creating path graph
    gdot = pydot.Dot(graph_type='digraph',  center='True') # graph_type='graph' for non-directional

    # to weed out taking the same reaction twice
    # e.g. R1: A + B -> C + D 
    # A (R1)-> C (R1)-> B is not allowed

    # also exclude the reaction that produces the straing species again
    # e.g. Path A -> B -> C 
    # the dominant reaction for B -> C to be B + X -> C + "A" should not be allowed!!!


    repet_re = []
    loop_indx = 1
    while True:
        round_trip = False
        path, path_re, time_lables = shortest_path(g, conv_sp[0], conv_sp[1], repet_re)
        print ('Try ' + str(loop_indx) + ':' )
        print (path_re)
        for re in path_re:
            if re % 2 ==1: # only checking the forward reactions is enough
                if re+1 in path_re:
                    # when a round trip has been taken i.e. R1 forward and R2 back
                    # Try remove R1 and keep R2, then remove R2 and keep R1
                    # If still casues a loop, remove both R1 and R2
                    
                    temp_skip = repet_re.copy()
                    path_1, path_re_1, time_lables_1 = shortest_path(g, conv_sp[0], conv_sp[1], temp_skip + [re])
                    path_2, path_re_2, time_lables_2 = shortest_path(g, conv_sp[0], conv_sp[1], temp_skip + [re+1])
                    
                    for re in path_re_1:
                        if re % 2 ==1: 
                            if re+1 in path_re_1: rm_1 = False
                            else: rm_1 = True
                            
                    for re in path_re_2:
                        if re % 2 ==1: 
                            if re+1 in path_re_2: rm_2 = False
                            else: rm_2 = True
                    
                    if path_re_1 and path_re_2:
                        min_rate_1, min_rate_2 = np.inf, np.inf
                        max_rate_1, max_rate_2 = 0., 0.
                        for sp in path_1:
                            if not sp==path_1[-1]:
                                pair = (sp, (path_1[path_1.index(sp)+1]))
                                if maxRate[pair][0] < min_rate_1:
                                    min_rate_1 = maxRate[pair][0]
                        for sp in path_2:
                            if not sp==path_2[-1]:
                                pair = (sp, (path_2[path_2.index(sp)+1]))
                                if maxRate[pair][0] < min_rate_2:
                                    min_rate_2 = maxRate[pair][0]
                        
                        if max_rate_1 > max_rate_2: path, path_re, time_lables = path_1, path_re_1, time_lables_1
                        else:  path, path_re, time_lables = path_2, path_re_2, time_lables_2
                        break
                        
                    elif path_re1:
                        path, path_re, time_lables = path_1, path_re_1, time_lables_1
                        break
                    elif path_re2:
                        path, path_re, time_lables = path_2, path_re_2, time_lables_2
                        break
                    else:
                        repet_re.append(re)
                        repet_re.append(re+1)
                        round_trip = True
            
        if round_trip == False: break
        
    pathway = {}
    min_rate = np.inf
    max_rate = 0.
    for sp in path:
        if not sp==path[-1]:
            #print ( (sp, (path[path.index(sp)+1]) ) )
            pair = (sp, (path[path.index(sp)+1]))
            pathway[pair] = maxRate[pair]
            if maxRate[pair][0] < min_rate:
                min_rate = maxRate[pair][0]
                rls = pair
                rls_re = maxRate[pair][1]
            
            # for scaling the plot
            if maxRate[pair][0] > max_rate:
                max_rate = maxRate[pair][0] 


    # maximum and minumum rates used to scale the edge width
    w_min = np.log10(min_rate)
    w_max = np.log10(max_rate)


    for sp in path:
        if sp == 'C4H2': sp = 'Haze'
        node = pydot.Node(sp, color="k",fontcolor='k',style='bold')
        gdot.add_node(node)

    # RLS
    if rls_re % 2 == 0: rls_re -= 1
    re_exp = ' ' + data['variable']['Rf'][rls_re] + ' '
    re_exp = re_exp.replace("->", "=")
    
    if 'C4H2' in rls:
        edge = pydot.Edge( rls[0], 'Haze', label="{:0.2E}".format(pathway[pair][0]) + '\nR' + str(pathway[pair][1])\
        + ' ('+"{:0.2f}".format(max_contri[pair]*100.) +'%)\n' + re_exp, penwidth=1, color='red', fontsize=13, fontcolor='red')
    else:
        edge = pydot.Edge( rls[0], rls[1], label="{:0.2E}".format(pathway[rls][0]) + '\nR' + str(pathway[rls][1])\
        + ' ('+"{:0.2f}".format(max_contri[rls]*100.) +'%)\n' + re_exp, penwidth=1, color='red', fontsize=13, fontcolor='red')
    gdot.add_edge(edge)

    for pair in pathway.keys():
        
        if pair is not rls:
            lwidth = 1. + 4.5*(np.log10(pathway[pair][0]) - w_min)/(w_max-w_min)
            main_re = pathway[pair][1]
            if main_re % 2 == 0: main_re -= 1
            re_exp = ' ' + data['variable']['Rf'][main_re] + " "
            re_exp = re_exp.replace("->", "=")
            
            # changing the printing of photodissociation re
            if main_re == 585: re_exp = 'CH4 + hv -> CH3 + H'
            
            if 'C4H2' in pair:
                edge = pydot.Edge( pair[0], 'Haze', label="{:0.2E}".format(pathway[pair][0]) + '\nR' + str(pathway[pair][1])\
                + ' ('+"{:0.2f}".format(max_contri[pair]*100.) +'%)\n' + re_exp, penwidth=lwidth, color='k', fontsize=12, fontcolor='k', nodecolor='red')
            else:
                edge = pydot.Edge( pair[0], pair[1], label="{:0.2E}".format(pathway[pair][0]) + '\nR' + str(pathway[pair][1])\
            + ' ('+"{:0.2f}".format(max_contri[pair]*100.) +'%)\n' + re_exp, penwidth=lwidth, color='k', fontsize=12, fontcolor='k', nodecolor='red')

            gdot.add_edge(edge)

    
    if not os.path.exists('plot/pathway/'):
        print ('The directory for pathway plots does not exist.')
        print( 'Directory ' , 'plot/pathway/',  " created.")
        os.mkdir('plot/pathway/')
        
    outfile = 'plot/pathway/dynamic_' + conv_sp[0]+"-"+conv_sp[1]+'_' + "{:0.1E}".format(p_ana) + '_bar_T' + str(int(T_ana))+'_path.png'
    gdot.write_png(outfile)
    plot = Image.open(outfile)
    plot.show()
    gdot.write_pdf('plot/pathway/dynamic_' + conv_sp[0]+"-"+conv_sp[1]+'_' + "{:0.1E}".format(p_ana) + '_bar_T' + str(int(T_ana))+'_path.pdf')

