import numpy as np
import scipy
from scipy import sparse
import matplotlib.pyplot as plt
import matplotlib.legend as lg
import time, os, pickle

import vulcan_cfg
try: from PIL import Image
except ImportError: 
    try: import Image
    except: vulcan_cfg.use_PIL = False

import build_atm
import chem_funs
from chem_funs import ni, nr  # number of species and reactions in the network
from phy_const import kb, Navo
from vulcan_cfg import nz

# imported functions 
chemdf = chem_funs.chemdf
achemjac = chem_funs.symjac
compo = build_atm.compo
compo_row = build_atm.compo_row

species = chem_funs.spec_list


class ReadRate(object):
    
    """
    to read in rate constants from the network file and compute the reaction rates for the corresponding Tco and pco 
    """
    
    def __init__(self):
        
        self.i = 1
        # flag of trimolecular reaction
        self.re_tri, self.re_tri_k0 = False, False
        self.list_tri = []

        
    def read_rate(self, var, atm, build_table=False):
        
        Rf, Rindx, a, n, E, a_inf, n_inf, E_inf, k, k_fun, k_inf, kinf_fun, k_fun_new = \
        var.Rf, var.Rindx, var.a, var.n, var.E, var.a_inf, var.n_inf, var.E_inf, var.k, var.k_fun, var.k_inf,  var.kinf_fun,  var.k_fun_new
        
        i = self.i
        # flags for starting 3-body, 3-body without high-P rates
        re_tri, re_tri_k0 = self.re_tri, self.re_tri_k0
        list_tri = self.list_tri
        
        Tco = atm.Tco.copy()
        M = atm.M.copy()
        
        special_re = False
               
        with open(vulcan_cfg.network) as f:
            for line in f.readlines():
                
                # switch to 3-body and dissociation reations 
                if line.startswith("# 3-body"): 
                    re_tri = True
                
                # switch to 3-body without high-pressure rates
                if line.startswith("# 3-body reactions without high-pressure rates"):
                    re_tri_k0 = True
                    
                if line.startswith("# special"): 
                    special_re = True # switch to reactions with special forms (hard coded)                   
                    
                # skip common lines and blank lines
                # ========================================================================================
                if not line.startswith("#") and line.strip() and special_re == False: # if not starts
                    
                    Rf[i] = line.partition('[')[-1].rpartition(']')[0].strip()
                    li = line.partition(']')[-1].strip()
                    columns = li.split()
                    Rindx[i] = int(line.partition('[')[0].strip())
                    
                    a[i] = float(columns[0])
                    n[i] = float(columns[1])
                    E[i] = float(columns[2])
                    
                    if build_table == True:
                        ofstr += re_label + str(i) + '\n'
                        ofstr +=  Rf[i] + '\n'
                
                    # switching to trimolecular reactions (re_tri_k0 == False for those with high-P limit rates)   
                    if re_tri == True and re_tri_k0 == False:
                        a_inf[i] = float(columns[3])
                        n_inf[i] = float(columns[4])
                        E_inf[i] = float(columns[5])
                        list_tri.append(i) 
                    
                    if columns[-1].strip() == 'He': re_He = i
                    elif columns[-1].strip() == 'ex1': re_CH3OH = i
                
                    # Note: make the defaut i=i
                    k_fun[i] = lambda temp, mm, i=i: a[i] *temp**n[i] * np.exp(-E[i]/temp)
                
                
                    if re_tri == False:
                        k[i] = k_fun[i](Tco, M)
                    
                    # for 3-body reactions, also calculating k_inf
                    elif re_tri == True and len(columns)>=6:
        
        
                        kinf_fun[i] = lambda temp, i=i: a_inf[i] *temp**n_inf[i] * np.exp(-E_inf[i]/temp)
                        k_fun_new[i] = lambda temp, mm, i=i: (a[i] *temp**n[i] * np.exp(-E[i]/temp))/(1 + (a[i] *temp**n[i] * np.exp(-E[i]/temp))*mm/(a_inf[i] *temp**n_inf[i] * np.exp(-E_inf[i]/temp)) ) 
        
                        #k[i] = k_fun_new[i](Tco, M)
                        k_inf = a_inf[i] *Tco**n_inf[i] * np.exp(-E_inf[i]/Tco)
                        k[i] = k_fun[i](Tco, M)
                        k[i] = k[i]/(1 + k[i]*M/k_inf )
        
        
                    else: # for 3-body reactions without high-pressure rates
                        k[i] = k_fun[i](Tco, M)
                                
                    ### TEST CAPPING
                    # k[i] = np.minimum(k[i],1.E-11)
                    ###
    
                    i += 2
                    # end if not 
                 # ========================================================================================    
                elif special_re == True and line.strip() and not line.startswith("#"):

                    Rindx[i] = int(line.partition('[')[0].strip())
                    Rf[i] = line.partition('[')[-1].rpartition(']')[0].strip()
                
                    if Rf[i] == 'OH + CH3 + M -> CH3OH + M':
                        print ('Using special form for the reaction: ' + Rf[i])
                    
                        k[i] = 1.932E3*Tco**-9.88 *np.exp(-7544./Tco) + 5.109E-11*Tco**-6.25 *np.exp(-1433./Tco)
                        k_inf = 1.031E-10 * Tco**-0.018 *np.exp(16.74/Tco)
                        k[i] = k[i]/(1 + k[i]*M/k_inf )
                    
                        k_fun[i] = lambda temp, mm, i=i: 1.932E3 *temp**-9.88 *np.exp(-7544./temp) + 5.109E-11*temp**-6.25 *np.exp(-1433./temp)
                        kinf_fun[i] = lambda temp, mm, i=i: 1.031E-10 * temp**-0.018 *np.exp(16.74/temp)
                        k_fun_new[i] = lambda temp, mm, i=i: (1.932E3 *temp**-9.88 *np.exp(-7544./temp) + 5.109E-11*temp**-6.25 *np.exp(-1433./temp))/\
                        (1 + (1.932E3 *temp**-9.88 *np.exp(-7544./temp) + 5.109E-11*temp**-6.25 *np.exp(-1433./temp)) * mm / (1.031E-10 * temp**-0.018 *np.exp(16.74/temp)) )
                    
                    i += 2

        k_fun.update(k_fun_new)
    
        # store k into data_var
        # remeber k_fun has not removed reactions from remove_list
        var.k = k
        var.k_fun = k_fun
        var.kinf_fun = kinf_fun

        return var
        
    def rev_rate(self, var, atm):
        
        rev_list = range(2,nr+1,2)
        Tco = atm.Tco.copy()
        
        # reversing rates and storing into data_var
        for i in rev_list: 
            var.k_fun[i] = lambda temp, mm, i=i: var.k_fun[i-1](temp, mm)/chem_funs.Gibbs(i-1,temp)
            var.k[i] = var.k[i-1]/chem_funs.Gibbs(i-1,Tco)
       
        return var
        
    
    def remove_rate(self, var):
        
        for i in vulcan_cfg.remove_list:
            var.k[i] = np.repeat(0,nz)
            var.k_fun[i] = lambda temp, mm, i=i: np.repeat(0,nz)
            
        return var
        

class Integration(object):
    """
    time-stepping until the stopping criteria (steady-state) is satisfied
    #all the operators required in the continuity equation: dn/dt + dphi/dz = P - L
    #or class incorporating the esential numerical operations?
    """
    
    def __init__(self, odesolver, output):

        self.mtol = vulcan_cfg.mtol
        self.atol = vulcan_cfg.atol
        self.output = output
         
        self.odesolver = odesolver
        
    def __call__(self, var, atm, para, make_atm):
        
        use_print_prog, use_live_plot = vulcan_cfg.use_print_prog, vulcan_cfg.use_live_plot
        
        while not self.stop(var, para): # Looping until the stop condition is satisfied
            
            var = self.backup(var)
            var, para = self.odesolver.one_step(var, atm, para)
           
            if para.count % vulcan_cfg.update_frq == 0: # updating mu and dz (dzi) due to diffusion
                atm = self.update_mu_dz(var, atm, make_atm)
            
            # MAINTAINING HYDROSTATIC BALANCE
            var.y = np.vstack(atm.n_0)*var.ymix
            
            # calculating the changing of y
            var = self.f_dy(var, para)
            
            # save values of the current step
            var, para = self.save_step(var, para)
            
            # adjusting the step-size
            var = self.odesolver.step_size(var, para)
            
            if use_print_prog == True and para.count % vulcan_cfg.print_prog_num==0:
                self.output.print_prog(var)
            if use_live_plot == True and para.count % vulcan_cfg.live_plot_frq ==0:
                self.output.plot_update(var, atm, para)
                          
        
    def backup(self, var):
        var.y_prev = np.copy(var.y)
        var.dy_prev = np.copy(var.dy)
        var.atom_loss_prev = var.atom_loss.copy()
        return var
        
    def update_mu_dz(self, var, data_atm, make_atm): #y, ni, spec, Tco, pco
        
        # gravity
        g = make_atm.g
        
        # calculating mu (mean molecular weight)
        data_atm = make_atm.mean_mass(var, data_atm, ni)
        data_atm.Hp = kb*data_atm.Tco/(data_atm.mu/Navo*g)
        
        # calculating the pressure at interface
        data_atm = make_atm.f_pico(data_atm)
        
        # updating dz, zco, and dzi
        dz, dzi, pico = data_atm.dz, data_atm.dzi, data_atm.pico
        dz = data_atm.Hp* np.log(pico[:-1]/np.roll(pico,-1)[:-1])
        data_atm.zco = np.insert(np.cumsum(dz),0, 0.) # cumulative sum of dz
        
        # for the j grid, dzi[j] counts the distance from the grid above and dz[j-1] from the grid below
        dzi = 0.5*(dz + np.roll(dz,1))
        dzi = dzi[1:]
        
        data_atm.dz, data_atm.dzi, data_atm.pico = dz, dzi, pico
        
        return data_atm
    

    # function calculating the change of y
    def f_dy(self, var, para): # y, y_prev, ymix, yconv, count, dt
        if para.count == 0: 
            var.dy, var.dydt = 1., 1.
            return var
        y, ymix, y_prev = var.y, var.ymix, var.y_prev    
        dy =  np.abs(y - y_prev)
        dy[ymix < vulcan_cfg.mtol] = 0   
        dy[y < vulcan_cfg.atol] = 0 
        dy = np.amax( dy[y>0]/y[y>0] )
        
        var.dy, var.dydt = dy, dy/var.dt
        
        return var
    
    
    def conv(self, var, para, out=False, print_freq=200):
        '''
        funtion returns TRUE if the convergence condition is satisfied
        '''
        st_factor, mtol_conv, atol, yconv_cri, slope_cri = vulcan_cfg.st_factor, vulcan_cfg.mtol_conv, vulcan_cfg.atol, vulcan_cfg.yconv_cri, vulcan_cfg.slope_cri
        y, ymix, ymix_time, t_time = var.y.copy(), var.ymix.copy(), var.ymix_time, var.t_time
        count = para.count
        
        # if t < trun_min: indx = -100
        indx = np.abs(t_time-var.t*st_factor).argmin()   
     
        if indx == para.count-1: indx-=1  #Important!! For dt larger than half of the runtime (count-1 is the last one) 
        longdy = np.abs(ymix_time[count-1] - ymix_time[indx])
        longdy[ymix < mtol_conv] = 0
        longdy[y < atol] = 0 
        indx_max = np.nanargmax(longdy/ymix)
        longdy = np.amax( longdy[ymix>0]/ymix[ymix>0] )
        longdydt = longdy/(t_time[-1]-t_time[indx])
        # store longdy and longdydt
        var.longdy, var.longdydt = longdy, longdydt

        if longdy < yconv_cri and longdydt < slope_cri: 
            return True
        return False
    
    def stop(self, var, para):
        '''
        This function is 
        '''
        if var.t > vulcan_cfg.trun_min and para.count > vulcan_cfg.count_min and self.conv(var, para):
            print ('Integration successful with ' + str(para.count) + ' steps and long dy, long dydt = ' + str(var.longdy) + ' ,' + str(var.longdydt) + '\nContinue with quasi-steady runs...') 
            para.end_case = 1
            return True
        elif var.t > vulcan_cfg.runtime:
            print ('Integration not completed...\nMaximal allowed runtime exceeded ('+ \
            str (vulcan_cfg.runtime) + ' sec)!')
            para.end_case = 2
            return True
        elif para.count > vulcan_cfg.count_max:
            print ('Integration not completed...\nMaximal allowed steps exceeded (' + \
            str (vulcan_cfg.count_max) + ')!')
            para.end_case = 3
            return True
    
    def save_step(self, var, para):
        '''
        save current values of y and add 1 to the counter
        '''
        tmp = list(var.y)
        # y_time is initially []   
        var.y_time.append(tmp)
        var.ymix_time.append(var.ymix)
        var.t_time.append(var.t)
                
        # var.dt_time.append(var.dt)
        # the dt_time can be post-calculated
        
        # only used in PI_control
        var.dy_time.append(var.y)
        var.dydt_time.append(var.dydt)
        var.atom_loss_time.append(var.atom_loss.values() )
        
        var.t += var.dt   
        para.count += 1

        return var, para
        
        
class QuasiSteady(Integration):
    
    def __init__(self, ros2, output):
        Integration.__init__(self, ros2, output)

    def __call__(self, var, atm, para, make_atm):
        #Integration.__call__(self, var, atm, para, make_atm)
        
        # the cases for t > runtime or count > count_max
        if para.end_case > 1: return
        
        use_print_prog, use_live_plot = vulcan_cfg.use_print_prog, vulcan_cfg.use_live_plot
        
        while not self.stop(var, para): 
            
            var = self.backup(var)
            var, para = self.odesolver.one_step(var, atm, para)
           
            if para.count % vulcan_cfg.update_frq == 0: 
                atm = self.update_mu_dz(var, atm, make_atm)
            
            # MAINTAINING HYDROSTATIC BALANCE
            var.y = np.vstack(atm.n_0)*var.ymix
            
            # calculating the changing of y
            var = self.f_dy(var, para)
            
            # save values of the current step
            var, para = self.save_step(var, para)
            
            # adjusting the step-size
            var = self.final_step_size(var, para)
            
            if use_print_prog == True and para.count % vulcan_cfg.print_prog_num==0:
                self.output.print_prog(var)
            if use_live_plot == True and para.count % vulcan_cfg.live_plot_frq ==0:
                self.output.plot_update(var, atm, para)
        
               
    # override stop in the parent class Integration
    def stop(self, var, para):
        
        if var.dt < vulcan_cfg.dt_std:
            print ('After the quasi-steady runs the integration has finished with ' + str(para.count) + ' steps and long dy, long dydt = ' + str(var.longdy) + ' ,' + str(var.longdydt) + ' in total.') 
            self.output.print_end_msg(var, para)
            return True
            
    
    # override step_size in the parent class Integration
    def final_step_size(self, var, para):
        '''
        step-size control that decreases the stepsize
        '''
        var.dt *= vulcan_cfg.dt_var_min 
        return var
        

class ODESolver(object):
    
    def __init__(self): # do I always need to update var, atm, para ?
        
        self.mtol = vulcan_cfg.mtol
        self.atol = vulcan_cfg.atol       
        
    def diffdf(self, var, atm):  # input y,dzi,Kzz
        """
        function of eddy diffusion with zero-flux boundary conditions and non-uniform grids (dzi)
        in the form of Aj*y_j + Bj+1*y_j+1 + Cj-1*y_j-1
        """
        
        y = var.y.copy()
        ysum = np.sum(y, axis=1)
        dzi = atm.dzi.copy()
        Kzz = atm.Kzz.copy()
           
        A, B, C = np.zeros(nz), np.zeros(nz), np.zeros(nz)

        A[0] = -1./(dzi[0])*(Kzz[0]/dzi[0]) *(ysum[1]+ysum[0])/2. /ysum[0]     
        B[0] = 1./(dzi[0])*(Kzz[0]/dzi[0]) *(ysum[1]+ysum[0])/2. /ysum[1] 
        C[0] = 0 
        A[nz-1] = -1./(dzi[nz-2])*(Kzz[nz-2]/dzi[nz-2]) *(ysum[nz-1]+ysum[nz-2])/2. /ysum[nz-1] 
        B[nz-1] = 0 
        C[nz-1] = 1./(dzi[nz-2])*(Kzz[nz-2]/dzi[nz-2]) *(ysum[nz-1]+ysum[nz-2])/2. /ysum[nz-2] 
       
        for j in range(1,nz-1):  
            A[j] = -2./(dzi[j-1] + dzi[j])* ( Kzz[j]/dzi[j]*(ysum[j+1]+ysum[j])/2. + Kzz[j-1]/dzi[j-1]*(ysum[j]+ysum[j-1])/2. ) /ysum[j]  
            B[j] = 2./(dzi[j-1] + dzi[j])*Kzz[j]/dzi[j] *(ysum[j+1]+ysum[j])/2. /ysum[j+1]
            C[j] = 2./(dzi[j-1] + dzi[j])*Kzz[j-1]/dzi[j-1] *(ysum[j]+ysum[j-1])/2. /ysum[j-1]
    
        tmp0 = A[0]*y[0] + B[0]*y[1]
        tmp1 = np.ndarray.flatten( (np.vstack(A[1:nz-1])*y[1:(nz-1)] + np.vstack(B[1:nz-1])*y[1+1:(nz-1)+1] + np.vstack(C[1:nz-1])*y[1-1:(nz-1)-1]) )
        tmp2 = (A[nz-1]*y[nz-1] +C[nz-1]*y[nz-2]) 
        diff = np.append(np.append(tmp0, tmp1), tmp2)
        diff = diff.reshape(nz,ni)

        return diff
            
    def jac_tot(self, var, atm): # input y,dzi,Kzz
        """
        jacobian matrix for dn/dt + dphi/dz = P - L (including diffusion)
        zero-flux BC:  1st derivitive of y is zero
        """
        
        y = var.y.copy()
        ysum = np.sum(y, axis=1)
        dzi = atm.dzi.copy()
        Kzz = atm.Kzz.copy()
        
        dfdy = achemjac(y, atm.M, var.k)
        j_indx = []
        
        for j in range(nz):
            j_indx.append( np.arange(j*ni,j*ni+ni) )

        for j in range(1,nz-1): 
            # excluding the buttom and the top cell
            # at j level consists of ni species 
            dfdy[j_indx[j], j_indx[j]] +=  -2./(dzi[j-1] + dzi[j])*( Kzz[j]/dzi[j]*(ysum[j+1]+ysum[j])/2. + Kzz[j-1]/dzi[j-1]*(ysum[j-1]+ysum[j])/2. ) /ysum[j] 
            dfdy[j_indx[j], j_indx[j+1]] += 2./(dzi[j-1] + dzi[j])*( Kzz[j]/dzi[j]*(ysum[j+1]+ysum[j])/(2.*ysum[j+1]) )
            dfdy[j_indx[j], j_indx[j-1]] += 2./(dzi[j-1] + dzi[j])*( Kzz[j-1]/dzi[j-1]*(ysum[j-1]+ysum[j])/(2.*ysum[j-1]) )
    
        dfdy[j_indx[0], j_indx[0]] += -1./(dzi[0])*(Kzz[0]/dzi[0]) * (ysum[1]+ysum[0])/(2*ysum[0])
        dfdy[j_indx[0], j_indx[1]] += 1./(dzi[0])*(Kzz[0]/dzi[0]) * (ysum[1]+ysum[0])/(2*ysum[1])   
        dfdy[j_indx[nz-1], j_indx[nz-1]] += -1./(dzi[nz-2])*(Kzz[nz-2]/dzi[nz-2]) *(ysum[(nz-1)-1]+ysum[nz-1])/(2.*ysum[nz-1])   
        dfdy[j_indx[nz-1], j_indx[(nz-1)-1]] += 1./(dzi[nz-2])*(Kzz[nz-2]/dzi[nz-2])* (ysum[(nz-1)-1]+ysum[nz-1])/(2.*ysum[(nz-1)-1])  

        return dfdy
    
          
    def clip(self, var, para, pos_cut = vulcan_cfg.pos_cut, nega_cut = vulcan_cfg.nega_cut):
        '''
        function to clip samll and negative values
        and to calculate the particle loss
        '''
        y = var.y
         
        para.small_y += np.abs(np.sum(y[np.logical_and(y<pos_cut, y>=0)]))
        para.nega_y += np.abs(np.sum(y[np.logical_and(y>nega_cut, y<=0)]))
        y[np.logical_and(y<pos_cut, y>=nega_cut)] = 0.
        
        var = self.loss(var)
        
        # store y
        var.y = y
        
        return var , para
        
    def loss(self, data_var): 
        
        y = data_var.y
        atom_list = vulcan_cfg.atom_list
        
        # changed atom_tot to dictionary atom_sum
        atom_sum = data_var.atom_sum
        
        for atom in atom_list:
            data_var.atom_sum[atom] = np.sum([compo[compo_row.index(species[i])][atom] * data_var.y[:,i] for i in range(ni)])
            data_var.atom_loss[atom] = (data_var.atom_sum[atom] - data_var.atom_ini[atom])/data_var.atom_ini[atom]

        return data_var
        
    def step_ok(self, var, para, loss_eps = vulcan_cfg.loss_eps, rtol = vulcan_cfg.rtol):
        if np.all(var.y>=0) and np.amax( np.abs( np.array(var.atom_loss.values()) - np.array(var.atom_loss_prev.values()) ) )<loss_eps and para.delta<=rtol:
            return True
        else:
            return False
            
    def step_reject(self, var, para, loss_eps = vulcan_cfg.loss_eps, rtol = vulcan_cfg.rtol):
        
        if para.delta > rtol: # truncation error larger than the tolerence value
            para.delta_count += 1
            
        elif np.any(var.y < 0):             
            para.nega_count += 1
            if vulcan_cfg.use_print_prog == True:
                self.print_nega(var,para) # print the info for the negative solutions (where y < 0)
            # print input: y, t, count, dt
            

        else: # meaning np.amax( np.abs( np.abs(y_loss) - np.abs(loss_prev) ) )<loss_eps
            para.loss_count +=1
            if vulcan_cfg.use_print_prog == True:
                self.print_lossBig(para)
        
        
        var = self.reset_y(var) # reset y and dt to the values at previous step
        
        if var.dt < vulcan_cfg.dt_min:
            var.dt = vulcan_cfg.dt_min
            var.y[var.y<0] = 0. # clipping of negative values
            if vulcan_cfg.use_print_prog == True:
                print ('Keep producing negative values! Clipping negative solutions and moving on!')
            return True
        
        return False
            
    def reset_y(self, var, dt_reduc = vulcan_cfg.dt_var_min):
        '''
        reset y and reduce dt by dt_reduc
        '''
        
        # reset and store y and dt
        var.y = var.y_prev
        var.dt *= dt_reduc
        # var.dt = np.maximum(var.dt, vulcan_cfg.dt_min)

        return var
        
    def print_nega(self, data_var, data_para): 
        
        nega_i = np.where(data_var.y<0)
        print ('Negative y at time ' + str("{:.2e}".format(data_var.t)) + ' and step: ' + str(data_para.count) )
        print ('Negative values:' + str(data_var.y[data_var.y<0]) )
        print ('from levels: ' + str(nega_i[0]) )
        print ('species: ' + str([species[_] for _ in nega_i[1]]) )
        print ('dt= ' + str(data_var.dt))
        print ('...reset dt to dt*0.2...')
        print ('------------------------------------------------------------------')
    
    def print_lossBig(self, para):
        
        print ('partical conservation is violated too large (by numerical errors)')
        print ('at step: ' + str(para.count))
        print ('------------------------------------------------------------------')
        
    def thomas_vec(a, b, c, d): 
        '''
        Thomas vectorized solver, a b c d refer to http://en.wikipedia.org/wiki/Tridiagonal_matrix_algorithm
        d is a matrix
        not used in this current version
        '''
        # number of equations
        nf = len(a) 
        aa, bb, cc, dd = map(np.copy, (a, b, c, d))  
        # d needs to reshape
        dd = dd.reshape(nf,-1)
        #C' and D'
        cp = [cc[0]/bb[0]]; dp = [dd[0]/bb[0]]  
        x = np.zeros((nf, np.shape(dd)[1]))
  
        for i in range(1, nf-1):
            cp.append( cc[i]/(bb[i] - aa[i]*cp[i-1]) ) 
            dp.append( (dd[i] - aa[i]*dp[i-1])/(bb[i] - aa[i]*cp[i-1]) )  
   
        dp.append( (dd[(nf-1)] - aa[(nf-1)]*dp[(nf-1)-1])/(bb[(nf-1)] - aa[(nf-1)]*cp[(nf-1)-1]) ) # nf-1 is the last element
        x[nf-1] = dp[nf-1]/1
        for i in range(nf-2, -1, -1):
            x[i] = dp[i] - cp[i]*x[i+1]
        
        return x


class Ros2(ODESolver):
    '''
    class inheritance from ODEsolver for 2nd order Rosenbrock solver 
    '''
    def __init__(self):
        ODESolver.__init__(self)
         
    def store_bandM(self, a, nb, nn):
        """
        store block-tridiagonal matrix(bandwidth=1) into diagonal ordered form 
        (http://docs.scipy.org/doc/scipy/reference/generated/scipy.linalg.solve_banded.html) 
        a : square block-tridiagonal matirx
        nb: size of the block matrix (number of species)
        nn: number of the block matrices (number of layers)
        """
    
        # band width (treat block-banded as banded matrix)
        bw = 2*nb-1 
        ab = np.zeros((2*bw+1,nb*nn))

        # first 2 columns
        for i in range(0,2*nb):
            ab[-(2*nb+i):,i] = a[0:2*nb+i,i]
    
        # middle
        for i in range(2*nb, nn*nb-2*nb):
            ab[:,i] = a[(i-2*nb+1):(i-2*nb+1)+(2*bw+1),i] 
    
        # last 2 columns
        for ne,i in enumerate(range(nn*nb-2*nb,nn*nb)):
            ab[:(2*bw+1 -ne),i] = a[-(2*bw+1 -ne):,i]

        return (ab, bw)
    
       
    def solver(self, var, atm, para):
        """
        2nd order Rosenbrock [Verwer et al. 1997] with banded-matrix solver
        """
        y, ymix, h, k = var.y, var.ymix, var.dt, var.k
        M, dzi, Kzz = atm.M, atm.dzi, atm.Kzz

        diffdf = self.diffdf
        jac_tot = self.jac_tot
    
        r = 1. + 1./2.**0.5
    
        df = chemdf(y,M,k).flatten() + diffdf(var, atm).flatten()
        dfdy = jac_tot(var, atm)
  
        lhs = 1./(r*h)*np.identity(ni*nz) - dfdy
        lhs_b, bw = self.store_bandM(lhs,ni,nz)
        k1_flat = scipy.linalg.solve_banded((bw,bw),lhs_b,df)
        k1 = k1_flat.reshape(y.shape)

        yk2 = y + k1/r
        df = chemdf(yk2,M,k).flatten() + diffdf(var, atm).flatten()

        rhs = df - 2./(r*h)*k1_flat
        k2 = scipy.linalg.solve_banded((bw,bw),lhs_b,rhs)
        k2 = k2.reshape(y.shape)

        sol = y + 3./(2.*r)*k1 + 1/(2.*r)*k2  

        delta = np.abs(sol-yk2)
        delta[ymix < self.mtol] = 0
        delta[sol < self.atol] = 0
        delta = np.amax( delta[sol>0]/sol[sol>0] )
    
        var.y = sol
        var.ymix = var.y/np.vstack(np.sum(var.y,axis=1))
        para.delta = delta
    
        return var, para  
        
    
    def one_step(self, var, atm, para):
        
        while True:
           var, para = self.solver(var, atm, para)
           
           # clipping small negative values and also calculating atomic loss (atom_loss)  
           var , para = self.clip(var, para) 
            
           if self.step_ok(var, para): break
           elif self.step_reject(var, para): break # giving up and moving on
               
        return var, para                    
        
    def step_size(self, var, para, dt_var_min = vulcan_cfg.dt_var_min, dt_var_max = vulcan_cfg.dt_var_max, dt_min = vulcan_cfg.dt_min, dt_max = vulcan_cfg.dt_max):  
        """
        step-size control by delta(truncation error) for the Rosenbrock method
        """
        y = var.y
        h = var.dt
        delta = para.delta
        rtol = vulcan_cfg.rtol
               
        if delta==0: delta = 0.01*rtol
        h_factor = 0.9*(rtol/delta)**0.5
        h_factor = np.maximum(h_factor, dt_var_min)    
        h_factor = np.minimum(h_factor, dt_var_max)    
        
        h *= h_factor
        h = np.maximum(h, dt_min)
        h = np.minimum(h, dt_max)
        
        # store the adopted dt
        var.dt = h
        
        return var
            
    
class Output(object):
    
    def __init__(self):
        
        output_dir, out_name, plot_dir = vulcan_cfg.output_dir, vulcan_cfg.out_name, vulcan_cfg.plot_dir

        if not os.path.exists(output_dir): os.makedirs(output_dir)
        if not os.path.exists(plot_dir): os.makedirs(plot_dir)
        
        if os.path.isfile(output_dir+out_name):
            # Fix Python 3.x and 2.x.
            try: input = raw_input
            except NameError: pass
            input("  The output file: " + str(out_name) + " already exists.\n"
                      "  Press enter to overwrite the existing file,\n"
                      "  or Ctrl-Z and Return to leave and choose a different out_name in vulcan_cfg.")
        
    def print_prog(self, var):
        print ('time: ' +str("{:.2e}".format(var.t)) + '  and dt= ' + str(var.dt) )
        print ('longdy=' + str(var.longdy) + ' and longdy/dt= ' + str(var.longdydt) )
        print ('------------------------------------------------------' )

    def print_end_msg(self, var, para ): 
        print ("After ------- %s seconds -------" % ( time.time()- para.start_time ) + ' s CPU time') 
        print ('VULCAN has sucesfully run to steady-state with ' + str(para.count) + ' steps and ' + str("{:.2e}".format(var.t)) + ' s' )
        print ('long dy = ' + str(var.longdy) + ' and long dy/dy = ' + str(var.longdydt) )
        
        print ('total atom loss:')
        for atom in vulcan_cfg.atom_list: print (atom + ': ' + str(var.atom_loss[atom]) + ' ')
      
        print ('negative solution counter:')
        print (para.nega_count)
        print ('loss rejected counter:')
        print (para.loss_count)
        print ('delta rejected counter:')
        print (para.delta_count)
        
        print ('------ Live long and prosper \V/ ------') 
    
    def save_out(self, var, atm, para, dname): 
        
        output_dir, out_name = vulcan_cfg.output_dir, vulcan_cfg.out_name
        fq = vulcan_cfg.out_y_time_freq
        
        # convert lists into numpy arrays
        # and slicing time-sequential data to reduce ouput filesize
        for key in ['dy_time', 'dydt_time', 'ymix_time',  'y_time', 't_time', 'dt_time', 'atom_loss_time' ]:
            as_nparray = np.array(getattr(var, key))
            setattr(var, key, as_nparray)
            # sort of like: var."key" = np.array(var."key")
        
        # plotting
        if vulcan_cfg.use_plot_evo == True: 
            self.plot_evo(var)
        if vulcan_cfg.use_plot_end == True:
            self.plot_end(var, atm, para)
        
        for key in ['dy_time', 'dydt_time', 'ymix_time',  'y_time', 't_time', 'dt_time', 'atom_loss_time' ]:
            np_slicing = getattr(var, key)[::fq]
            setattr(var, key, np_slicing)
            
        var_save = vars(var).copy()
        # include the numerical order of species
        var_save['species'] = species
        
        # delete the lambda funtions
        del var_save['k_fun']; del var_save['kinf_fun']; del var_save['k_fun_new']
        
        output_file = dname + '/' + output_dir + out_name
        with open(output_file, 'wb') as outfile:
            if vulcan_cfg.output_humanread == True: # human-readable form, less efficient 
                outfile.write(str({'variable': var_save, 'atm': vars(atm), 'parameter': vars(para)}))
            else:
                pickle.dump( {'variable': var_save, 'atm': vars(atm), 'parameter': vars(para) }, outfile, protocol=pickle.HIGHEST_PROTOCOL)
                # how to add  'config': vars(vulcan_cfg) ?
        
            
    def plot_update(self, var, atm, para):
        
        images = []
        colors = ['b','g','r','c','m','y','k','orange','pink', 'grey',\
        'darkred','darkblue','salmon','chocolate','mediumspringgreen','steelblue','plum','hotpink']
        plt.figure('live mixing ratios')
        plt.ion()
        color_index = 0
        for sp in vulcan_cfg.live_plot_spec:
            line, = plt.plot(var.ymix[:,species.index(sp)], atm.pco/1.e6, color = colors[color_index], label=sp)
            color_index +=1
            images.append((line,))
        
        plt.title(str(para.count)+' steps and ' + str("{:.2e}".format(var.t)) + ' s' )
        plt.gca().set_xscale('log')       
        plt.gca().set_yscale('log') 
        plt.gca().invert_yaxis() 
        plt.xlim(1.E-20, 1.)
        plt.ylim((1.E3,1.E-4))
        plt.legend(frameon=0, prop={'size':14}, loc=3)
        plt.xlabel("Mixing Ratios")
        plt.ylabel("Pressure (bar)")
        plt.show(block=0)
        plt.pause(0.001)
        plt.clf()
    
    def plot_end(self, var, atm, para):
        
        plot_dir = vulcan_cfg.plot_dir
        colors = ['b','g','r','c','m','y','k','orange','pink', 'grey',\
        'darkred','darkblue','salmon','chocolate','mediumspringgreen','steelblue','plum','hotpink']
        
        plt.figure('live mixing ratios')
        color_index = 0
        for sp in vulcan_cfg.live_plot_spec:
            line, = plt.plot(var.ymix[:,species.index(sp)], atm.pco/1.e6, color = colors[color_index], label=sp)
            color_index +=1
                  
        plt.title(str(para.count)+' steps and ' + str("{:.2e}".format(var.t)) + ' s' )
        plt.gca().set_xscale('log')       
        plt.gca().set_yscale('log') 
        plt.gca().invert_yaxis() 
        plt.xlim(1.E-20, 1.)
        plt.ylim((1.E3,1.E-4))
        plt.legend(frameon=0, prop={'size':14}, loc=3)
        plt.xlabel("Mixing Ratios")
        plt.ylabel("Pressure (bar)")
        plt.savefig(plot_dir + 'mix.png')       
        if vulcan_cfg.use_live_plot == True:
            # plotting in the same window of real-time plotting
            plt.draw()
        elif vulcan_cfg.use_PIL == True: # plotting in a new window with PIL package            
            plot = Image.open(plot_dir + 'mix.png')
            plot.show()
            plt.close()
            
    def plot_evo(self, var, plot_j=-1, dn=1):
        
        plot_spec = vulcan_cfg.plot_spec
        plot_dir = vulcan_cfg.plot_dir
        plt.figure('evolution')
    
        for i,sp in enumerate(vulcan_cfg.plot_spec):
            plt.plot(var.t_time[::dn], var.ymix_time[::dn,plot_j,species.index(sp)],c = plt.cm.rainbow(float(i)/len(plot_spec)),label=sp)

        plt.gca().set_xscale('log')       
        plt.gca().set_yscale('log') 
        plt.xlabel('time')
        plt.ylabel('mixing ratios')
        plt.ylim((1.E-30,1.))
        plt.legend(frameon=0, prop={'size':14}, loc='best')
        plt.savefig(plot_dir + 'evo.png')
        if vulcan_cfg.use_PIL == True:
            plot = Image.open(plot_dir + 'evo.png')
            plot.show()
            plt.close()
        # else: plt.show(block = False)
        
    def plot_TP(self, atm):
        
        plot_dir = vulcan_cfg.plot_dir
        
        plt.figure('TP')
        plt.semilogy( atm.Tco, atm.pco/1.e6, c='black')
        plt.gca().invert_yaxis()
        plt.xlabel("Temperature (K)")
        plt.ylabel("Pressure (bar)")
        plot_name = plot_dir + 'TP.png'
        plt.savefig(plot_name)
        if vulcan_cfg.use_PIL == True:        
            plot = Image.open(plot_name)
            plot.show()
            # close the matplotlib window
            plt.close()
        else: plt.show(block = False)
            
class SemiEU(ODESolver):
    '''
    class inheritance from ODEsolver for semi-implicit Euler solver 
    '''
    def __init__(self):
        ODESolver.__init__(self)
           
    def solver(self, var, atm):
        """
        semi-implicit Euler solver (1st order) 
        """
        y, ymix, h, k = var.y, var.ymix, var.dt, var.k
        M, dzi, Kzz = atm.M, atm.dzi, atm.Kzz

        diffdf = self.diffdf
        jac_tot = self.jac_tot
        
        df = chemdf(y,M,k).flatten() + diffdf(var, atm).flatten()
        dfdy = jac_tot(var, atm)        
        aa = np.identity(ni*nz) - h*dfdy
        aa = scipy.linalg.solve(aa,df)
        aa = aa.reshape(y.shape)
        y = y + aa*h
        
        var.y = y
        var.ymix = var.y/np.vstack(np.sum(var.y,axis=1))
        
        return var

    def step_ok(self, var, para, loss_eps = vulcan_cfg.loss_eps):
        if np.all(var.y>=0) and np.amax( np.abs( np.array(var.atom_loss.values()) - np.array(var.atom_loss_prev.values()) ) )<loss_eps:
            return True
        else:
            return False    
    
    def one_step(self, var, atm, para):
        
        while True:
           var = self.solver(var, atm)
           
           # clipping small negative values and also calculating atomic loss (atom_loss)  
           var , para = self.clip(var, para) 
            
           if self.step_ok(var, para): break
           elif self.step_reject(var, para): break # giving up and moving on
               
        return var, para                    
    
    def step_size(self, var, para):
        '''
        PID control required for all semi-Euler like methods
        '''
        dt_var_min, dt_var_max, dt_min, dt_max = vulcan_cfg.dt_var_min, vulcan_cfg.dt_var_max, vulcan_cfg.dt_min, vulcan_cfg.dt_max
        PItol = vulcan_cfg.PItol
        dy, dy_prev, h = var.dy, var.dy_prev, var.dt
        
        if dy == 0 or dy_prev == 0: 
            var.dt = np.minimum(h*2.,dt_max)
            return var
            
        if para.count > 0:
            
            h_factor = (dy_prev/dy)**0.075 * (PItol/dy)**0.175
            h_factor = np.maximum(h_factor, dt_var_min)    
            h_factor = np.minimum(h_factor, dt_var_max)
            h *= h_factor
            h = np.maximum(h, dt_min)
            h = np.minimum(h, dt_max)
        
        # store the adopted dt
        var.dt = h

        return var


class SparSemiEU(SemiEU):
    '''
    class inheritance from SemiEU.
    It is the same semi-implicit Euler solver except for utilizing sparse-matrix solvers
    '''
    def __init__(self):
        SemiEU.__init__(self)
    
    # override solver
    def solver(self, var, atm):
        """
        sparse-matrix semi-implicit Euler solver (1st order) 
        """
        y, ymix, h, k = var.y, var.ymix, var.dt, var.k
        M, dzi, Kzz = atm.M, atm.dzi, atm.Kzz

        diffdf = self.diffdf
        jac_tot = self.jac_tot
        
        df = chemdf(y,M,k).flatten() + diffdf(var, atm).flatten()
        dfdy = jac_tot(var, atm)        
                
        aa = sparse.csc_matrix( np.identity(ni*nz) - h*dfdy )
        aa = sparse.linalg.spsolve(aa,df)
        aa = aa.reshape(y.shape)
        y = y + aa*h
        
        var.y = y
        var.ymix = var.y/np.vstack(np.sum(var.y,axis=1))

        return var 
    

### back-up methods: extrapolation semi_implicit Euler ###