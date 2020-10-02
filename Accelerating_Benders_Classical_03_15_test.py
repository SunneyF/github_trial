# -*- coding: utf-8 -*-
"""
Created on Wed Jan 29 09:39:16 2020

@author: sunney
"""
# Pseudocode for CCB algorithm 
# source:G. K. D. Saharidis et al. / Intl. Trans. in Op. Res. 17 (2010) 221–237
#tic = time.time()        
def switchLabel(x):
    return{
                1:'LT',
                2:'ST',
                3:'T',
                4:'LM',
                5:'SM',
                6:'MTC',
                7:'D'
             }[x]

def unique(list1): 
      
    # insert the list to the set 
    list_set = set(list1) 
    # convert the set to the list 
    unique_list = (list(list_set)) 
    return unique_list


import gurobipy as grb 
import numpy as np
import csv
T=3
A=480
# Extracting data from txt files ############################################

# Inputs required  ---------------------------------------------------------- 
cost_vector = [1,2,3,4] # qualification cost levels
# -------- Machines ---------------------------------------------------------

# 1. KO (ComponentArea) # dictionary
# 2. Work_Center_ID (code for the work center)
# 3. ComponentAreaToMachine till exempel: {1:[1,3,5,10], 2: [97, 74,...],....}
     # key is the code of the component area and values contain a vector of
     # all the machines in respective component area
     
# 4. Categories = [1,2]  refers to the machine category
# 5. CategoriesToMachine = {1:[1,3,5,10], 2: [97, 74,...],....}
# 6. Number of equipments per work center = {1:1,...., 5:21, ...}
     
f= open("MachineData.txt","r")
f1 = f.readlines()
KO ={}
Work_Center_ID = []
ComponentAreaToMachine = {}
Categories = []
CategoriesToMachine = {}
MachineToCatTuple = grb.tuplelist()
NumberofEquipments = {}
NumberofComponentAreas = 10

with open('MachineData.txt') as csv_file:
    csv_reader = csv.reader(csv_file, delimiter=',')
    line_count = 0
    for row in csv_reader:
#       print(row[3])
        Work_Center_ID.append(int(row[0]))
        ComponentAreaToMachine[int(row[0])] = int(row[3])
        CategoriesToMachine[int(row[0])] = int(row[1])
        MachineToCatTuple+=[(int(row[0]),int(row[1]))] 
        Categories.append(int(row[1]))
        NumberofEquipments[int(row[0])] = int(row[2])
Categories = unique(Categories) 
ComponentAreas = range(1,NumberofComponentAreas+1)
  
# --------Part type--------------------------------------------------------
# 6. PartType = [1:...]
#    JobTypeData = {'1': [100,400,500,600], '2':[100,200,300], '3':[100,200],...}     
#    Pijk = grb.tuplelist([(1,100,work_center_id, 10.0),..,..]}
PartType = []

#dictLists = dict((key, []) for key in ["xcg", "bsd", ...])

with open('ProcessingTime.txt') as csv_file:
    csv_reader = csv.reader(csv_file, delimiter=',')
    line_count = 0
    for row in csv_reader:
      PartType.append(int(row[0]))
csv_file.close()      
PartType = unique(PartType)      
JobTypeData = dict((key, []) for key in PartType) 
Pijk = grb.tuplelist()  

with open('ProcessingTime.txt') as csv_file:
    csv_reader = csv.reader(csv_file, delimiter=',')
    line_count = 0
    for row in csv_reader:
     
      JobTypeData[int(row[0])].append(int(row[1]))
      Pijk+= [(int(row[0]),int(row[1]),int(row[2]),float(row[3]))]
csv_file.close()

for key in JobTypeData.keys():
    JobTypeData[key]= unique(JobTypeData[key])
#-----------Qualification cost --------------------------------------------

# set the qualification cost to zero for machines which have been used before
     
Omega = grb.tuplelist() 
QualCost={}
JobsToCategories = dict(((i,j), []) for i in PartType for j in JobTypeData[i] )

for i in PartType:
    for j in JobTypeData[i]:
        for k in Pijk.select(i,j,'*','*'):
            
            JobsToCategories[i,j].append(CategoriesToMachine[k[2]])
            Omega+=[(i,j,k[2])]  # qualified machines(i,j,k)
            QualCost[i,j,k[2]] =0 
for key,value in JobsToCategories.items():
    
        JobsToCategories[key] = unique(value)

# find set of machines for each i,j where Lambda_ijk =1 and Omega_ijk =0
NotQualified = dict(((key1,key2), []) for key1 in PartType for key2 in JobTypeData[key1]) 

for i in PartType:
    for j in JobTypeData[i]:
        K_qualified=[]
        for kk in Omega.select(i,j,'*'):
            #print(kk[2])
            K_qualified.append(kk[2])
        for k in Work_Center_ID:
            d = MachineToCatTuple.select(k,'*')
           
            if (d[0][1] in JobsToCategories[i,j]):
                                
                if not(k in K_qualified):
                    NotQualified[i,j].append(k) 
                    NotQualified[i,j] = unique(NotQualified[i,j])

for i,j in NotQualified.keys():
    for k in NotQualified[i,j]:
        QualCost[i,j,k]= np.random.permutation(cost_vector)[0]

FeasibleMachines = dict(((key1,key2), []) for key1 in PartType for key2 in JobTypeData[key1])

for i in PartType:
    for j in JobTypeData[i]:
        for k in Pijk.select(i,j,'*','*'):
            FeasibleMachines[i,j].append(k[2])
            for k2 in NotQualified[i,j]:
                FeasibleMachines[i,j].append(k2)
for i in PartType:
    for j in JobTypeData[i]:
        FeasibleMachines[i,j] = unique(FeasibleMachines[i,j])
# Demand data
Demand = dict(((key1,key2,key3), []) for key1 in PartType for key2 in JobTypeData[key1] for key3 in range(1,T+1)) 
cnt=0     

# create a dictionary for categories to machines
CatMachList = dict(((key1), []) for key1 in Categories)
for c in Categories:
    for k in MachineToCatTuple.select('*',c):
        CatMachList[c].append(k[0])
   
with open('DemandData.txt') as csv_file:
    csv_reader = csv.reader(csv_file, delimiter=',')
    line_count = 0
    
    for row in csv_reader:
      if ((int(row[0]) not in PartType) or (int(row[1]) not in JobTypeData[int(row[0])]) or (int(row[4]) <= 0)): # for those where processing time too less
          
          cnt+=1
          continue
      if(int(row[3]) <= 3):
          if (int(row[2]) == 2016):
              Demand[int(row[0]),int(row[1]),1].append(int(row[4]))
          elif (int(row[2]) == 2017):
              Demand[int(row[0]),int(row[1]),1+4].append(int(row[4]))
          else:
              Demand[int(row[0]),int(row[1]),1+8].append(int(row[4]))
              
      elif((int(row[3]) >= 4) & (int(row[3]) <= 6)):
          if (int(row[2]) == 2016):
              Demand[int(row[0]),int(row[1]),2].append(int(row[4]))
          elif (int(row[2]) == 2017):
              Demand[int(row[0]),int(row[1]),2+4].append(int(row[4]))
          else:
              Demand[int(row[0]),int(row[1]),2+8].append(int(row[4]))
              
      elif((int(row[3]) >= 7) & (int(row[3]) <= 9)):
          
          if (int(row[2]) == 2016):
              Demand[int(row[0]),int(row[1]),3].append(int(row[4]))
          elif (int(row[2]) == 2017):
              Demand[int(row[0]),int(row[1]),3+4].append(int(row[4]))
          else:
              Demand[int(row[0]),int(row[1]),3+8].append(int(row[4]))
      elif((int(row[3]) >= 10) & (int(row[3]) <= 12)):
          if (int(row[2]) == 2016):
              Demand[int(row[0]),int(row[1]),4].append(int(row[4]))
          elif (int(row[2]) == 2017):
              Demand[int(row[0]),int(row[1]),4+4].append(int(row[4]))
          else:
              Demand[int(row[0]),int(row[1]),4+8].append(int(row[4]))          
                  

        
csv_file.close()

  
for  key in Demand.keys():
    Demand[key] = sum(Demand[key])
  
# Capacity

Capacity = {}  
factor=1 #  % additional time
for t in range(1,T+1):
    Capacity[1,t] = 5*5000*(1/4)
    Capacity[9,t] = 21*5000*(1/4)
    Capacity[116,t]= 3*5000*(1/4)
    Capacity[49,t] = 2*5000*(1/4)

for i in set(set(Work_Center_ID)-{1,9,116,49}):
    for j in range(1,T+1):
        Capacity[i,j] = 5000*(1/4)

#Set the value of Processing time for not qualified machines as the max processing time +4 hours
for i in PartType:
    for j in JobTypeData[i]:
        maxs =0
        for k2 in Pijk.select(i,j,'*','*'):
                if k2[3] > maxs:
                    maxs = k2[3]*factor
        for k in NotQualified[i,j]:
            Pijk+= [(i,j,k,maxs)]

################################ End of Data  ###############################
            
            
################################ Parameters   ##############################
            
Tau=3
Theta=0.5
gamma=10
T=10
threshold = .7
e= ((1-Theta)/(T*(1.0-threshold)))
A=480
#############################################################################


#%%
import time

# your code


def subproblem_t(t,s_var,x_var_t):
    # Input: s:= array of s[i,j,k,t] values and t is the time period
    model = grb.Model("TRA_SubProblem") 
    x_t ={}
    for i in PartType:
        for j in JobTypeData[i]:
            for k in FeasibleMachines[i,j]:
                if Demand[i,j,t] !=0:
                    x_t[i,j,k] = model.addVar(vtype= grb.GRB.CONTINUOUS,lb=0.0)
    n_t = model.addVar(vtype = grb.GRB.CONTINUOUS,lb=0.0)
    model.setObjective(((Theta)/(T*(1.0-threshold)))*n_t,grb.GRB.MINIMIZE)
    Assignments_t = {}
    for i in PartType:
        for j in JobTypeData[i]:
            if Demand[i,j,t]!=0:
                Assignments_t[i,j] = model.addConstr(grb.quicksum(x_t[i,j,k] for k in FeasibleMachines[i,j]) >= max(0,Demand[i,j,t]),name="Assignments[%s]" % (str(i) + "," + str(j)))

    sijkvariable_t = {}
 
    for i in PartType:
        for j in JobTypeData[i]:
            for k in NotQualified[i,j]:
                if Demand[i,j,t]!=0:
                    sijkvariable_t[i,j,k] = model.addConstr(-x_t[i,j,k] >= -Demand[i,j,t]*max(0,s_var[i,j,k,t]),name= "sijkvariable_t[%s]" % (str(i) + "," + str(j) + "," + str(k)))
                    
                    
    MinmaxCategory_t={}

    for k in Work_Center_ID: 
        MinmaxCategory_t[k] = model.addConstr(Capacity[k,t]*n_t - grb.quicksum(Pijk.select(i,j,k,'*')[0][3]*x_t[i,j,k] for i in PartType for j in JobTypeData[i] if (i,j,k) in x_t) >= -threshold*Capacity[k,t], name= "MinmaxCategory_t[%s]" % str(k))
    
    thres = model.addConstr(-n_t >= -1 + threshold, name= "thres")
    
    
    model.Params.InfUnbdInfo = 1
    model.Params.Method=1
    model.Params.DualReductions = 0
    model.Params.PrePasses = 1
    
    
    model.optimize()
    
    pi_t_Dual={}
    y_t_Dual ={}
    mu_t_Dual={}
    d_t_Dual=0
    Stage2_Cost_t=0
    pi_t_Dual_ray ={}
    y_t_Dual_ray ={}
    mu_t_Dual_ray = {}
    d_t_Dual_ray = 0
    status= model.status
    rays =0
    
    if (status==grb.GRB.Status.INFEASIBLE): # https://www.gurobi.com/documentation/9.0/refman/optimization_status_codes.html
           rays=1
           for i in PartType:
             for j in JobTypeData[i]:
                if (i,j) in Assignments_t:
                    pi_t_Dual_ray[i,j] = -Assignments_t[i,j].FarkasDual
                
           for i in PartType:
            for j in JobTypeData[i]:
                for k in FeasibleMachines[i,j]:
                    if (i,j,k) in sijkvariable_t:
                      y_t_Dual_ray[i,j,k]= -sijkvariable_t[i,j,k].FarkasDual
           
           for k in Work_Center_ID:
             if k in MinmaxCategory_t:   
                mu_t_Dual_ray[k]= -MinmaxCategory_t[k].FarkasDual
        
           d_t_Dual_ray = -thres.FarkasDual
#            print('The model is infeasible; computing IIS')
#            model.computeIIS()
#            print('\nThe following constraint(s) cannot be satisfied:')
#            
#            for c in model.getConstrs():
#                if c.IISConstr:
#                    print('%s' % c.constrName)
#           exit(0)          

                             
        
    else:
            
        
        for i in PartType:
            for j in JobTypeData[i]:
                if (i,j) in Assignments_t:
                    pi_t_Dual[i,j] =  Assignments_t[i,j].Pi
        
        
        
        for i in PartType:
            for j in JobTypeData[i]:
                for k in FeasibleMachines[i,j]:
                    if (i,j,k) in sijkvariable_t:
                      y_t_Dual[i,j,k] = sijkvariable_t[i,j,k].Pi
        
                
        for k in Work_Center_ID:
            if k in MinmaxCategory_t:
              mu_t_Dual[k] = MinmaxCategory_t[k].Pi
        
        d_t_Dual = thres.Pi
        Stage2_Cost_t = model.objVal
        for i in PartType:
            for j in JobTypeData[i]:
                for k in FeasibleMachines[i,j]:
                    if (i,j,k) in x_t:
                      x_var_t[i,j,k] = x_t[i,j,k].X 
    # find extreme rays of the dual
    # output scenarios
    # 1. Primal feasible => rays are empty since dual is bounded; update Stage2_Cost
    # 2. Primal infeasible => 1 - Dual is unbounded => find the extreme rays which leads to unboundedness (Farkas Dual)
    #                         2 - Dual is infeasible => we assume primal sub-problem is always bounded and thus dual is feasible
    # rays ==1 means primal sub-problem is infeasible thus, we generate feasibility cut
    return  rays,pi_t_Dual,mu_t_Dual,y_t_Dual,d_t_Dual,pi_t_Dual_ray,mu_t_Dual_ray,y_t_Dual_ray,d_t_Dual_ray,Stage2_Cost_t,x_var_t # \pi^{t}_{ij}, \mu^{t}_{k\in \mathcal{K}}, y^{t}_{ijk} (i,j);k\in K_j, d^{t}_k
##############-----------------Dual Sub-problem End---------------------------------------------------------
# custom optimize() function that uses callback
    
def mycallback(model, where): # https://www.gurobi.com/documentation/9.0/examples/cb_py.html

    if where == grb.GRB.Callback.MIPSOL:
        # General MIP callback
        objbnd = model.cbGet(grb.GRB.Callback.MIPSOL_OBJBND)  # Current best objective bound. MIPSOL_OBJBST
        objbst = model.cbGet(grb.GRB.Callback.MIPSOL_OBJBST)
        if (objbnd > -grb.GRB.INFINITY and model._lowerbounds >= objbst) :
            
            model.terminate()
    return
BestLowerBounds={}
BestUpperBounds = {}
TimeLapse= {}
storeQual = {}
for iter in range(1,2): 
    print('New Iteration-%d' %(iter))
    start_time = time.time()
    Theta=0.9
    nCUT=0
    gamma=10
    T=12
    A=gamma*T
    threshold=0.7
 
    # Start
    model = grb.Model('Master-Problem-TRA')
    s = {}
    z= {}
    o ={}
    
    
    for i in PartType:
        for j in JobTypeData[i]:
            for k in NotQualified[i,j]:
                for t in range(1,T+1):
                    if Demand[i,j,t]!=0:
                      s[i,j,k,t] = model.addVar(vtype=grb.GRB.BINARY,lb=0.0)
                    
                    z[i,j,k,t] = model.addVar(vtype=grb.GRB.BINARY)
                    
    for t in range(1,T+1):
        o[t]= model.addVar(vtype=grb.GRB.CONTINUOUS,lb=0.0)
        
    #model.setObjective(((1-Theta)/A)*grb.quicksum(QualCost[i,j,k]*z[i,j,k,t] 
    #    for i in PartType for j in JobTypeData[i] for k in FeasibleMachines[i,j] for t in range(1,T+1)) + sum(o[t] for t in range(1,T+1))
    #                        - eps*sum(s[i,j,k,t] for i in PartType for j in JobTypeData[i] for k in FeasibleMachines[i,j] for t in range(1,T+1)),grb.GRB.MINIMIZE)
    
    
    model.setObjective(((1-Theta)/A)*grb.quicksum(QualCost[i,j,k]*z[i,j,k,t] 
        for i in PartType for j in JobTypeData[i] for k in NotQualified[i,j] for t in range(1,T+1) ) + sum(o[t] for t in range(1,T+1)),grb.GRB.MINIMIZE)
    
        
    temp={}
    for i in PartType:
        for j in JobTypeData[i]:
            for k in NotQualified[i,j]:
                for t in range(1,T+1):
                    if Demand[i,j,t]!=0:
                        temp[i,j,k,t] = model.addConstr(sum(z[i,j,k,l] for l in range(1,t+1))>= s[i,j,k,t])
    Budget={}                
    for t in range(1,T+1):
        Budget[t] = model.addConstr(sum(z[i,j,k,t] for i in PartType for j in JobTypeData[i] for k in FeasibleMachines[i,j] if k in NotQualified[i,j]) <= gamma)
    
    #for i in PartType:
    #    for j in JobTypeData[i]:
    #        for k in FeasibleMachines[i,j]:
    #            for t in range(1,T+1):
    #                s[i,j,k,t].start=newS[i,j,k,t] # heuristic solution
    #                z[i,j,k,t].start = newZ[i,j,k,t]
    model.Params.PrePasses = 1
    model.Params.MIPFocus = 1
    model.Params.Heuristics = 0.05    # 10% time spent on Heuristics, default is 0.05
    model.Params.Method=2
    
    model.optimize()
    
    s_var={}
    z_var={}
    for i in PartType:
        for j in JobTypeData[i]:
            for k in NotQualified[i,j]:
                for t in range(1,T+1):
                    if Demand[i,j,t]!=0:
                        s_var[i,j,k,t] = s[i,j,k,t].X
                
                    z_var[i,j,k,t] = z[i,j,k,t].X
    
    pi_Dual = {} 
    mu_Dual ={} #
    y_Dual = {}
    d_Dual = {}
    pi_Dual_ray={}
    mu_Dual_ray = {}
    y_Dual_ray={}
    d_Dual_ray={}
    
    Stage2Total=12
    Stage2 = {}
    queue_not_alpha_Covered = grb.tuplelist()
    CCB_counter = 1
    lowerbound={}
    lowerbound[0]=0
    UpperBound = {}
    UpperBound[0]=100
    breaks = 0
    lazyConst={}
    for r in range(1,30):
        print("\nIteration %d \n\n" %(nCUT+1))
        s2 =0
        nCUT=nCUT+1
        infeas=0
        infeasi = {}
        x_var={}
        for t in range(1,T+1):
            print("\nSubproblem- %d \n\n" %(t))
            x_var_t = {}
            if r==1:
                
                 rays,pi_t_Dual,mu_t_Dual,y_t_Dual,d_t_Dual,pi_t_Dual_ray,mu_t_Dual_ray,y_t_Dual_ray,d_t_Dual_ray,Stage2_Cost_t,x_var_t = subproblem_t(t,s_var,x_var_t)
            elif (r-1,t) in Stage2:
                 rays,pi_t_Dual,mu_t_Dual,y_t_Dual,d_t_Dual,pi_t_Dual_ray,mu_t_Dual_ray,y_t_Dual_ray,d_t_Dual_ray,Stage2_Cost_t,x_var_t = subproblem_t(t,s_var,x_var_t)
                
            else:
                continue
            Defn_Cut={}
            
            
            if ((rays==0)):
                if (Stage2_Cost_t==0):
                    continue
                
                for i in PartType:
                    for j in JobTypeData[i]:
                        if (i,j) in pi_t_Dual:
                         pi_Dual[i,j,t,nCUT] = pi_t_Dual[i,j]
                for k in Work_Center_ID:
                   if k in mu_t_Dual: 
                    mu_Dual[k,t,nCUT] = mu_t_Dual[k]
                for i in PartType:
                    for j in JobTypeData[i]:
                        for k in FeasibleMachines[i,j]:
                            if (i,j,k) in y_t_Dual:  
                                y_Dual[i,j,k,t,nCUT] = y_t_Dual[i,j,k]
                d_Dual[t,nCUT] = d_t_Dual
                s2 = s2 +  Stage2_Cost_t
                Stage2[r,t] = Stage2_Cost_t
                for i in PartType:
                    for j in JobTypeData[i]:
                        for k in FeasibleMachines[i,j]:
                            if Demand[i,j,t] !=0:
                              x_var[i,j,k,t] = x_var_t[i,j,k] # current value of x
            else:
                  
                  for i in PartType:
                    for j in JobTypeData[i]:
                        if (i,j) in pi_t_Dual_ray:
                         pi_Dual_ray[i,j,t,nCUT] = pi_t_Dual_ray[i,j]
                  for k in Work_Center_ID:
                    mu_Dual_ray[k,t,nCUT] = mu_t_Dual_ray[k]
                  for i in PartType:
                    for j in JobTypeData[i]:
                       for k in FeasibleMachines[i,j]:
                           if (i,j) in y_t_Dual_ray:
                             y_Dual_ray[i,j,k,t,nCUT] = y_t_Dual_ray[i,j,k]
                  d_Dual_ray[t,nCUT] = d_t_Dual_ray
                  infeas = 1
                  infeasi[r,t] = 1
        
        #####
        
        if infeas ==1:
            Stage2Total= 12
        else:
            Stage2Total=s2
            
        
        f_y = ((1-Theta)/A)*grb.quicksum(QualCost[i,j,k]*z_var[i,j,k,t] for i in PartType for j in JobTypeData[i] for k in FeasibleMachines[i,j] for t in range(1,T+1) if (i,j,k,t) in z_var) 
        
        if (min(UpperBound.values()) > max(lowerbound.values()))  :
                # add cuts to RMP    
                #print('continue')
                for t in range(1,T+1):
                    
                        
                     if (nCUT,t) in Stage2:
                          
                          a1 = sum((-threshold*Capacity[k,t])*mu_Dual[k,t,nCUT] for k in Work_Center_ID if (k,t,nCUT) in mu_Dual)
                          a2 = sum(Demand[i,j,t]*pi_Dual[i,j,t,nCUT] for i in PartType for j in JobTypeData[i] if (i,j,t,nCUT) in pi_Dual)
                          a3= (-1+threshold)*d_Dual[t,nCUT]
                          lazyConst[nCUT,t]=model.addConstr(o[t] >= a1 + a2 + grb.quicksum((-Demand[i,j,t]*y_Dual[i,j,k,t,nCUT]*s[i,j,k,t]) for i in PartType for j in JobTypeData[i] for k in FeasibleMachines[i,j] if (i,j,k,t,nCUT) in y_Dual)+ a3)
                          lazyConst[nCUT,t].lazy=1
                          
                          print("Optimality cut-%d,%d" %(t,nCUT))
            
        else:
              breaks = 1  
              break
        print("\n")
        print("\nRE-Solving Master Problem\n\n")
        
#        for i in PartType:
#            for j in JobTypeData[i]:
#                for k in FeasibleMachines[i,j]:
#                    for t in range(1,T+1):
#                        if (i,j,k,t) in s_var:
#                            s[i,j,k,t].start=  s_var[i,j,k,t] # heuristic solution
#                        if(i,j,k,t) in z_var:
#                         z[i,j,k,t].start = z_var[i,j,k,t] # heuristic solution
#        model.Params.PrePasses = 1
#        model.Params.MIPFocus = 1        # MIP Focus on the bounds
#        model.Params.Heuristics = 0.1    # 10% time spent on Heuristics, default is 0.05
#        #model.Params.Method=2
#        model.Params.MIPGap=.05
#        model.Params.TimeLimit = 150
#        model._lowerbounds = max(lowerbound.values())
        dd = ((1-Theta)/A)*sum(QualCost[i,j,k]*z_var[i,j,k,t] for i in PartType for j in JobTypeData[i] for k in FeasibleMachines[i,j] for t in range(1,T+1) if (i,j,k,t) in z_var)
        sf=  sum(o[t].X for t in range(1,T+1))
        lowerbound[r]= model.objbound
        print("\nLower Bound = %f" %lowerbound[r])
        print("Upper Bound = %f" % (dd+Stage2Total))
        UpperBound[r] = (dd+Stage2Total)
        
        if dd >0:
            storeQual[r] = sum(z_var.values())
        
        if breaks==0:
            for i in PartType:
                for j in JobTypeData[i]:
                    for k in FeasibleMachines[i,j]:
                        for t in range(1,T+1):
                            if (i,j,k,t) in s_var:
                                s[i,j,k,t].start=  s_var[i,j,k,t] # heuristic solution
                            if(i,j,k,t) in z_var:
                                z[i,j,k,t].start = z_var[i,j,k,t] # heuristic solution
            model.Params.PrePasses = 1
            model.Params.MIPFocus = 2        # MIP Focus on the bounds
            model.Params.Heuristics = 0.05    # 10% time spent on Heuristics, default is 0.05
            #model.Params.Method=2
            model.Params.MIPGap=.005
            #model.Params.lazyConstraints = 1
            if r >=20:
                model.Params.TimeLimit = 200
            else:
                model.Params.TimeLimit = 200
            model._lowerbounds = max(lowerbound.values())
            model.optimize(mycallback)
            elapsed_time = time.time() - start_time
            TimeLapse[r] = elapsed_time
            print("\n")
            s_var={}
            z_var={}
            for i in PartType:
                for j in JobTypeData[i]:
                    for k in FeasibleMachines[i,j]:
                        for t in range(1,T+1):
                            if (i,j,k,t) in s:
                              s_var[i,j,k,t] = s[i,j,k,t].X
                            if (i,j,k,t) in z:
                                z_var[i,j,k,t]= z[i,j,k,t].X
#            dd = ((1-Theta)/A)*sum(QualCost[i,j,k]*z_var[i,j,k,t] for i in PartType for j in JobTypeData[i] for k in FeasibleMachines[i,j] for t in range(1,T+1) if (i,j,k,t) in z_var)
#            sf=  sum(o[t].X for t in range(1,T+1))
#            lowerbound[r]= model.objbound
#            print("\nLower Bound = %f" %lowerbound[r])
#            print("Upper Bound = %f" % (dd+Stage2Total))
#            UpperBound[r] = (dd+Stage2Total)
#            if dd >0:
#                storeQual[r] = sum(z_var.values())
    
    BestLowerBounds[iter] = max(lowerbound.values())
    BestUpperBounds[iter] = min(UpperBound.values())
    
    

#
# Copyright 2020, Gurobi Optimization, LLC
#
# Interactive shell customization example
#
# Define a set of customizations for the Gurobi shell.
# Type 'from custom import *' to import them into your shell.
#
# Custom termination criterion: Quit optimization
# - after  5s if a high quality (1% gap) solution has been found, or
# - after 10s if a feasible solution has been found.

#%%
# Number of job types which are above \tau
counters=dict((key,0) for key in range(1,T+1))
totaljobs= dict((key2,0) for key2 in range(1,T+1))
for i in PartType:
    for j in JobTypeData[i]:
        for t in range(1,T+1):
           if Demand[i,j,t] >0: 
            totaljobs[t]=totaljobs[t]+1
            if sum(s[i,j,k,t].X for k in FeasibleMachines[i,j]) > Tau:
              counters[t]= counters[t]+1
              
              
ratios = {}
for t in range(1,T+1):
    ratios[t] = counters[t]/totaljobs[t]    

Utilization = dict(((key1,key3), 0) for key1 in Work_Center_ID for key3 in range(1,T+1)) 
   
above = {}
store=[]
counts=0
for k in Work_Center_ID:
    for t in range(1,T+1):
        Utilization[k,t] = round(sum((1/Capacity[k,t])*Pijk.select(i,j,k,'*')[0][3]*x_var[i,j,k,t] for i in PartType for j in JobTypeData[i] if (i,j,k,t) in x_var),2)
        if Utilization[k,t] > 1.0:
            above[k,t] = 1
            store.append(k)
            counts+=1   
maxUtlization = dict(((key1), 0) for key1 in range(1,T+1))

for t in range(1,T+1):
    for k in Work_Center_ID:
    
        if Utilization[k,t] > maxUtlization[t]:
            maxUtlization[t]= Utilization[k,t]
            maxUtlization[t] = max(maxUtlization[t] -0.7,0)
# chech demand satisfaction
            
for i in PartType:
    for j in JobTypeData[i]:
            if Demand[i,j,t]!=0:
                d =sum(x_var[i,j,k,t] for k in FeasibleMachines[i,j])
                if round(d) < round(Demand[i,j,t]):
                    print(i,j,t)
                    print(Demand[i,j,t])
                    print(d)
                    break

#%%
import matplotlib.pyplot as plt
lowerbounds_X=[]
lowerbounds_Y=[]
upperbounds_Y = []
for ss in range(1,r):
    lowerbounds_X.append(TimeLapse[ss])
    lowerbounds_Y.append(lowerbound[ss])
    upperbounds_Y.append(UpperBound[ss])
    
plt.plot(lowerbounds_X, lowerbounds_Y, 'r--', lowerbounds_X, upperbounds_Y, 'bs--')    
plt.xlabel('time (seconds)')
Draws_Y =[]
minUpp=100
counts=0
Draws_X = []
for iter in range(1,r):
    if UpperBound[iter] <= minUpp:
        counts=counts+1
        minUpp = UpperBound[iter]
        Draws_Y.append(minUpp)
        Draws_X.append(TimeLapse[counts])
        #plt.plot(Draws_X, Draws_Y, 'b',label="Bend-ub")
    else:
        counts=counts+1
        Draws_Y.append(minUpp)
        Draws_X.append(TimeLapse[counts])
plt.plot(Draws_X, Draws_Y, 'b',label="Bend-ub")
        
Draws_Y =[]
maxLow=0
counts=0
Draws_X = []
for iter in range(1,r):
    if lowerbound[iter] >= maxLow:
        counts=counts+1
        maxLow = lowerbound[iter]
        Draws_Y.append(maxLow)
        Draws_X.append(TimeLapse[counts])
#        plt.plot(Draws_X, Draws_Y, 'r',label="Bend-lb") 
    else:
        counts=counts+1
        Draws_Y.append(maxLow)
        Draws_X.append(TimeLapse[counts])
plt.plot(Draws_X, Draws_Y, 'r',label="Bend-lb") 

Draws_Y =[0.0781326,0.0497663, 0.0491283, .04913, .04913, .04913]
maxLow=0
counts=0
Draws_X = [60.4642,349.206,358.426, 1431, 1659, 2200]     
plt.plot(Draws_X, Draws_Y, 'x-',label="Gur-ub")    
Draws_Y =[0,0.046357, 0.046357, .04762, .04782, .04795]
maxLow=0
counts=0
Draws_X = [45,349.206,358.426, 1431, 1659, 2200]     
plt.plot(Draws_X, Draws_Y, marker='^', label="Gur-lb") 
plt.legend(loc='upper right')

#%%
lowerbounds_X=[]
lowerbounds_Y=[]
upperbounds_Y = []
for ss in range(1,r):
    lowerbounds_X.append(TimeLapse[ss])
    lowerbounds_Y.append(lowerbound[ss])
    upperbounds_Y.append(UpperBound[ss])
    
plt.plot(lowerbounds_X, lowerbounds_Y, 'r--', lowerbounds_X, upperbounds_Y, 'bs--')    
plt.xlabel('time (seconds)')
Draws_Y =[]
minUpp=100
counts=0
Draws_X = []
for iter in range(1,r):
    if UpperBound[iter] <= minUpp:
        counts=counts+1
        minUpp = UpperBound[iter]
        Draws_Y.append(minUpp)
        Draws_X.append(TimeLapse[counts])
        plt.plot(Draws_X, Draws_Y, 'b')
    else:
        counts=counts+1
        Draws_Y.append(minUpp)
        Draws_X.append(TimeLapse[counts])
        plt.plot(Draws_X, Draws_Y, 'b')
        
Draws_Y =[]
maxLow=0
counts=0
Draws_X = []
for iter in range(1,r):
    if lowerbound[iter] >= maxLow:
        counts=counts+1
        maxLow = lowerbound[iter]
        Draws_Y.append(maxLow)
        Draws_X.append(TimeLapse[counts])
        plt.plot(Draws_X, Draws_Y, 'r') 
    else:
        counts=counts+1
        Draws_Y.append(maxLow)
        Draws_X.append(TimeLapse[counts])
        plt.plot(Draws_X, Draws_Y, 'r') 

