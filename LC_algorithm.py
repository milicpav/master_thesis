# -*- coding: utf-8 -*-
"""
Created on Sun Oct 21 17:36:33 2018

@author: pmilicka
"""
import gurobipy as g
from RCPSP import RCPSP
import numpy as np
import os
import time
import pandas as pd

'''
Lazy callback
'''
def bilevel_feasibility_callback(model, where):
    lc_alg = model._lc_alg 
    rcpsp = lc_alg.rcpsp
    y_ka = []
    if  where  == g.GRB.Callback.MIPSOL:
        #print('#####')
        start_time = time.time() 
        """ get the hired resources from the leader -> y_ka """
        for k in range(rcpsp.n_resources):
            y_ka.append(round(model.cbGetSolution(model._y_k[k])))
        """ Get the bit representation of the hired resources """
        y_ka_b = [''.join(reversed(format(y_ka[k], 'b').zfill(rcpsp.resource_nbits[k]))) for k in range(rcpsp.n_resources)]
        b_0 = []
        b_1 = []
        for k in range(rcpsp.n_resources):
            b_0.append([i for i, x in enumerate(y_ka_b[k]) if x == '0'])
            b_1.append([i for i, x in enumerate(y_ka_b[k]) if x == '1'])
        y_ka = np.array(y_ka, dtype = int)
        """Get the duration of the high point solution"""
        f_hps = None
        for t in range(rcpsp.T):
            if round(model.cbGetSolution(model._v_ti[(t, rcpsp.n_activities - 1)])) == 1:
                f_hps = t
                break
        """ Try to obtain bound on the project makespan by examining previous lc """
        bound = []
        for y_vec in model._br_dict.keys():
            y_vec = np.array(y_vec, dtype = int)
            if np.all(y_ka == y_vec):
                #print("Type XXX")
                return
            if np.all(y_ka >= y_vec): # all resources in current solution are bigger or equal the stored one
                bound.append(model._br_dict[tuple(y_vec)])
        if len(bound) > 0:
            bound_min = min(bound)
            if f_hps > bound_min: 
                model.cbLazy(1 - (g.quicksum([1 - model._y_kb[(k, b)] for k in range(rcpsp.n_resources) for b in b_1[k]]) + g.quicksum([model._y_kb[(k,b)] for k in range(rcpsp.n_resources) for b in b_0[k]])) \
                             <= g.quicksum([model._v_ti[(t, rcpsp.n_activities - 1)] for t in range(0, bound_min + 1)]))
                model._n_lc2 += 1
            #print("Type two implanted")
            return
        """ Solve the follower problem and retrieve the follower makespan"""
        lc_alg.fm.reset(0)  # reset the follower model    
        if lc_alg.resource_constraints is not None: # Delete previous resource constraints from the model
            lc_alg.fm.remove(lc_alg.resource_constraints)
            lc_alg.fm.update()
        # Add the resource constraints to the follower problem
        lc_alg.resource_constraints = lc_alg.fm.addConstrs((g.quicksum([rcpsp.activity_resource[i][k] * lc_alg.fm._v_ti[(t_, i)]
                            for i in range(rcpsp.n_activities)
                            for t_ in range(max(0, t - rcpsp.activity_duration[i] + 1), t + 1)])
                       <= rcpsp.resource_available[k] + y_ka[k]
                       for t in range(rcpsp.T) for k in range(rcpsp.n_resources)), name = 'res_req')
        # Optimize the follower problem and get the follower project duration -> f_bfs
        lc_alg.fm.optimize()
        f_bfs = None
        for t in range(rcpsp.T):
            if lc_alg.fm._v_ti[(t, rcpsp.n_activities-1)].x > 0.9:
                f_bfs = t
                break
        # Output intermediate results
        #print("Obj val: %f" % model.cbGet(g.GRB.Callback.MIPSOL_OBJ))
        #print("HPS: %d, BFS: %d" % (f_hps, f_bfs))
        #print("Hired resources: " + str(y_ka))
        
        if f_hps  < f_bfs:
            lc_alg.fm.write("Error.lp")
            input("ERROR!")
        """ Add lazy constraints if follower has shorter makespan """
        if f_hps != f_bfs:
            model.cbLazy(1 - (g.quicksum([1 - model._y_kb[(k, b)] for k in range(rcpsp.n_resources) for b in b_1[k]]) + g.quicksum([model._y_kb[(k,b)] for k in range(rcpsp.n_resources) for b in b_0[k]])) <= model._v_ti[(f_bfs, rcpsp.n_activities - 1)])
            model._br_dict[tuple(y_ka)] = f_bfs
            #for DEBUG purposes
            start_times_hpm = []
            for i in range(rcpsp.n_activities):
                for t in range(rcpsp.T):
                    if round(model.cbGetSolution(model._v_ti[(t, i)])) == 1:
                        start_times_hpm.append(t)
                        break  
            start_times_fol = []
            for i in range(rcpsp.n_activities):
                for t in range(rcpsp.T):
                    if round(lc_alg.fm._v_ti[(t, i)].x) == 1:
                        start_times_fol.append(t)
                        break
            model._n_lc1 += 1
        model._br_dict[tuple(y_ka)] = f_bfs
        end_time = time.time()
        #print("LC_time %f", end_time-start_time)
        model._lc_time += (end_time - start_time)
#    elif where == g.GRB.Callback.MIPNODE:
#        status = model.cbGet(g.GRB.Callback.MIPNODE_STATUS)
#        if status == g.GRB.OPTIMAL:
#            #print("mipnode optimal")
#            for t in range(rcpsp.T):
#                for i i n range(rcpsp.n_resources):
#                    #print(model.cbGetSolution(model._v_ti[(t, rcpsp.n_activities - 1)])) == 1
#                    if abs(round(model.cbGetNodeRel(model._v_ti[(t,i)])) - model.cbGetNodeRel(model._v_ti[(t,i)])) > model.Params.IntFeasTol:
#                        return
#            print("mipnode optimal is integer")
                    
        
    elif where == g.GRB.Callback.MESSAGE:
        # Message callback
        msg = model.cbGet(g.GRB.Callback.MSG_STRING)
        if any(name in msg for name in ['Lazy constraints:']):
            print("Lazy constr: ", int(msg.split(':')[1]))
    
    
class LC_alg:
    '''
    Initializes the LC_algorithm, prepares the High point model structure etc.
    '''
    def __init__(self, rcpsp, tricks_used = ""):
        self.tricks_used = tricks_used
        self.start_time = time.time()
        self.M = 1000
        self.hpm = g.Model('') # high point solution model
        self.rcpsp = rcpsp   # problem instance
        self.hpm._v_ti = self.hpm.addVars(self.rcpsp.T, self.rcpsp.n_activities, vtype = g.GRB.BINARY, name = 'v_ti') #
        self.hpm._y_k = self.hpm.addVars(self.rcpsp.n_resources, vtype = g.GRB.INTEGER, lb = 0, name= 'y_k')
#        self.hpm._y_kb = self.hpm.addVars(rcpsp.n_resources, 8, vtype = g.GRB.BINARY, name = 'y_kb') # TODO: 8 -> RCPSP.B
        self.hpm._y_kb = self.hpm.addVars([(k, b) for k in range(self.rcpsp.n_resources) for b in range(self.rcpsp.resource_nbits[k])], vtype = g.GRB.BINARY, name = "y_kb")
        #self.g_tk = self.hpm.addVars(self.rcpsp.T, self.rcpsp.n_resources, vtype = g.GRB.INTEGER, name = 'g_tk', lb = 0)
        self.hpm._j_tk_p = self.hpm.addVars(list(range(1, self.rcpsp.T)), self.rcpsp.n_resources, vtype = g.GRB.CONTINUOUS, lb = 0, name ="j_tk_p")
        self.hpm._j_tk_m = self.hpm.addVars(list(range(1, self.rcpsp.T)), self.rcpsp.n_resources, vtype = g.GRB.CONTINUOUS, lb = 0, name ="j_tk_m")
        # Set the objective
        self.hpm.setObjective(self.rcpsp.alfa * g.quicksum([self.hpm._j_tk_p[(t, k)] + self.hpm._j_tk_m[(t, k)] for t in range(1, self.rcpsp.T) for k in range(self.rcpsp.n_resources)]) \
                       + g.quicksum([self.rcpsp.resource_cost[k] * self.hpm._y_k[(k)] for k in range(self.rcpsp.n_resources)]) \
                       #+ 0.01 * g.quicksum([t * self.v_ti[(t, rcpsp.n_activities-1)] for t in range(rcpsp.T)])
                       , g.GRB.MINIMIZE)
        # Set the constraints
        # 1. Every activity must start
        self.hpm.addConstrs(g.quicksum([self.hpm._v_ti[(t, i)] for t in range(self.rcpsp.T)]) == 1 for i in range(self.rcpsp.n_activities))
        # 2. Precedence relations
        self.hpm.addConstrs((g.quicksum([t * self.hpm._v_ti[(t, j)] for t in range(self.rcpsp.T)]) \
                     - g.quicksum([t* self.hpm._v_ti[(t, i)] for t in range(self.rcpsp.T)]) \
                     >= self.rcpsp.activity_duration[i] \
                     for i in range(self.rcpsp.n_activities) for j in self.rcpsp.activity_successor[i]), name = 'precedence')
        # (Optional): set start times before the earliest start time of an activity to 0
        if "1" in self.tricks_used:
            self.hpm.addConstrs(self.hpm._v_ti[(t, i)] == 0 for i in range(self.rcpsp.n_activities) for t in range(0, self.rcpsp.activity_est[i]))
       #self.hpm.addConstr(self.v_ti[(0, 0)] == 1)
        # 3. Bit representation of the y_k variables
        self.hpm.addConstrs(self.hpm._y_k[k] == g.quicksum([(2 ** b) * self.hpm._y_kb[(k, b)]
                for b in range(self.rcpsp.resource_nbits[k])]) for k in range(self.rcpsp.n_resources)) #TODO: 8 -> RCPSP.B
        # 4. The resource requirements per resource type per time unit
        for t in range(self.rcpsp.T):
            for k in range(self.rcpsp.n_resources):
                self.hpm.addConstr(g.quicksum([self.rcpsp.activity_resource[i][k] * self.hpm._v_ti[(t_, i)]
                                for i in range(self.rcpsp.n_activities)
                                for t_ in range(max(0, t - self.rcpsp.activity_duration[i] + 1), t+1)])
                            <= self.rcpsp.resource_available[k] + self.hpm._y_k[k], name = 'resources_t=%d,k=%d' %(t, k))
        # Get the upper bounds for resources 
        self.hpm.addConstrs(self.hpm._y_k[k] <= round(self.rcpsp.resource_ub[k]) for k in range(self.rcpsp.n_resources))    
        # 5. The definition of jumps
        for t in range(0, self.rcpsp.T-1):
            for k in range(self.rcpsp.n_resources):
#                self.hpm.addConstr(self.g_tk[(t, k)] - self.j_tk_p[(t+1,k)] + self.j_tk_m[(t+1,k)] \
#                                             <= self.g_tk[(t+1, k)] + self.M * g.quicksum([self.v_ti[(t_, rcpsp.n_activities-1)] for t_ in range(0, t+2)])
#                                             , name = "jumps1_t%d_k%d" % (t, k))
#                self.hpm.addConstr(self.M * g.quicksum([self.v_ti[(t_, rcpsp.n_activities-1)] for t_ in range(0, t+2)]) \
#                                   + self.g_tk[(t, k)] - self.j_tk_p[(t+1,k)] + self.j_tk_m[(t+1,k)] \
#                                             >= self.g_tk[(t+1, k)]
#                                             , name = "jumps2_t%d_k%d" % (t, k))
                self.hpm.addConstr(g.quicksum([self.rcpsp.activity_resource[i][k] * self.hpm._v_ti[(t_, i)]
                                for i in range(self.rcpsp.n_activities)
                                for t_ in range(max(0, t - self.rcpsp.activity_duration[i] + 1), t+1)]) + self.hpm._j_tk_p[(t+1,k)] - self.hpm._j_tk_m[(t+1,k)] \
                                             <= g.quicksum([self.rcpsp.activity_resource[i][k] * self.hpm._v_ti[(t_, i)]
                                for i in range(self.rcpsp.n_activities)
                                for t_ in range(max(0, t + 1- self.rcpsp.activity_duration[i] + 1), t+2)]) + self.M * g.quicksum([self.hpm._v_ti[(t_, rcpsp.n_activities-1)] for t_ in range(0, t+2)])
                                             , name = "jumps1_t%d_k%d" % (t, k))
                self.hpm.addConstr(self.M * g.quicksum([self.hpm._v_ti[(t_, rcpsp.n_activities-1)] for t_ in range(0, t+2)]) \
                                   + g.quicksum([self.rcpsp.activity_resource[i][k] * self.hpm._v_ti[(t_, i)]
                                for i in range(self.rcpsp.n_activities)
                                for t_ in range(max(0, t - self.rcpsp.activity_duration[i] + 1), t+1)]) + self.hpm._j_tk_p[(t+1,k)] - self.hpm._j_tk_m[(t+1,k)] \
                                             >= g.quicksum([self.rcpsp.activity_resource[i][k] * self.hpm._v_ti[(t_, i)]
                                for i in range(self.rcpsp.n_activities)
                                for t_ in range(max(0, t - self.rcpsp.activity_duration[i] + 2), t+2)])
                                             , name = "jumps2_t%d_k%d" % (t, k))
        self.hpm.Params.lazyConstraints = 1
        self.hpm.setParam( 'OutputFlag', True )
        self.hpm.update()
        self.hpm._lc_time = 0
        self.hpm._n_lc1 = 0
        self.hpm._n_lc2 = 0
        # Define the follower problem without resource constraints and store it as a class attribute
        self.initFollowerModel()
    
    """
    Initializes the follower model structure 
    """    
    def initFollowerModel(self):
        self.fm = g.Model('Follower')
        self.fm._v_ti = self.fm.addVars(self.rcpsp.T, self.rcpsp.n_activities, vtype = g.GRB.BINARY, name = 'v_ti') #
        # Set the objective
        self.fm.setObjective(g.quicksum([t * self.fm._v_ti[(t, self.rcpsp.n_activities-1)] for t in range(self.rcpsp.T)]), g.GRB.MINIMIZE)
        # Set the constraints 
        # 1. Every activity must start
        self.fm.addConstrs(g.quicksum([self.fm._v_ti[(t, i)] for t in range(self.rcpsp.T)]) == 1 for i in range(self.rcpsp.n_activities))
        # 2. Precedence relations
        self.fm.addConstrs((g.quicksum([t * self.fm._v_ti[(t, j)] for t in range(self.rcpsp.T)]) \
                     - g.quicksum([t* self.fm._v_ti[(t, i)] for t in range(self.rcpsp.T)]) \
                     >= self.rcpsp.activity_duration[i] \
                     for i in range(self.rcpsp.n_activities) for j in self.rcpsp.activity_successor[i]), name = 'precedence')
        # (Optional): set start times before the earliest start time of an activity to 0
        if "2" in self.tricks_used:
            self.fm.addConstrs(self.fm._v_ti[(t, i)] == 0 for i in range(self.rcpsp.n_activities) for t in range(0, self.rcpsp.activity_est[i]))
        # 4. The resource requirements per resource type per time unit -> these are added in the lazy contraint
        self.resource_constraints = None
        self.fm.setParam('OutputFlag', False)
    
    def getProjectMakespan(self, v_ti_variables):
        for t in range(self.rcpsp.T):
            if v_ti_variables[(t, self.rcpsp.n_activities - 1)].x > 0.9:
                return t
    '''
    Solve the current High Point Solution
    '''
    def solve(self):
        self.hpm._lc_alg = self
        self.hpm._br_dict = {}
        self.hpm.update()
        self.hpm.setParam('TimeLimit', 30*60)
        self.hpm.optimize(bilevel_feasibility_callback)
        self.end_time = time.time()
        
        print("Problem instance solved: ")
        print("alpha:", self.rcpsp.alfa)
        print("resource cost: " + str(self.rcpsp.resource_cost))
        tmp = []
        for i in range(len(self.hpm._y_k)):
            tmp.append(round(self.hpm._y_k[i].x))
        for t in range(self.rcpsp.T):
            if self.hpm._v_ti[(t, self.rcpsp.n_activities-1)].x > 0.9:
                print("shortest possible makespan %d" % self.rcpsp.shortest_makespan)
                print("project makespan: %d" % t)
                break
        print("resource ub: " + str(self.rcpsp.resource_ub))
        print("resource hired: " + str(tmp))
        print("Time required to solve: " + str(self.hpm.Runtime))
        print("obj val: %.2f" % self.hpm.objVal)
        self.hpm.write("file.lp")
        self.solution_start_times = [] 
        for i in range(self.rcpsp.n_activities):
            for t in range(self.rcpsp.T):
                if self.hpm._v_ti[(t, i)].x > 0.9:
                    self.solution_start_times.append(t)
                    break
        self.solution_resource_hired = []
        if self.solution_start_times[self.rcpsp.n_activities-1] != self.getProjectMakespan(self.hpm._v_ti):
            raise Exception("Something is wrong")
        for k in range(self.rcpsp.n_resources):
            self.solution_resource_hired.append(round(self.hpm._y_k[k].x))
            
        self.solution_j_tk_p = np.zeros((self.rcpsp.T, self.rcpsp.n_resources), dtype = int)
        self.solution_j_tk_m = np.zeros((self.rcpsp.T, self.rcpsp.n_resources), dtype = int) 
        for t in range(1, self.rcpsp.T):
            for k in range(self.rcpsp.n_resources):
                self.solution_j_tk_p[t, k] = round(self.hpm._j_tk_p[(t, k)].x)
                self.solution_j_tk_m[t, k] = round(self.hpm._j_tk_m[(t, k)].x)
        #self.solution_zsum = np.sum(self.solution_ztk)
        self.processing_time = self.end_time - self.start_time
        if self.hpm.status == g.GRB.OPTIMAL:
            self.model_status = "Optimal"
        elif self.hpm.status == g.GRB.TIME_LIMIT:
            self.model_status = "Time Limit Reached"
        self.model_nodecount = self.hpm.getAttr("NodeCount")
        print("Model status: %s" % self.model_status)
        print("Time elapsed %.2f s" % self.processing_time)
        print("Time init %.2f s" % self.rcpsp.init_time)
        print("LC time %.2f s" % self.hpm._lc_time)
        print("LC1 # %d" % self.hpm._n_lc1)
        print("LC2 # %d" % self.hpm._n_lc2)

if __name__ == '__main__':
    dataset_folder = "datasets/patterson/"
    dataset_folder = "datasets/RG10/rg10_02/"
    problem_start = 1
    nproblems_solve = 100
    TEST_RUN = True
    TRICKS_USED = ""
    ZERO_RK = False
    ONE_RESOURCE_ONLY = False
    approach_rub = "sum"
    stat_dict = {
            "alfa": [],
            "problem_name": [],
            "zero_rk": [],
            "one_resource_only": [],
            "n_activities": [],
            "n_edges": [],
            "resource_ub": [],
            "resource_ub_approach": [],
            "sol_resources_hired": [],
            "shortest_makespan": [],
            "sol_project_makespan": [],
            "processing_time" : [],
            "init_time": [],
            "n_lc1": [],
            "n_lc2": [],
            "lc_time": [],
            "model_status": [],
            "nodes_explored": [],
            "objval": []
            }
    
    if TEST_RUN:
        problem_name = "datasets/patterson/pat8.rcp"
        problem_name = "datasets/test_data.rcp"
        #dataset_folder = ""
        print("xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx")
        print(problem_name)
        rcpsp = RCPSP(problem_name, 
                      alfa = 250, 
                      resource_cost = 1000, 
                      approach_rub = approach_rub, 
                      zero_rk = ZERO_RK,
                      one_resource_only = ONE_RESOURCE_ONLY)
        lc = LC_alg(rcpsp, tricks_used = TRICKS_USED)
        lc.solve()
        stat_dict["alfa"].append(250)
        stat_dict["problem_name"].append(problem_name)
        stat_dict["zero_rk"].append(str(ZERO_RK))
        stat_dict["one_resource_only"].append(str(ONE_RESOURCE_ONLY))
        stat_dict["n_activities"].append(rcpsp.n_activities)
        stat_dict["n_edges"].append(rcpsp.n_edges)
        stat_dict["resource_ub"].append(lc.rcpsp.resource_ub)
        stat_dict["resource_ub_approach"].append(str(approach_rub))
        stat_dict["sol_resources_hired"].append(lc.solution_resource_hired)
        stat_dict["shortest_makespan"].append(lc.rcpsp.shortest_makespan)
        stat_dict["sol_project_makespan"].append(lc.solution_start_times[rcpsp.n_activities-1])
        stat_dict["processing_time"].append(lc.processing_time)
        stat_dict["init_time"].append(lc.rcpsp.init_time)
        stat_dict["n_lc1"].append(lc.hpm._n_lc1)
        stat_dict["n_lc2"].append(lc.hpm._n_lc2)
        stat_dict["lc_time"].append(lc.hpm._lc_time)
        stat_dict["model_status"].append(lc.model_status)
        stat_dict["nodes_explored"].append(lc.model_nodecount)
        stat_dict["objval"].append(lc.hpm.objVal)
        
    else:
        for alfa in [250]:
            for i in range(problem_start, problem_start + min(nproblems_solve, len(os.listdir(dataset_folder)))):
                problem_name = "Pat%d.rcp" % i 
                print("############################################################")
                print(problem_name)
                init_time_start = time.time()
                rcpsp = RCPSP(dataset_folder + problem_name, 
                      alfa = alfa, 
                      resource_cost = 1000, 
                      approach_rub = approach_rub, 
                      zero_rk = ZERO_RK,
                      one_resource_only = ONE_RESOURCE_ONLY)
                init_time_end = time.time()
                init_time = init_time_end - init_time_start
                lc = LC_alg(rcpsp, tricks_used = TRICKS_USED)
                lc.solve()
                stat_dict["alfa"].append(alfa)
                stat_dict["problem_name"].append(problem_name)
                stat_dict["zero_rk"].append(str(ZERO_RK))
                stat_dict["one_resource_only"].append(str(ONE_RESOURCE_ONLY))
                stat_dict["n_activities"].append(rcpsp.n_activities)
                stat_dict["n_edges"].append(rcpsp.n_edges)
                stat_dict["resource_ub"].append(lc.rcpsp.resource_ub)
                stat_dict["resource_ub_approach"].append(str(approach_rub))
                stat_dict["sol_resources_hired"].append(lc.solution_resource_hired)
                stat_dict["shortest_makespan"].append(lc.rcpsp.shortest_makespan)
                stat_dict["sol_project_makespan"].append(lc.solution_start_times[rcpsp.n_activities-1])
                stat_dict["processing_time"].append(lc.processing_time)
                stat_dict["init_time"].append(lc.rcpsp.init_time)
                stat_dict["n_lc1"].append(lc.hpm._n_lc1)
                stat_dict["n_lc2"].append(lc.hpm._n_lc2)
                stat_dict["lc_time"].append(lc.hpm._lc_time)
                stat_dict["model_status"].append(lc.model_status)
                stat_dict["nodes_explored"].append(lc.model_nodecount)
                stat_dict["objval"].append(lc.hpm.objVal)
    df = pd.DataFrame(stat_dict)
    df=df[["alfa","problem_name", "zero_rk", "one_resource_only", "n_activities", "n_edges", 
           "resource_ub", "resource_ub_approach", "sol_resources_hired", "shortest_makespan", 
           "sol_project_makespan", "processing_time", "init_time", "n_lc1", "n_lc2", "lc_time", 
           "model_status", "nodes_explored", "objval"]]
    df.to_csv(dataset_folder + "x.csv", sep = ";")
         
        
        