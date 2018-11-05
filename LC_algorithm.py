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
        print('#####')
        start_time = time.time() 
        lc_alg.fm.reset(0)  # reset the follower model
        # get the hired resources from the leader -> y_ka
        for k in range(rcpsp.n_resources):
            y_ka.append(round(model.cbGetSolution(model._y_k[k])))
        # Get the project duration from HPS:
        f_hps = None
        for t in range(rcpsp.T):
            if round(model.cbGetSolution(model._v_ti[(t, rcpsp.n_activities - 1)])) == 1:
                f_hps = t
                break
        # Delete previous resource constraints from the model
        if lc_alg.resource_constraints is not None:
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
        print("Obj val: %f" % model.cbGet(g.GRB.Callback.MIPSOL_OBJ))
        print("HPS: %d, BFS: %d" % (f_hps, f_bfs))
        print("Hired resources: " + str(y_ka))
        #lc_alg.fm.write("%f_file.lp" % start_time)
        # Get the bit representation of the hired resources
        y_ka_b = [''.join(reversed(format(y_ka[k], 'b').zfill(8))) for k in range(rcpsp.n_resources)]
        b_0 = []
        b_1 = []
        for k in range(rcpsp.n_resources):
            b_0.append([i for i, x in enumerate(y_ka_b[k]) if x == '0'])
            b_1.append([i for i, x in enumerate(y_ka_b[k]) if x == '1'])
        if f_hps  < f_bfs:
            lc_alg.fm.write("Error.lp")
            input("ERROR!")
        # Add lazy contstraint when the follower project legth is not equal the leader's
        if f_hps != f_bfs:
            model.cbLazy(1 - (g.quicksum([1 - model._y_kb[(k, b)] for k in range(rcpsp.n_resources) for b in b_1[k]]) + g.quicksum([model._y_kb[(k,b)] for k in range(rcpsp.n_resources) for b in b_0[k]])) <= model._v_ti[(f_bfs, rcpsp.n_activities - 1)])
            #DEBUG
            start_times_hpm = []
            for i in range(rcpsp.n_activities):
                for t in range(rcpsp.T):
                    if round(model.cbGetSolution(model._v_ti[(t, i)])) == 1:
                        start_times_hpm.append(t)
                        break
            z_sum = 0
            zs = []
            gs = []
            for t in range(rcpsp.T):
                for k in range(rcpsp.n_resources):
                    zs.append(round(model.cbGetSolution(model._z_tk[(t,k)])))
            for t in range(rcpsp.T):
                for k in range(rcpsp.n_resources):
                    gs.append(round(model.cbGetSolution(model._g_tk[(t,k)])))
            z_sum = sum(zs)        
            start_times_fol = []
            for i in range(rcpsp.n_activities):
                for t in range(rcpsp.T):
                    if round(lc_alg.fm._v_ti[(t, i)].x) == 1:
                        start_times_fol.append(t)
                        break
            model._n_lc += 1
        end_time = time.time()
        print("LC_time %f", end_time-start_time)
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
        #self.y_ka = [None] * self.rcpsp.n_resources #
        self.v_ti = self.hpm.addVars(self.rcpsp.T, self.rcpsp.n_activities, vtype = g.GRB.BINARY, name = 'v_ti') #
        self.y_k = self.hpm.addVars(self.rcpsp.n_resources, vtype = g.GRB.CONTINUOUS, lb = 0, name= 'y_k')
        self.y_kb = self.hpm.addVars(rcpsp.n_resources, 8, vtype = g.GRB.BINARY, name = 'y_kb') # TODO: 8 -> RCPSP.B
        self.g_tk = self.hpm.addVars(self.rcpsp.T, self.rcpsp.n_resources, vtype = g.GRB.INTEGER, name = 'g_tk', lb = 0)
        self.z_tk = self.hpm.addVars(self.rcpsp.T, self.rcpsp.n_resources, vtype = g.GRB.CONTINUOUS, name = "z_tk", lb = 0)
        # Set the objective
        self.hpm.setObjective(self.rcpsp.alfa * g.quicksum([self.z_tk[(t, k)] for t in range(self.rcpsp.T) for k in range(self.rcpsp.n_resources)]) \
                       + g.quicksum([self.rcpsp.resource_cost[k] * self.y_k[(k)] for k in range(self.rcpsp.n_resources)]) \
                       #+ 0.01 * g.quicksum([t * self.v_ti[(t, rcpsp.n_activities-1)] for t in range(rcpsp.T)])
                       , g.GRB.MINIMIZE)
        # Set the constraints
        # 1. Every activity must start
        self.hpm.addConstrs(g.quicksum([self.v_ti[(t, i)] for t in range(self.rcpsp.T)]) == 1 for i in range(self.rcpsp.n_activities))
        # 2. Precedence relations
        self.hpm.addConstrs((g.quicksum([t * self.v_ti[(t, j)] for t in range(self.rcpsp.T)]) \
                     - g.quicksum([t* self.v_ti[(t, i)] for t in range(self.rcpsp.T)]) \
                     >= self.rcpsp.activity_duration[i] \
                     for i in range(self.rcpsp.n_activities) for j in self.rcpsp.activity_successor[i]), name = 'precedence')
        # (Optional): set start times before the earliest start time of an activity to 0
        self.hpm.addConstrs(self.v_ti[(t, i)] == 0 for i in range(self.rcpsp.n_activities) for t in range(0, self.rcpsp.activity_est[i]))
        # 3. Bit representation of the y_k variables
        self.hpm.addConstrs(self.y_k[k] == g.quicksum([(2 ** b) * self.y_kb[(k, b)]
                for b in range(8)]) for k in range(self.rcpsp.n_resources)) #TODO: 8 -> RCPSP.B
        # 4. The resource requirements per resource type per time unit
        for t in range(self.rcpsp.T):
            for k in range(self.rcpsp.n_resources):
                self.hpm.addConstr(g.quicksum([self.rcpsp.activity_resource[i][k] * self.v_ti[(t_, i)]
                                for i in range(self.rcpsp.n_activities)
                                for t_ in range(max(0, t - self.rcpsp.activity_duration[i] + 1), t+1)])
                           + self.g_tk[(t,k)] == self.rcpsp.resource_available[k] + self.y_k[k], name = 'resources_t=%d,k=%d' %(t, k))
        for t in range(self.rcpsp.T):
            for k in range(self.rcpsp.n_resources):
                self.hpm.addConstr(self.g_tk[(t,k)] <= self.z_tk[(t,k)] + self.M * g.quicksum([self.v_ti[(t_, rcpsp.n_activities-1)] for t_ in range(0, t+1)]), name = "z_t%d_k%d" % (t, k))
        for k in range(self.rcpsp.n_resources):
            self.hpm.addConstr(self.y_k[k] <= 250)
        self.hpm.Params.lazyConstraints = 1
        self.hpm.setParam( 'OutputFlag', False )s
        self.hpm.update()
        self.hpm._lc_time = 0
        self.hpm._n_lc = 0
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
        start_time_init = time.time() 
        # Get the smallest possible project makespan
        self.fm.optimize()
        self.shortest_makespan = self.getProjectMakespan(self.fm._v_ti)
        self.fm.reset(1)
        # Get the upper bounds for resources by solving resource availibility model
        if "1" in self.tricks_used:
            y_ka = self.fm.addVars(self.rcpsp.n_resources, vtype = g.GRB.CONTINUOUS, lb = 0, name="y_ka")
            self.fm.setObjective(g.quicksum([y_ka[k]*self.rcpsp.resource_cost[k] for k in range(self.rcpsp.n_resources)]), g.GRB.MINIMIZE)
            resource_constraints = self.fm.addConstrs((g.quicksum([rcpsp.activity_resource[i][k] * self.fm._v_ti[(t_, i)]
                            for i in range(rcpsp.n_activities)
                            for t_ in range(max(0, t - rcpsp.activity_duration[i] + 1), t + 1) ])
                       <= rcpsp.resource_available[k] + y_ka[k]
                       for t in range(rcpsp.T) for k in range(self.rcpsp.n_resources)), name = 'res_req')
            res_av_deadline = self.fm.addConstr(self.fm._v_ti[(self.shortest_makespan, self.rcpsp.n_activities-1)] == 1)
            self.fm.optimize()
            self.fm.write("res_av.lp")
            self.resource_ub = []
            for k in range(self.rcpsp.n_resources):
                self.resource_ub.append(round(y_ka[k].x))
            self.hpm.addConstrs(self.y_k[k] <= round(y_ka[k].x) for k in range(self.rcpsp.n_resources))
            # Adjust the model so it would correspond again with the follower
            self.fm.remove(y_ka)
            self.fm.remove(resource_constraints)
            self.fm.remove(res_av_deadline)
            self.fm.setObjective(g.quicksum([t * self.fm._v_ti[(t, self.rcpsp.n_activities-1)] for t in range(self.rcpsp.T)]), g.GRB.MINIMIZE)
            self.fm.reset(0)
        else:
            self.resource_ub = ["NA"] * self.rcpsp.n_resources
        end_time_init = time.time() 
        self.init_time = end_time_init - start_time_init
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
        self.hpm._y_k = self.y_k
        self.hpm._y_kb = self.y_kb
        self.hpm._v_ti = self.v_ti
        self.hpm._z_tk = self.z_tk
        self.hpm._g_tk = self.g_tk
        self.hpm._lc_num = 0
        self.hpm.update()
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
                print("shortest possible makespan %d" % self.shortest_makespan)
                print("project makespan: %d" % t)
                break
        print("resource ub: " + str(self.resource_ub))
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
            self.solution_resource_hired.append(round(self.y_k[k].x))
        self.solution_ztk = np.ones((self.rcpsp.T, self.rcpsp.n_resources), dtype = int) * -1 
        for t in range(self.rcpsp.T):
            for k in range(self.rcpsp.n_resources):
                self.solution_ztk[t, k] = round(self.z_tk[(t, k)].x)
        self.solution_zsum = np.sum(self.solution_ztk)
        self.processing_time = self.end_time - self.start_time
        print("Time elapsed %.2f s" % self.processing_time)
        print("Time init %.2f s" % self.init_time)
        print("LC time %.2f s" % self.hpm._lc_time)
        print("LC # %d" % self.hpm._n_lc)

if __name__ == '__main__':
    dataset_folder = "datasets/patterson/"
    nproblems_solve = 51
    TEST_RUN = True
    TRICKS_USED = "12 "
    stat_dict = {
            "alfa": [],
            "problem_name": [],
            "n_activities": [],
            "n_edges": [],
            "resource_ub": [],
            "sol_resources_hired": [],
            "shortest_makespan": [],
            "sol_project_makespan": [],
            "processing_time" : [],
            "init_time": [],
            "n_lc": [],
            "lc_time": []
            }
    
    if TEST_RUN:
        #problem_name = "datasets/patterson/pat8.rcp"
        problem_name = "test_data.rcp"
        dataset_folder = ""
        print("xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx")
        print(problem_name)
        rcpsp = RCPSP(dataset_folder + problem_name, alfa = 1000, resource_cost = 1000)
        lc = LC_alg(rcpsp, tricks_used = TRICKS_USED)
        lc.solve()
        stat_dict["alfa"].append(100)
        stat_dict["problem_name"].append(problem_name)
        stat_dict["n_activities"].append(rcpsp.n_activities)
        stat_dict["n_edges"].append(rcpsp.n_edges)
        stat_dict["resource_ub"].append(lc.resource_ub)
        stat_dict["sol_resources_hired"].append(lc.solution_resource_hired)
        stat_dict["shortest_makespan"].append(lc.shortest_makespan)
        stat_dict["sol_project_makespan"].append(lc.solution_start_times[rcpsp.n_activities-1])
        stat_dict["processing_time"].append(lc.processing_time)
        stat_dict["init_time"].append(lc.init_time)
        stat_dict["n_lc"].append(lc.hpm._n_lc)
        stat_dict["lc_time"].append(lc.hpm._lc_time)
    else:
        for alfa in [1000]:
            for i in range(1, min(nproblems_solve, len(os.listdir(dataset_folder)))):
                problem_name = "pat%d.rcp" % i 
                print("xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx")
                print(problem_name)                
                rcpsp = RCPSP(dataset_folder + problem_name, alfa = alfa, resource_cost = 1000)
                lc = LC_alg(rcpsp, tricks_used = TRICKS_USED)
                lc.solve()
                stat_dict["alfa"].append(alfa)
                stat_dict["problem_name"].append(problem_name)
                stat_dict["n_activities"].append(rcpsp.n_activities)
                stat_dict["n_edges"].append(rcpsp.n_edges)
                stat_dict["resource_ub"].append(lc.resource_ub)
                stat_dict["sol_resources_hired"].append(lc.solution_resource_hired)
                stat_dict["shortest_makespan"].append(lc.shortest_makespan)
                stat_dict["sol_project_makespan"].append(lc.solution_start_times[rcpsp.n_activities-1])
                stat_dict["processing_time"].append(lc.processing_time)
                stat_dict["init_time"].append(lc.init_time)
                stat_dict["n_lc"].append(lc.hpm._n_lc)
                stat_dict["lc_time"].append(lc.hpm._lc_time)
    df = pd.DataFrame(stat_dict)
    df=df[["alfa","problem_name", "n_activities", "n_edges", "resource_ub", "sol_resources_hired", "shortest_makespan", "sol_project_makespan", "processing_time", "init_time", "n_lc", "lc_time"]]
    df.to_csv("x.csv", sep = ";")
         
        
        