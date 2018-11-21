# -*- coding: utf-8 -*-
"""
Created on Thu Nov 15 18:28:01 2018

@author: pmilicka
"""
import gurobipy as g
from RCPSP import RCPSP
import time

class Seq_algorithm():
    
        def __init__(self, rcpsp):
            self.solutions = {"racp_solution": [],
                              "deadlines": [],
                              "activity_start": [],
                              "objective_values": [],
                              "objective_value_total": [],
                              "runtimes": [],
                              "bilevel_feasible": []}


            self.M = 1000
            self.rcpsp = rcpsp   # problem instance
            '''Init the RACP'''
            self.m1 = g.Model('') # high point solution model
            self.m1._v_ti = self.m1.addVars(self.rcpsp.T, self.rcpsp.n_activities, vtype = g.GRB.BINARY, name = 'v_ti') #
            self.m1._y_k = self.m1.addVars(self.rcpsp.n_resources, vtype = g.GRB.INTEGER, lb = 0, name= 'y_k')
            # Set the objective
            self.m1.setObjective(g.quicksum([self.rcpsp.resource_cost[k] * self.m1._y_k[(k)] for k in range(self.rcpsp.n_resources)]) \
                           , g.GRB.MINIMIZE)
            # Set the constraints
            # 1. Every activity must start
            self.m1.addConstrs(g.quicksum([self.m1._v_ti[(t, i)] for t in range(self.rcpsp.T)]) == 1 for i in range(self.rcpsp.n_activities))
            # 2. Precedence relations
            self.m1.addConstrs((g.quicksum([t * self.m1._v_ti[(t, j)] for t in range(self.rcpsp.T)]) \
                         - g.quicksum([t* self.m1._v_ti[(t, i)] for t in range(self.rcpsp.T)]) \
                         >= self.rcpsp.activity_duration[i] \
                         for i in range(self.rcpsp.n_activities) for j in self.rcpsp.activity_successor[i]), name = 'precedence')
            self.m1.addConstr(self.m1._v_ti[(0,0)] == 1)
            # 4. The resource requirements per resource type per time unit
            for t in range(self.rcpsp.T):
                for k in range(self.rcpsp.n_resources):
                    self.m1.addConstr(g.quicksum([self.rcpsp.activity_resource[i][k] * self.m1._v_ti[(t_, i)]
                                    for i in range(self.rcpsp.n_activities)
                                    for t_ in range(max(0, t - self.rcpsp.activity_duration[i] + 1), t+1)])
                                <= self.rcpsp.resource_available[k] + self.m1._y_k[k], name = 'resources_t=%d,k=%d' %(t, k))
    
            self.m2 = g.Model('')
            self.m2._v_ti = self.m2.addVars(self.rcpsp.T, self.rcpsp.n_activities, vtype = g.GRB.BINARY, name = 'v_ti') #
            self.m2._j_tk_p = self.m2.addVars(list(range(1, self.rcpsp.T)), self.rcpsp.n_resources, vtype = g.GRB.CONTINUOUS, lb = 0, name ="j_tk_p")
            self.m2._j_tk_m = self.m2.addVars(list(range(1, self.rcpsp.T)), self.rcpsp.n_resources, vtype = g.GRB.CONTINUOUS, lb = 0, name ="j_tk_m")
            # Set the objective
            self.m2.setObjective(self.rcpsp.alfa * g.quicksum([self.m2._j_tk_p[(t, k)] + self.m2._j_tk_m[(t, k)] for t in range(1, self.rcpsp.T) for k in range(self.rcpsp.n_resources)])\
                           , g.GRB.MINIMIZE)
            # 1. Every activity must start
            self.m2.addConstrs(g.quicksum([self.m2._v_ti[(t, i)] for t in range(self.rcpsp.T)]) == 1 for i in range(self.rcpsp.n_activities))
            # 2. Precedence relations
            self.m2.addConstrs((g.quicksum([t * self.m2._v_ti[(t, j)] for t in range(self.rcpsp.T)]) \
                         - g.quicksum([t* self.m2._v_ti[(t, i)] for t in range(self.rcpsp.T)]) \
                         >= self.rcpsp.activity_duration[i] \
                         for i in range(self.rcpsp.n_activities) for j in self.rcpsp.activity_successor[i]), name = 'precedence')
            # 5. The definition of jumps
            for t in range(0, self.rcpsp.T-1):
                for k in range(self.rcpsp.n_resources):
                    self.m2.addConstr(g.quicksum([self.rcpsp.activity_resource[i][k] * self.m2._v_ti[(t_, i)]
                                    for i in range(self.rcpsp.n_activities)
                                    for t_ in range(max(0, t - self.rcpsp.activity_duration[i] + 1), t+1)]) + self.m2._j_tk_p[(t+1,k)] - self.m2._j_tk_m[(t+1,k)] \
                                                 <= g.quicksum([self.rcpsp.activity_resource[i][k] * self.m2._v_ti[(t_, i)]
                                    for i in range(self.rcpsp.n_activities)
                                    for t_ in range(max(0, t + 1- self.rcpsp.activity_duration[i] + 1), t+2)]) + self.M * g.quicksum([self.m2._v_ti[(t_, rcpsp.n_activities-1)] for t_ in range(0, t+2)])
                                                 , name = "jumps1_t%d_k%d" % (t, k))
                    self.m2.addConstr(self.M * g.quicksum([self.m2._v_ti[(t_, rcpsp.n_activities-1)] for t_ in range(0, t+2)]) \
                                       + g.quicksum([self.rcpsp.activity_resource[i][k] * self.m2._v_ti[(t_, i)]
                                    for i in range(self.rcpsp.n_activities)
                                    for t_ in range(max(0, t - self.rcpsp.activity_duration[i] + 1), t+1)]) + self.m2._j_tk_p[(t+1,k)] - self.m2._j_tk_m[(t+1,k)] \
                                                 >= g.quicksum([self.rcpsp.activity_resource[i][k] * self.m2._v_ti[(t_, i)]
                                    for i in range(self.rcpsp.n_activities)
                                    for t_ in range(max(0, t - self.rcpsp.activity_duration[i] + 2), t+2)])
                                                 , name = "jumps2_t%d_k%d" % (t, k))
            self.m1.setParam( 'OutputFlag', True )
            """ Init the BFC"""
            self.bfc = g.Model('Follower')
            self.bfc._v_ti = self.bfc.addVars(self.rcpsp.T, self.rcpsp.n_activities, vtype = g.GRB.BINARY, name = 'v_ti') #
            # Set the objective
            self.bfc.setObjective(g.quicksum([t * self.bfc._v_ti[(t, self.rcpsp.n_activities-1)] for t in range(self.rcpsp.T)]), g.GRB.MINIMIZE)
            # Set the constraints 
            # 1. Every activity must start
            self.bfc.addConstrs(g.quicksum([self.bfc._v_ti[(t, i)] for t in range(self.rcpsp.T)]) == 1 for i in range(self.rcpsp.n_activities))
            # 2. Precedence relations
            self.bfc.addConstrs((g.quicksum([t * self.bfc._v_ti[(t, j)] for t in range(self.rcpsp.T)]) \
                         - g.quicksum([t* self.bfc._v_ti[(t, i)] for t in range(self.rcpsp.T)]) \
                         >= self.rcpsp.activity_duration[i] \
                         for i in range(self.rcpsp.n_activities) for j in self.rcpsp.activity_successor[i]), name = 'precedence')
            # 4. The resource requirements per resource type per time unit -> these are added in the lazy contraint
            self.bfc.setParam('OutputFlag', False)
            
            
        
        
        def solve(self):
            m1_deadline_constr = None
            m2_resource_constr = None
            m2_deadline_constr = None
            m3_resource_constr = None
            for deadline in range(self.rcpsp.shortest_makespan, self.rcpsp.T):
                start_time = time.time()
                self.solutions["deadlines"].append(deadline)
                self.m1.reset()
                if m1_deadline_constr is not None:
                    self.m1.remove(m1_deadline_constr)
                m1_deadline_constr = self.m1.addConstr(self.m1._v_ti[(deadline, self.rcpsp.n_activities-1)] == 1)
                self.m1.optimize()
                # Get the resources from the racp
                resources = []
                for k in range(self.rcpsp.n_resources):
                    resources.append(round(self.m1._y_k[k].x))
                self.solutions["racp_solution"].append(resources)
                racp_objval = round(self.m1.objVal)
                # 
                self.m2.reset()
                if m2_resource_constr is not None:
                    self.m2.remove(m2_resource_constr)
                if m2_deadline_constr is not None:
                    self.m2.remove(m2_deadline_constr)
                m2_deadline_constr = self.m2.addConstr(self.m2._v_ti[(deadline, self.rcpsp.n_activities-1)] == 1)
                m2_resource_constr = self.m2.addConstrs((g.quicksum([self.rcpsp.activity_resource[i][k] * self.m2._v_ti[(t_, i)]
                                for i in range(self.rcpsp.n_activities)
                                for t_ in range(max(0, t - self.rcpsp.activity_duration[i] + 1), t+1)])
                            <= self.rcpsp.resource_available[k] + resources[k] for k in range(self.rcpsp.n_resources) for t in range(self.rcpsp.T)), name = 'resources')
                self.m2.optimize()
                self.bfc.reset()
                if m3_resource_constr is not None:
                    self.bfc.remove(m3_resource_constr)
                m3_resource_constr = self.bfc.addConstrs((g.quicksum([self.rcpsp.activity_resource[i][k] * self.bfc._v_ti[(t_, i)]
                                for i in range(self.rcpsp.n_activities)
                                for t_ in range(max(0, t - self.rcpsp.activity_duration[i] + 1), t+1)])
                            <= self.rcpsp.resource_available[k] + resources[k] for k in range(self.rcpsp.n_resources) for t in range(self.rcpsp.T)), name = 'resources')
                self.bfc.optimize()
                bfc_pd = -1
                for t in range(self.rcpsp.T):
                    if round(self.bfc._v_ti[(t, self.rcpsp.n_activities-1)].x) == 1:
                        bfc_pd = t
                        break
                if bfc_pd < deadline:
                    self.solutions["bilevel_feasible"].append("False")
                else:
                    self.solutions["bilevel_feasible"].append("True")
                end_time = time.time()
                self.solutions["runtimes"].append(end_time - start_time)
                activity_start = []
                for i in range(self.rcpsp.n_activities):
                   for t in range(self.rcpsp.T):
                       if round(self.m2._v_ti[(t, i)].x) == 1:
                           activity_start.append(t)
                           break
                reslev_objval = round(self.m2.objVal)
                self.solutions["activity_start"].append(activity_start)
                self.solutions["objective_value_total"].append(racp_objval + reslev_objval)
                self.solutions["objective_values"].append([racp_objval, reslev_objval])
        
        
        def getBestSolutionIdx(self):
            obj_val_best = float("inf")
            best_idx = -1
            for i in range(len(self.solutions["deadlines"])):
                if self.solutions["objective_value_total"][i] < obj_val_best:
                    best_idx = i
                    obj_val_best = self.solutions["objective_value_total"][i]
            return best_idx
                
                
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
#    stat_dict = {
#            "alfa": [],
#            "problem_name": [],
#            "zero_rk": [],
#            "one_resource_only": [],
#            "n_activities": [],
#            "n_edges": [],
#            "resource_ub": [],
#            "resource_ub_approach": [],
#            "sol_resources_hired": [],
#            "shortest_makespan": [],
#            "sol_project_makespan": [],
#            "processing_time" : [],
#            "init_time": [],
#            "model_status": [],
#            "nodes_explored": []
#            }
    
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
        lc = Seq_algorithm(rcpsp)
        lc.solve()
#        stat_dict["problem_name"].append(problem_name)
#        stat_dict["zero_rk"].append(str(ZERO_RK))
#        stat_dict["one_resource_only"].append(str(ONE_RESOURCE_ONLY))
#        stat_dict["n_activities"].append(rcpsp.n_activities)
#        stat_dict["n_edges"].append(rcpsp.n_edges)
#        stat_dict["resource_ub"].append(lc.rcpsp.resource_ub)
#        stat_dict["resource_ub_approach"].append(str(approach_rub))
#        stat_dict["sol_resources_hired"].append(lc.solution_resource_hired)
#        stat_dict["shortest_makespan"].append(lc.rcpsp.shortest_makespan)
#        stat_dict["sol_project_makespan"].append(lc.solution_start_times[rcpsp.n_activities-1])
#        stat_dict["processing_time"].append(lc.processing_time)
#        stat_dict["init_time"].append(lc.rcpsp.init_time)
#        stat_dict["n_lc1"].append(lc.hpm._n_lc1)
#        stat_dict["n_lc2"].append(lc.hpm._n_lc2)
#        stat_dict["lc_time"].append(lc.hpm._lc_time)
#        stat_dict["model_status"].append(lc.model_status)
#        stat_dict["nodes_explored"].append(lc.model_nodecount)
#    else:
#        for alfa in [250]:
#            for i in range(problem_start, problem_start + min(nproblems_solve, len(os.listdir(dataset_folder)))):
#                problem_name = "Pat%d.rcp" % i 
#                print("############################################################")
#                print(problem_name)
#                init_time_start = time.time()
#                rcpsp = RCPSP(dataset_folder + problem_name, 
#                      alfa = alfa, 
#                      resource_cost = 1000, 
#                      approach_rub = approach_rub, 
#                      zero_rk = ZERO_RK,
#                      one_resource_only = ONE_RESOURCE_ONLY)
#                init_time_end = time.time()
#                init_time = init_time_end - init_time_start
#                lc = LC_alg(rcpsp, tricks_used = TRICKS_USED)
#                lc.solve()
#                stat_dict["alfa"].append(alfa)
#                stat_dict["problem_name"].append(problem_name)
#                stat_dict["zero_rk"].append(str(ZERO_RK))
#                stat_dict["one_resource_only"].append(str(ONE_RESOURCE_ONLY))
#                stat_dict["n_activities"].append(rcpsp.n_activities)
#                stat_dict["n_edges"].append(rcpsp.n_edges)
#                stat_dict["resource_ub"].append(lc.rcpsp.resource_ub)
#                stat_dict["resource_ub_approach"].append(str(approach_rub))
#                stat_dict["sol_resources_hired"].append(lc.solution_resource_hired)
#                stat_dict["shortest_makespan"].append(lc.rcpsp.shortest_makespan)
#                stat_dict["sol_project_makespan"].append(lc.solution_start_times[rcpsp.n_activities-1])
#                stat_dict["processing_time"].append(lc.processing_time)
#                stat_dict["init_time"].append(lc.rcpsp.init_time)
#                stat_dict["n_lc1"].append(lc.hpm._n_lc1)
#                stat_dict["n_lc2"].append(lc.hpm._n_lc2)
#                stat_dict["lc_time"].append(lc.hpm._lc_time)
#                stat_dict["model_status"].append(lc.model_status)
#                stat_dict["nodes_explored"].append(lc.model_nodecount)
#    df = pd.DataFrame(stat_dict)
#    df=df[["alfa","problem_name", "zero_rk", "one_resource_only", "n_activities", "n_edges", 
#           "resource_ub", "resource_ub_approach", "sol_resources_hired", "shortest_makespan", 
#           "sol_project_makespan", "processing_time", "init_time", "n_lc1", "n_lc2", "lc_time", 
#           "model_status", "nodes_explored"]]
#    df.to_csv(dataset_folder + "y.csv", sep = ";")

