# -*- coding: utf-8 -*-
"""
Created on Fri Nov 23 09:19:51 2018

@author: pmilicka
"""
from RCPSP import RCPSP
from LC_algorithm import LC_alg
import gurobipy as g
import sys

class IntegratedApproach:
    def __init__(self, rcpsp):
        self.M = 1000
        self.hpm = g.Model('Integrated') # high point solution model
        self.rcpsp = rcpsp   # problem instance
        self.hpm._v_ti = self.hpm.addVars(self.rcpsp.T, self.rcpsp.n_activities, vtype = g.GRB.BINARY, name = 'v_ti') #
        self.hpm._y_k = self.hpm.addVars(self.rcpsp.n_resources, vtype = g.GRB.INTEGER, lb = 0, name= 'y_k')
        self.hpm._j_tk_p = self.hpm.addVars(list(range(1, self.rcpsp.T)), self.rcpsp.n_resources, vtype = g.GRB.CONTINUOUS, lb = 0, name ="j_tk_p")
        self.hpm._j_tk_m = self.hpm.addVars(list(range(1, self.rcpsp.T)), self.rcpsp.n_resources, vtype = g.GRB.CONTINUOUS, lb = 0, name ="j_tk_m")
        # Set the objective
        # Set the constraints
        # 1. Every activity must start
        self.hpm.addConstrs(g.quicksum([self.hpm._v_ti[(t, i)] for t in range(self.rcpsp.T)]) == 1 for i in range(self.rcpsp.n_activities))
        # 2. Precedence relations
        self.hpm.addConstrs((g.quicksum([t * self.hpm._v_ti[(t, j)] for t in range(self.rcpsp.T)]) \
                     - g.quicksum([t* self.hpm._v_ti[(t, i)] for t in range(self.rcpsp.T)]) \
                     >= self.rcpsp.activity_duration[i] \
                     for i in range(self.rcpsp.n_activities) for j in self.rcpsp.activity_successor[i]), name = 'precedence')
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
        #self.hpm.addConstr(self.hpm._v_ti[(21, self.rcpsp.n_activities-1)] == 1)
        self.hpm.setParam('OutputFlag', False)
    
    def solve(self, alfa, beta, gamma):
        self.hpm.reset(clearall = 1)
        self.hpm.setObjective(alfa * g.quicksum([self.hpm._j_tk_p[(t, k)] + self.hpm._j_tk_m[(t, k)] for t in range(1, self.rcpsp.T) for k in range(self.rcpsp.n_resources)]) \
               + g.quicksum([self.rcpsp.resource_cost[k] * self.hpm._y_k[(k)] for k in range(self.rcpsp.n_resources)]) \
               + gamma * g.quicksum([t * self.hpm._v_ti[(t, self.rcpsp.n_activities-1)] for t in range(self.rcpsp.T)])
               , g.GRB.MINIMIZE)
        self.hpm.optimize()
        # get project makespan
        sol_makespan = -1
        for t in range(self.rcpsp.T):
            if round(self.hpm._v_ti[(t, self.rcpsp.n_activities-1)].x) == 1:
                sol_makespan = t
                break
        # get n.jumps
        sol_jumps = 0
        for t in range(1, self.rcpsp.T):
            for k in range(self.rcpsp.n_resources):
                sol_jumps += round(self.hpm._j_tk_p[(t, k)].x)
                sol_jumps += round(self.hpm._j_tk_m[(t, k)].x)
        # get resources_hired
        sol_resources = 0
        for k in range(self.rcpsp.n_resources):
           sol_resources += self.hpm._y_k[k].x
        return sol_makespan, sol_resources, sol_jumps
    
    def getSolutionStartTimes(self):
        solution_start_times = [] 
        for i in range(self.rcpsp.n_activities):
            for t in range(self.rcpsp.T):
                if self.hpm._v_ti[(t, i)].x > 0.9:
                    solution_start_times.append(t)
                    break
        return solution_start_times
    
    def checkBilevelFeasibility(self):
        resources = []
        for k in range(self.rcpsp.n_resources):
            resources.append(round(self.hpm._y_k[k].x))
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
        self.bfc.addConstrs((g.quicksum([self.rcpsp.activity_resource[i][k] * self.bfc._v_ti[(t_, i)]
                                for i in range(self.rcpsp.n_activities)
                                for t_ in range(max(0, t - self.rcpsp.activity_duration[i] + 1), t+1)])
                            <= self.rcpsp.resource_available[k] + resources[k] for k in range(self.rcpsp.n_resources) for t in range(self.rcpsp.T)), name = 'resources')
        self.bfc.setParam('OutputFlag', False)
        self.bfc.optimize()
        bfc_pd = -1
        for t in range(self.rcpsp.T):
            if round(self.bfc._v_ti[(t, self.rcpsp.n_activities-1)].x) == 1:
                bfc_pd = t
                break
        return bfc_pd
#problem_name = "datasets/patterson/pat8.rcp"
#problem_name = "datasets/RG5/Pat20.rcp"
#print("xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx")
for i in range(55, 158):
    problem_name = "datasets/RG10/rg10_06/Pat%d.rcp" % i 
    problem_name = "datasets/RG6/Pat%d.rcp" % i
    print("xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx")
    #input()
    alfa = 25
    beta = 100
    rcpsp = RCPSP(problem_name, 
                  alfa = alfa, 
                  resource_cost = beta, 
                  approach_rub = "sum", 
                  zero_rk = False,
                  one_resource_only = False)
    lc = LC_alg(rcpsp, tricks_used = "")
    lc.solve()
    bl_sol_makespan = lc.solution_start_times[rcpsp.n_activities-1]
    bl_sol_objval = round(lc.hpm.objVal)
    bl_sol_resources = sum(lc.solution_resource_hired)
    bl_sol_jumps = 0
    for t in range(1, rcpsp.T):
        for k in range(rcpsp.n_resources):
            bl_sol_jumps += lc.solution_j_tk_p[t, k]
            bl_sol_jumps += lc.solution_j_tk_m[t, k]
    integrated = IntegratedApproach(rcpsp)
    gamma_bounds = [-10000, 10000]
    iter_n = 0
    print(problem_name)
    while iter_n <= 100:
        gamma = (gamma_bounds[1]+gamma_bounds[0]) / 2
        int_sol_makespan, int_sol_resources, int_sol_jumps = integrated.solve(alfa, beta, gamma)
        int_objval = alfa * int_sol_jumps + beta * int_sol_resources
        print("Gamma %.2f, int_makespan %d (bl %d), int_resources %d (bl %d), int_jumps %d (bl %d), int_objval %d (bl %d), int_objval: %.2f, Gamma bounds: %s"\
              % (gamma, int_sol_makespan, bl_sol_makespan, int_sol_resources, bl_sol_resources, int_sol_jumps, bl_sol_jumps,
                int_objval,  bl_sol_objval, integrated.hpm.objVal, str(gamma_bounds)))
        if int_sol_makespan > bl_sol_makespan:
            gamma_bounds[0] = gamma
        elif int_sol_makespan < bl_sol_makespan:
            gamma_bounds[1] = gamma
        else:
            #print("Mario was right on this one :( ")
            break
        if abs(((gamma_bounds[1]+gamma_bounds[0]) / 2) - gamma) < 0.01:
            print("I won!!!")
            sys.exit(0)
            break
    if int_objval != bl_sol_objval:
        print("Heureka!")
        input()