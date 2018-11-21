# -*- coding: utf-8 -*-
"""
Created on Sun Nov  4 09:31:47 2018

@author: pmilicka
"""

class BILP:
    def __init__(self, RCPSP):
        self.rcpsp = RCPSP

    def solve(self):
        bb = BranchAndBound(self, self.rcpsp)
        
        
class BranchAndBound:

    def __init__(self, rcpsp):
        self.root = Node(k = 0, h_lead = [float('inf')] * rcpsp.n_resources, h_foll = [], s_lead = set(), s_foll = set(), rcpsp = rcpsp)
        self.curr_node = self.root
        self.rcpsp = rcpsp
        f_ub = float('inf')
    
class Node:
    alfa = 1
    
    def __init__(self, k, h_lead, h_foll, s_lead, s_foll, rcpsp):
        self.k = k
        self.h_lead = h_lead
        self.h_foll = h_foll
        self.s_lead = s_lead
        self.s_foll = s_foll
        self.rcpsp = rcpsp
        self.m = g.Model('')
        
        
    def solveContinousBLP(self):
        pass
        
    def solveHPS(self):
        m = g.Model('')
        v = m.addVars(self.rcpsp.T, self.rcpsp.n_activities, vtype = g.GRB.BINARY, name = 'v')
        y = m.addVars(self.rcpsp.n_resources, vtype = g.GRB.INTEGER, lb = 0)
        m.setObjective(self.rcpsp.alfa * g.quicksum([t * v[(t, self.rcpsp.n_activities - 1)] for t in range(self.rcpsp.T)]) \
                       + g.quicksum([self.rcpsp.resource_cost[k] * y[(k)] for k in range(self.rcpsp.n_resources)]), g.GRB.MINIMIZE)
#        m.addConstrs((v[t,i] == 0 for t in range(self.rcpsp.T)
#                            for i in range(self.rcpsp.n_activities)))
        
        m.addConstrs(g.quicksum([v[(t, i)] for t in range(self.rcpsp.T)]) == 1 for i in range(self.rcpsp.n_activities)) #all activities have to start
        m.addConstrs(g.quicksum([t * v[(t, j)] for t in range(self.rcpsp.T)]) \
                     - g.quicksum([t* v[(t, i)] for t in range(self.rcpsp.T)]) \
                     >= self.rcpsp.activity_duration[i] \
                     for i in range(self.rcpsp.n_activities) for j in self.rcpsp.activity_successor[i])
        m.addConstrs(v[(t, i)] == 0 for i in range(self.rcpsp.n_activities) for t in range(0, self.rcpsp.activity_est[i]))
        for t in range(self.rcpsp.T):
            for k in range(self.rcpsp.n_resources):
                m.addConstr(g.quicksum([self.rcpsp.activity_resource[i][k] * v[(t_, i)] 
                                for i in range(self.rcpsp.n_activities) 
                                for t_ in range(max(0, t - self.rcpsp.activity_duration[i]), t)])
                            <= self.rcpsp.resource_available[k] + y[k], name = 't=%d,k=%d' %(t, k)) #TODO considered more skilled workers)
        
        m.relax()    
        m.optimize()
        m.write("file.lp")
        print("obj val: %.2f", m.objVal)
        start_times = {}
        hired_res = []
        for (t, i) in v:
            if abs(v[(t, i)].x - round(v[(t,i)].x)) <= m.Params.IntFeasTol \
                and v[(t,i)].x > 0.9:
                print("%s val: %.2f " %(v[(t,i)].varName, v[(t,i)].x))
                start_times[i] = t
        for k in range(self.rcpsp.n_resources):
            print("y_%d %.1f" % (k, y[k].x))
            hired_res.append(round(y[k].x))
#        print(start_times)
#        self.rcpsp.plotSchedule(start_times, hired_res)