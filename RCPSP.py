# -*- coding: utf-8 -*-
"""
Created on Sat Oct 13 10:52:02 2018

@author: pmilicka
"""
import numpy as np
import gurobipy as g
from matplotlib import pyplot as plt
from matplotlib.patches import Rectangle
import time

class RCPSP:
    def __init__(self, filename, alfa = 1, resource_cost = 1, input_format = 'patterson', approach_rub = None, one_resource_only = False, zero_rk = False):
        start_time = time.time()
        self.one_resource_only = one_resource_only
        self.zero_rk = zero_rk
        if input_format == 'patterson':
            self.readPatterson(filename)
        else:
            raise Exception('Unknown format to be read')
        self.alfa = alfa    
        
        self.resource_cost = [resource_cost] * self.n_resources
        self.computeEarliestStartTimes()
        self.computeMaxDuration()
        self.computeShortestMakespan()
        self.computeResourceUB(approach_rub)
        self.init_time = time.time() - start_time
        
        
    def readPatterson(self, filename):
        self.activity_duration = []
        self.activity_successor = []
        self.activity_resource = []
        self.n_edges = 0
        act_idx = 0
        with open(filename) as fin:
            tmp = fin.readline().split()
            if len(tmp) != 2:
                tmp = fin.readline().split()
            [self.n_activities, self.n_resources] = list(map(int, tmp))
            self.resource_available = list(map(int, fin.readline().split()))
            for line in fin:
                tmp = line.split()
                if len(tmp) < 2:
                    continue
                line = list(map(int, tmp))
                self.activity_duration.append(line[0])
                self.activity_resource.append(line[1:1+self.n_resources])
                self.activity_successor.append([x-1 for x in line[2+self.n_resources:]]) 
                act_idx += 1
                self.n_edges += len(self.activity_successor[-1])
            if self.one_resource_only:
                self.activity_resource = [[sum(lst)] for lst in self.activity_resource]
                self.resource_available = [sum(self.resource_available)]
                self.n_resources = 1
            if self.zero_rk:
                for k in range(self.n_resources):
                    self.resource_available[k] = 0
        
    def computeShortestMakespan(self):
        pm = g.Model("Shortest Makespan")
        pm._v_ti = pm.addVars(self.T, self.n_activities, vtype = g.GRB.BINARY, name = 'v_ti') #
        # Set the objective
        pm.setObjective(g.quicksum([t * pm._v_ti[(t, self.n_activities-1)] for t in range(self.T)]), g.GRB.MINIMIZE)
        # Set the constraints 
        # 1. Every activity must start
        pm.addConstrs(g.quicksum([pm._v_ti[(t, i)] for t in range(self.T)]) == 1 for i in range(self.n_activities))
        # 2. Precedence relations
        pm.addConstrs((g.quicksum([t * pm._v_ti[(t, j)] for t in range(self.T)]) \
                     - g.quicksum([t* pm._v_ti[(t, i)] for t in range(self.T)]) \
                     >= self.activity_duration[i] \
                     for i in range(self.n_activities) for j in self.activity_successor[i]), name = 'precedence')
        # (Optional): set start times before the earliest start time of an activity to 0
        pm.addConstrs(pm._v_ti[(t, i)] == 0 for i in range(self.n_activities) for t in range(0, self.activity_est[i]))
        # Get the smallest possible project makespan
        pm.optimize()
        for t in range(self.T):
            if pm._v_ti[(t, self.n_activities-1)].x > 0.9:
                self.shortest_makespan = t
                break
        
    def computeResourceUB(self, approach = None):
        self.resource_ub = []
        if approach is None:
            self.resource_ub = [250] * 4 # corresponds approx with 8 bits
            return
        elif approach == "sum":
            tmp = np.array(self.activity_resource)
            for k in range(self.n_resources):
                self.resource_ub.append(np.sum(tmp[:, k]))
        elif approach == "rap":
            rap = g.Model('Resource Availability Problem')
            rap._v_ti = rap.addVars(self.T, self.n_activities, vtype = g.GRB.BINARY, name = 'v_ti') #
            # Set the objective
            rap.setObjective(g.quicksum([t * rap._v_ti[(t, self.n_activities-1)] for t in range(self.T)]), g.GRB.MINIMIZE)
            # Set the constraints 
            # 1. Every activity must start
            rap.addConstrs(g.quicksum([rap._v_ti[(t, i)] for t in range(self.T)]) == 1 for i in range(self.n_activities))
            # 2. Precedence relations
            rap.addConstrs((g.quicksum([t * rap._v_ti[(t, j)] for t in range(self.T)]) \
                         - g.quicksum([t* rap._v_ti[(t, i)] for t in range(self.T)]) \
                         >= self.activity_duration[i] \
                         for i in range(self.n_activities) for j in self.activity_successor[i]), name = 'precedence')
            # (Optional): set start times before the earliest start time of an activity to 0
            rap.addConstrs(rap._v_ti[(t, i)] == 0 for i in range(self.n_activities) for t in range(0, self.activity_est[i]))
            # Get the upper bounds for resources by solving resource availibility model
            y_ka = rap.addVars(self.n_resources, vtype = g.GRB.CONTINUOUS, lb = 0, name="y_ka")
            rap.setObjective(g.quicksum([y_ka[k]*self.resource_cost[k] for k in range(self.n_resources)]), g.GRB.MINIMIZE)
            rap.addConstrs((g.quicksum([self.activity_resource[i][k] * rap._v_ti[(t_, i)]
                            for i in range(self.n_activities)
                            for t_ in range(max(0, t - self.activity_duration[i] + 1), t + 1) ])
                       <= self.resource_available[k] + y_ka[k]
                       for t in range(self.T) for k in range(self.n_resources)), name = 'res_req')
            rap.addConstr(rap._v_ti[(self.shortest_makespan, self.n_activities-1)] == 1)
            rap.optimize()
            self.resource_ub = []
            for k in range(self.n_resources):
                self.resource_ub.append(round(y_ka[k].x))
        self.resource_nbits = []
        for k in range(self.n_resources):
            self.resource_nbits.append(len(format(self.resource_ub[k], "b")))
            
            
            
        
    def computeMaxDuration(self):
        T = sum(self.activity_duration) + 1
        fm = g.Model('tmp')
        v_ti = fm.addVars(T, self.n_activities, vtype = g.GRB.BINARY, name = 'v_ti') #
        # Set the objective
        fm.setObjective(g.quicksum([t * v_ti[(t, self.n_activities-1)] for t in range(T)]), g.GRB.MINIMIZE)
        # Set the constraints 
        # 1. Every activity must start
        fm.addConstrs(g.quicksum([v_ti[(t, i)] for t in range(T)]) == 1 for i in range(self.n_activities))
        # 2. Precedence relations
        fm.addConstrs((g.quicksum([t * v_ti[(t, j)] for t in range(T)]) \
                     - g.quicksum([t* v_ti[(t, i)] for t in range(T)]) \
                     >= self.activity_duration[i] \
                     for i in range(self.n_activities) for j in self.activity_successor[i]), name = 'precedence')
        # (Optional): set start times before the earliest start time of an activity to 0
        fm.addConstrs(v_ti[(t, i)] == 0 for i in range(self.n_activities) for t in range(0, self.activity_est[i]))
        # Get the upper bounds for resources by solving resource availibility model
        resource_min = []
        if self.zero_rk:
            resource_min = np.max(np.array(self.activity_resource), axis = 0)
        else:
            resource_min = self.resource_available
        fm.addConstrs((g.quicksum([self.activity_resource[i][k] * v_ti[(t_, i)]
                            for i in range(self.n_activities)
                            for t_ in range(max(0, t - self.activity_duration[i] + 1), t + 1) ])
                       <= resource_min[k]
                       for t in range(T) for k in range(self.n_resources)), name = 'res_req')
        # Get the shortest makespan with no additional resources hired
        fm.optimize()
        for t in range(T):
            if round(v_ti[(t, self.n_activities-1)].x) == 1:
                self.T = t + 1
        #self.T = sum(self.activity_duration) + 1

    def computeEarliestStartTimes(self):
        self.activity_est = [float('-inf')] * self.n_activities
        self.activity_est[0] = 0
        for i in range(self.n_activities):
            for j in range(self.n_activities):
                for succ in self.activity_successor[i]:
                    if self.activity_duration[i] + self.activity_est[i] > self.activity_est[succ]:
                        self.activity_est[succ] = self.activity_duration[i] + self.activity_est[i]
    
    def plotSchedule(self, start_times):
        plt.figure(figsize = (15, 6))
        ax = plt.gca()
        plt.xlim([0, start_times[self.n_activities-1]])
        plt.ylim([0, self.n_activities])
        for i in range(self.n_activities):
            ax.add_patch(Rectangle((start_times[i], i - 0.5), height =1, width = self.activity_duration[i], linewidth = 1, fill = True, zorder = 4, linestyle = "-", edgecolor = "black"))
            if i > 0 and i < self.n_activities - 1:
                plt.text(0.5*(start_times[i]+start_times[i]+self.activity_duration[i]), i, str(i) + " : " + str(self.activity_resource[i]), horizontalalignment='center',verticalalignment='center',fontsize=10, color='black',  zorder = 6, fontweight = 'bold')
        plt.grid(axis = 'xy', which = 'both')
        plt.xticks(list(range(start_times[self.n_activities-1] +1)))
        plt.yticks(list(range(self.n_activities)))
        plt.title("Project Schedule")
        plt.grid()
        plt.show()
    
    def plotResourceUtilization(self, start_times, hired_res):
        res_demand = np.zeros((self.n_resources, self.T), dtype = int)
        for i in range(self.n_activities):
            for k in range(self.n_resources):
                for t in range(start_times[i], start_times[i] + self.activity_duration[i]):
                    res_demand[k, t] += self.activity_resource[i][k]
        for k in range(self.n_resources):
            plt.figure(figsize = (20, 6))
#            res_req = [0] * (start_times[self.n_activities-1] + 1)
#            for i in range(self.n_activities):
#                for t in range(self.activity_duration[i]):
#                    res_req[t + start_times[i]] += self.activity_resource[i][k]
            for t in range(self.T):
                plt.gca().add_patch(Rectangle((t, 0), height = res_demand[k, t], width = 1, alpha = 1,zorder = 2))
            #plt.plot(list(range(0,start_times[self.n_activities-1] + 1)), res_req, 'b')
            plt.plot([0, start_times[self.n_activities-1]], [self.resource_available[k] + hired_res[k]] * 2, 'r')
            plt.xlabel('Time')
            plt.xlim([0, start_times[self.n_activities-1]])
            plt.ylim([0, self.resource_available[k] + hired_res[k] + 0.3])
            plt.title('Resource type %d use' % k)
            plt.ylabel("Resource units used")
            plt.xticks(list(range(start_times[self.n_activities-1] +1)))
            plt.yticks(list(range(0,self.resource_available[k] + hired_res[k] + 1)))
            plt.grid(axis = 'both')
            plt.show()
    
if __name__ == "__main__":
    pass    
    

    
    
#self = self('test_data.rcp')
#bb = BranchAndBound(self)
#bb.root.solveHPS()