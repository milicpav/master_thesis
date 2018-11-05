# -*- coding: utf-8 -*-
"""
Created on Sat Oct 13 10:52:02 2018

@author: pmilicka
"""
import numpy as np
import gurobipy as g
from matplotlib import pyplot as plt
from matplotlib.patches import Rectangle

class RCPSP:
    def __init__(self, filename, alfa = 1, resource_cost = 1, input_format = 'patterson'):
        if input_format == 'patterson':
            self.readPatterson(filename)
        else:
            raise Exception('Unknown format to be read')
        self.alfa = alfa    
        self.resource_cost = [resource_cost] * self.n_resources
        self.computeEarliestStartTimes()
                        
    def readPatterson(self, filename):
        self.activity_duration = []
        self.activity_successor = []
        self.activity_resource = []
        self.n_edges = 0
        act_idx = 0
        with open(filename) as fin:
            [self.n_activities, self.n_resources] = list(map(int,fin.readline().split()))
            self.resource_available = list(map(int, fin.readline().split()))
            for line in fin:
                line = list(map(int, line.split()))
                self.activity_duration.append(line[0])
                self.activity_resource.append(line[1:1+self.n_resources])
                self.activity_successor.append([x-1 for x in line[2+self.n_resources:]]) 
                act_idx += 1
                self.n_edges += len(self.activity_successor[-1])
        self.T = sum(self.activity_duration) + 1

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
        
        for i in range(self.n_activities):
            #plt.plot([start_times[i], start_times[i] + self.activity_duration[i]], [i, i])
            #Rectangle((start_times[i], i - 0.5), height =1, width = self.activity_duration[i], linewidth = 1, fill = True, zorder = 4)
            ax.add_patch(Rectangle((start_times[i], i - 0.5), height =1, width = self.activity_duration[i], linewidth = 1, fill = True, zorder = 4, linestyle = "-", edgecolor = "black"))
            if i > 0 and i < self.n_activities - 1:
                plt.text(0.5*(start_times[i]+start_times[i]+self.activity_duration[i]), i, str(i) + " : " + str(self.activity_resource[i]), horizontalalignment='center',verticalalignment='center',fontsize=10, color='black',  zorder = 6, fontweight = 'bold')
        plt.grid(axis = 'y', which = 'both')
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
    

    
    
#rcpsp = RCPSP('test_data.rcp')
#bb = BranchAndBound(rcpsp)
#bb.root.solveHPS()