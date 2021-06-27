## This code is a Grid Search generator for experimental settings for the Switching Bandits.
import numpy as np
import io
from matplotlib.ticker import ScalarFormatter, EngFormatter, FuncFormatter
import matplotlib.pyplot as plt

#Utility functions for the generation of the strings---

def createEnvironment(cp, t='bernoulli'):
    '''This function generates string in the experiment format in order to create the environment
        "cp" is a dictionary of the type {arm:[(mean, until),...]}'''
    env_string = ""
    
    for arm in cp.keys():
        if t=='bernoulli':
            env_string += "BernoulliNonStationary " + str(len(cp[arm])) + " "
        else:
            env_string += t + " " + str(len(cp[arm])) + " "

        for el in cp[arm]:
            env_string += str(el[0]) + ' ' + str(el[1]) + ' '
        
        env_string += '\n'
    
    return env_string

def printEnvironment(dictEnv, T):
    for arm in dictEnv.keys():
        y = np.zeros(T)
        prev = 0
        for piece in dictEnv[arm]:
            y[prev:piece[1]] += [piece[0] for i in range(prev, piece[1])]
            prev = piece[1]
        plt.plot(y)
    plt.show()

def get_timestep(dictEnv):
    maxtime = 0
    for arm,cps in dictEnv.items():
        maxtime = max([t[1] for t in cps])
    return maxtime

def concatEnvironments(env1, env2):
    ts = get_timestep(env1)
    newEnv = {k: env1[k]+list(map(lambda x: (x[0],x[1]+ts), env2[k])) for k in env1.keys()}
    return newEnv

def computeMaxCohen(dictEnv):
    min_delta = np.inf
    for arm in dictEnv.keys():
        for i in range(len(dictEnv[arm])):
            for j in range(i, len(dictEnv[arm])):
                if dictEnv[arm][i][0] != dictEnv[arm][j][0]:
                    min_delta = min(min_delta, abs(dictEnv[arm][i][0] - dictEnv[arm][j][0]))
    return min_delta/2

def maxStationaryPeriod(dictEnv):
    max_period = 0
    for arm, cps in dictEnv.items():
        l = [0]+[x[1] for x in cps]
        max_diff = max([l[i+1]-l[i] for i in range(len(l)-1)])
        max_period = max(max_period,max_diff)
    return max_period

def minArmsDifference(dictEnv):
    min_difference = 1.0
    stationary_periods = sorted(set([x[1] for a,cp in dictEnv.items() for x in cp]))
    print(stationary_periods)
    for cp in stationary_periods:
        means = []
        for arm,cps in dictEnv.items():
            for x in cps:
                if x[1] <= cp:
                    means.append(x[0])
                    break
        diffs = [abs(means[i]-means[j]) for i in range(len(means)) for j in range(i+1,len(means))]
        min_difference = min(min_difference, min(diffs))
    return min_difference

#Experimental settings
exp4 = {'a1':[(0.5,20000),(0.9,40000),(0.5,60000),(0.9,80000),(0.5,100000)],
        'a2':[(0.8,20000),(0.6,40000),(0.8,60000),(0.6,80000),(0.8,100000)]} #Toy example
        
exp12_var1 = {'a1':[(0.60,35000),(0.40,65000),(0.60,100000)],
              'a2':[(0.30,35000),(0.70,65000),(0.30,100000)]}

exp12_var2 = {'a1':[(0.60,35000),(0.40,65000),(0.60,100000)],
              'a2':[(0.40,35000),(0.60,65000),(0.40,100000)]} #Simple Setting

exp15 = {'a1':[(0.8,10000),(0.6,20000),(0.8,30000),(0.6,40000),(0.8,50000),(0.6,60000),(0.8,70000),(0.6,80000),(0.8,90000),(0.6,100000)],
        'a2':[(0.7,10000),(0.9,20000),(0.7,30000),(0.9,40000),(0.7,50000),(0.9,60000),(0.7,70000),(0.9,80000),(0.7,90000),(0.9,100000)],
        'a3':[(0.6,10000),(0.4,20000),(0.6,30000),(0.4,40000),(0.6,50000),(0.4,60000),(0.6,70000),(0.4,80000),(0.6,90000),(0.4,100000)],
        'a4':[(0.5,10000),(0.7,20000),(0.5,30000),(0.7,40000),(0.5,50000),(0.7,60000),(0.5,70000),(0.7,80000),(0.5,90000),(0.7,100000)],
        'a5':[(0.4,10000),(0.2,20000),(0.4,30000),(0.2,40000),(0.4,50000),(0.2,60000),(0.4,70000),(0.2,80000),(0.4,90000),(0.2,100000)]} #Complex Setting

#Yahoo envs (use these to run)
yahoo1 = {'a1': [(0.06111736032995876, 5000), (0.06934306569343066, 10000), (0.05570175438596491, 15000), (0.054087277197295634, 20000), (0.05944986690328305, 25000), (0.04927884615384615, 30000), (0.038461538461538464, 35000), (0.05399792315680166, 40000), (0.0617906683480454, 45000), (0.06551141166525781, 50000), (0.05309090909090909, 55000), (0.05614881157423355, 60000), (0.051426426426426426, 65000), (0.05683403068340307, 70000), (0.04964308890330954, 75000), (0.04804327076041998, 80000), (0.04949095022624434, 85000), (0.06199460916442048, 90000)], 
          'a2': [(0.06996974281391831, 5000), (0.05352223534041716, 10000), (0.05785491766800178, 15000), (0.041120381406436236, 20000), (0.03932082216264522, 25000), (0.04250295159386069, 30000), (0.045572916666666664, 35000), (0.04381161007667032, 40000), (0.05044136191677175, 45000), (0.05857203933304831, 50000), (0.054596688974991196, 55000), (0.05890804597701149, 60000), (0.04587495375508694, 65000), (0.050129008477699966, 70000), (0.05357142857142857, 75000), (0.04732172198488334, 80000), (0.05676979846721544, 85000), (0.050666666666666665, 90000)],
          'a3': [(0.04326186830015314, 5000), (0.03628560280920796, 10000), (0.03429844097995546, 15000), (0.04589963280293757, 20000), (0.026006711409395974, 25000), (0.021810250817884406, 30000), (0.032520325203252036, 35000), (0.02967032967032967, 40000), (0.030637254901960783, 45000), (0.031166518254674976, 50000), (0.024991074616208496, 55000), (0.03739316239316239, 60000), (0.030775075987841946, 65000), (0.03259688518652662, 70000), (0.03171390013495277, 75000), (0.03625081645983018, 80000), (0.03066242406711021, 85000), (0.028645833333333332, 90000)], 
          'a4': [(0.026158445440956652, 5000), (0.026558891454965358, 10000), (0.020463847203274217, 15000), (0.023215322112594312, 20000), (0.023478260869565216, 25000), (0.012235817575083427, 30000), (0.016574585635359115, 35000), (0.021052631578947368, 40000), (0.021011162179908074, 45000), (0.020858895705521473, 50000), (0.0226628895184136, 55000), (0.02695134756737837, 60000), (0.02098164106406894, 65000), (0.02148997134670487, 70000), (0.02012072434607646, 75000), (0.024739583333333332, 80000), (0.022132796780684104, 85000), (0.029925187032418952, 90000)], 
          'a5': [(0.03575461046292811, 5000), (0.03257328990228013, 10000), (0.040253279059249206, 15000), (0.0367913148371532, 20000), (0.025984911986588432, 25000), (0.023282887077997673, 30000), (0.009497964721845319, 35000), (0.023835319609967497, 40000), (0.03225806451612903, 45000), (0.034908471690080885, 50000), (0.031487513572204126, 55000), (0.03145917001338688, 60000), (0.03424153166421208, 65000), (0.0340036563071298, 70000), (0.03506450545815415, 75000), (0.033920417482061316, 80000), (0.030826851051570153, 85000), (0.029649595687331536, 90000)]}

filename = "gridSearchYahoo.txt"
path = "exp_config/"
choosen_env = exp15 #Here you can select the environment you want to try
history_version = 1 # 1: logarithmic / 2: infinite
#The line to be put inside every 
regret_type = 0
n_of_arms = len(choosen_env.keys())
n_of_algs = 7
T = get_timestep(choosen_env) #100000

description_line = str(regret_type) + ' ' + str(n_of_arms) + ' ' + str(n_of_algs) + ' ' + str(T)

other_algs = """UCB1 ucb1
D-UCB round 2 d_ucb _d_ucb
SW-UBC round 2 sw_ucb _sw_ucb
AdaptEvE _adapt_eve
CUSUM-UCB _cd_algorithm 0 1 0
BR-CUSUM-UCB(B,/) _cd_algorithm """+str(history_version)+""" 1 0"""

alphas = [0.01, 0.05, 0.1]
ds = [max(0.01,min(0.1,computeMaxCohen(choosen_env)))]

#For the first line of the experimental file
n_exp = len(alphas)*len(ds)
simulations = 200
seed = 1

first_string = str(n_exp) + " " + str(simulations) + " " + str(seed)

env = createEnvironment(choosen_env)
printEnvironment(choosen_env,T)

with open(path+filename,'w') as f:
    f.write(first_string+'\n')
    exp_n = 0
    for a in alphas:
        for d in ds:
            exp_name = 'alpha'+str(a)+'d'+str(d)
            f.write(str(exp_name) + " " + description_line+'\n')
            f.write(env)
            f.write(other_algs+'\n')
            f.write("BR-CUSUM-UCB(B,E) _rec_cd_algorithm "+ str(a) + " " + str(d) +" "+str(history_version)+" 1 0\n")
            exp_n += 1