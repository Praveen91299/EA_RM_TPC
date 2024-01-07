#QRM TPC rates
#Priya Nadkarni and Praveen Jayakumar

import math
import matplotlib.pyplot as plt

def binom_sum(m,start,end):
    '''
    Binomial sum (m C start) + (m C start+1) + ... + (m C end)
    '''
    sum_binom = 0
    for i in range(end-start+1):
        sum_binom += math.comb(m,start+i)
    return sum_binom
 
def catalytic_rate(r, m): #TPC
    '''
    Catalytic rate of a EA RM TPC constructed from [RM(r, m), RM(r, m)]
    '''
    m1=m
    m2=m 
    r1 = r
    r2 = r
    return (pow(2, m1+m2) - (2*binom_sum(m1,0,m1-r1-1)*binom_sum(m2,0,m2-r2-1)))/pow(2, m1+m2)

def EA_rate(r, m):
    '''
    Code rate for EA RM TPC constructed from [RM(r, m), RM(r, m)]
    '''
    m1 =m
    m2 =m
    r1 =r
    r2 =r
    return (pow(2, m1+m2) - (2*binom_sum(m1,0,m1-r1-1)*binom_sum(m2,0,m2-r2-1)) + (binom_sum(m1,r1+1,m1-r1-1)*binom_sum(m2,r2+1,m2-r2-1)))/pow(2, m1+m2)

def EA_rate_gen(r1, r2, m1, m2):
    '''
    Code rate for EA RM TPC constructed from [RM(r1, m1), RM(r2, m2)]
    '''
    return (pow(2, m1+m2) - (2*binom_sum(m1,0,m1-r1-1)*binom_sum(m2,0,m2-r2-1)) + (binom_sum(m1,r1+1,m1-r1-1)*binom_sum(m2,r2+1,m2-r2-1)))/pow(2, m1+m2)

def catalytic_rate_gen(r1, r2, m1, m2):
    '''
    Catalytic rate for EA RM TPC constructed from [RM(r1, m1), RM(r2, m2)]
    
    '''
    return (pow(2, m1+m2) - (2*binom_sum(m1,0,m1-r1-1)*binom_sum(m2,0,m2-r2-1)))/pow(2, m1+m2)

# theorem 3 bound
def lr(r):
    '''
    lr, as defined in theorem 3
    '''
    diff = 1
    i = 0
    prods = [math.comb(2*r, u) for u in range(r+1)]
    while diff > 0:
        i+=1
        prods = [p*(2*r + i)/(2*r + i -u) for u, p in enumerate(prods)]
        diff = sum(prods) - ((2**(2*r + i))/(2 + math.sqrt(2)))
    return i-1

# obtaining lr (Table 1)
r_list = [1, 5, 6, 13, 14, 24, 25, 39, 40, 57, 58, 78, 79, 103, 104, 131, 132, 162]
print('\nTable 1 bounds:\n')
for r in r_list:
    print('r = {}, l(r) = {}'.format(r, lr(r)))

#obtaining catalytic rates (Table 2)
code_params = [(1, 4), (2, 6), (3, 8), (4, 10), (5, 12)]
print('\nTable 2 Catalytic rates:\n')
for params in code_params:
    r, m = params[0], params[1]
    print('(r, m) = ({}, {}), Catalytic rate: {}'.format(r, m, catalytic_rate(r, m)))


# obtaining Figure 2 (m-r plane)
points = []
rs = []
ms = []
m_list = list(range(1, 100))
for m in m_list:
    for r in range(int((m-1)/2)+1):
        if m-r-1 < m/2: #checking not self-dual containing.
            continue
        cr = catalytic_rate(r, m)
        if cr > 0:
            #print(r, m, cr)
            points.append([r, m])
            rs.append(r)
            ms.append(m)

rs_thm = []
ms_thm = []
points_thm = []

for r in range(max(rs)+1):
    l =  lr(r)
    m_allowed = [2*r + s for s in range(1, l+1)]
    for m in m_allowed:
        if m-r-1 < m/2: #checking not self-dual containing.
            continue
        cr = catalytic_rate(r, m) #sanity check
        if cr > 0:
            #print(r, m, cr)
            points_thm.append([r, m])
            rs_thm.append(r)
            ms_thm.append(m)

r_lemma = list(range(max(rs)+1))
m_lemma = [r*2 + 2 for r in r_lemma]

r_lemma_n = list(range(max(rs)+1))
m_lemma_n = [r*3 + 1 for r in r_lemma]

plt.scatter(rs_thm, ms_thm, s = 15, label = 'Theorem 3')
plt.scatter(rs, ms, s=2, label = 'found by brute force')
plt.plot(r_lemma, m_lemma, label = 'lemma 6, m = 2r+2', color='r')
#plt.plot(r_lemma_n, m_lemma_n, label = 'm = 3r+1', color='g')
plt.legend()
plt.xlabel('r')
plt.ylabel('m')
plt.savefig('EATPRM.png')