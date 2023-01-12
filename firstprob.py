def f(x):
    ac = 0
    cc = 0
    gc = 0
    tc = 0
    for elem in x:
        if (elem=="A"):
            ac+=1
        elif (elem=="C"):
            cc+=1
        elif (elem == "G"):
            gc+=1
        else:
            tc+=1
    print("{0} {1} {2} {3}".format(ac,cc,gc,tc))
    
def transcribe(str_dna):
    str_rna = []
    for elem in str_dna:
        if (elem == "T"):
            str_rna.append("U")
        else:
            str_rna.append(elem)
    print(''.join(str_rna))
            
    
def f_2(x_l):
    neu = list(x_l)
    neu.reverse()
    fin_list = []
    for elem in neu:
        if (elem=="A"):
            fin_list.append("T")
        elif (elem=="C"):
            fin_list.append("G")
        elif (elem == "G"):
            fin_list.append("C")
        else:
            fin_list.append("A")
    print(''.join(fin_list))
    

import itertools

def mend1_wrapper(a,b,k,m,n):
    tot = k + m + n
    if (a == b):
        if (a == 'm'):
            return ((m-1)/(tot-1))*(m/tot)*(3/4)
        elif (a == 'n'):
            return 0
        else:
            return (k/tot)*((k-1)/(tot-1))
    else:
        if (a == 'k'):
            if (b == 'm'):
                return (k/tot)*(m/(tot-1))
            else:
                return (k/tot)*(n/(tot-1))
    
        elif (a == 'm'):
            if (b == 'k'):
                return (m/tot)*(k/(tot-1))
            else:
                return (m/tot)*(n/(tot-1))*(1/2)
        else:
            if (b == 'k'):
                return (n/tot)*(k/(tot-1))
            else:
                return (n/tot)*(m/(tot-1))*(1/2)

def mend1(k,m,n):
    dict1 = {'k':k, 'm':m, 'n':n}
    extracomb = [('k','k'), ('m','m'), ('n','n')]
    anslst = list(itertools.permutations(dict1.keys(),2))
    anslst.extend(extracomb)
    print(anslst)
    sumlst = []
    for elem in anslst:
        sumlst.append(mend1_wrapper(elem[0],elem[1],k,m,n))
    print(sumlst)
    return sum(sumlst)

def exp_off_prac(n):
    ans=0
    k=0
    while (k<n):
        ans += ((k+1)*(1/n))
        print(ans)
        k += 1
    return ans

def exp_off(a,b,c,d,e,f):
    p2 = [a,b,c,d,e,f]
    ans=0
    q=0
    while (q < 6):
        if (q < 3):
            ans += (2 * p2[q])
            q+=1
        elif (q == 3):
            ans += (2 * .75 * p2[q])
            q+=1
        elif (q == 4):
            ans += (2 * .5 * p2[q])
            q+=1
        else:
            q+=1  
    return ans
        
"""
prac
"""
def factorial(x):
    if x==1:
        print(x)
        return 1
    else:
        print(x)
        return(x*factorial(x-1))
"""
pend
"""
def fib(n):
    if (n==1):
        return 1
    elif (n==2):
        return 1
    else:
        return fib(n-1)+fib(n-2)
    
def rec_rel(n_1,k_1):
    if (n_1==1):
        return 1
    elif (n_1==2):
        return 1
    else:
        return rec_rel(n_1-1,k_1)+ (rec_rel(n_1-2,k_1)*(k_1))  

def GCperc(dna):
    answer = len(list(filter(lambda x: x == 'G', list(dna))))+len(list(filter(lambda x: x == 'C', list(dna))))
    return (answer / len(dna))*100      
        

def str_clean(s):
    splt_s=s.split("\n")
    splt_s.remove("")
    maindict = {}
    k=0
    while (k<len(splt_s)):
        if (list(splt_s[k])[0]==">"):
            curr_val = ""
            q = k+1
            while (q < len(splt_s)):
                if (list(splt_s[q])[0]==">"):
                    break
                else:
                    curr_val+=splt_s[q]
                    q+=1
            maindict[splt_s[k]]=curr_val
            k+=1
        else:
            k+=1
    perdict = {}
    for elem in maindict:
        perdict[elem] = GCperc(maindict[elem])
    print(maindict)
    print(perdict)
    print((list(perdict.keys())[(list(perdict.values()).index(max(list(perdict.values()))))])+"\n"+ str(max(list(perdict.values()))))
        
def translate(strand):
    rna_dict = {"UUU":'F', "UUC":'F',"UUA":'L',"UUG":'L',"UCU":'S',
                "UCC":'S', "UCA":'S',"UCG":'S',"UAU":'Y',"UAC":'Y',
                "UAA":"Stop", "UAG":"Stop", "UGU":'C', "UGC":'C',
                "UGA":"Stop", "UGG":'W', "CUU": 'L', "CUC": 'L', "CUA": 'L',
                "CUG":'L', "CCU":'P', "CCC": 'P', "CCA":'P', "CCG":'P',
                "CAU": 'H', "CAC":'H', "CAA": 'Q', "CAG": 'Q', "CGU":'R',
                "CGC":'R', "CGA":'R', "CGG":'R', "AUU":'I', "AUC":'I', 
                "AUA":'I', "AUG":'M', "ACU": 'T', "ACC":'T', "ACA":'T',
                "ACG":'T', "AAU": 'N', "AAC":'N', "AAA":'K', "AAG" :'K',
                "AGU":'S',"AGC":'S',"AGA":'R',"AGG":'R',"GUU":'V',"GUC":'V',
                "GUA":'V', "GUG":'V', "GCU":'A',"GCC":'A', "GCA":'A',
                "GCG":'A',"GAU":'D',"GAC":'D',"GAA":'E',"GAG":'E',
                "GGU":'G',"GGC":'G',"GGA":'G',"GGG":'G'}
    k=0
    cod_list = []
    while (k<(len(strand)-2)):
        if (k == (len(strand)-3)):
            cod_list.append(''.join(list(strand)[k:len(strand)]))
            k+=1
        else:
            cod_list.append(''.join(list(strand)[k:k+3]))
            k+=3
    polypep = ""
    for elem in cod_list:
        if (rna_dict[elem] == "Stop"):
            break
        else:
            polypep += rna_dict[elem]
    return polypep
        
            
def Hammingdist(a,b):
    ind = 0
    count = 0
    while (ind < len(a)):
        if ((list(a)[ind]) == (list(b)[ind])):
            ind+=1
        else:
            count+=1
            ind+=1
    return count

def motif_find(s,t):
    ind_1=0
    pos_str=""
    while (ind_1 < len(s)):
        if (ind_1 == (len(s)-len(t))):
            if (s[ind_1:]==t):
                pos_str += str(ind_1+1)
                ind_1+=1
            else:
                ind_1+=1
        elif (s[ind_1:len(t)+ind_1] == t):
            pos_str += (str(ind_1+1) + " ")
            ind_1+=1
        else:
            ind_1+=1
    if (pos_str[(len(pos_str)-1):]==" "):
        pos_str = pos_str[:len(pos_str)-1]
    return pos_str
            
            
def mort_fib(n,m):
    if (n==1):
        return 1
    elif (n==2):
        return 1
    elif (n==3):
        return mort_fib(n-1,m)+mort_fib(n-2,m)    
    else:
        return (mort_fib(n-1,m)+mort_fib(n-2,m))-(2**(n-(m-1)))
                                                        
def mflist(n,m):
    newlst = []
    for elem in range(7):
        newlst.append(mort_fib(elem+1,m))
    print(newlst)
    
def prot_mass(s):
    pmt = {"A" : 71.03711,
           "C" : 103.00919,
           "D" : 115.02694,
           "E" : 129.04259,
           "F" : 147.06841,
           "G" : 57.02146,
           "H" : 137.05891,
           "I" : 113.08406,
           "K" : 128.09496,
           "L" : 113.08406,
           "M" : 131.04049,
           "N" : 114.04293,
           "P" : 97.05276,
           "Q" : 128.05858,
           "R" : 156.10111,
           "S" : 87.03203,
           "T" : 101.04768,
           "V" : 99.06841,
           "W" : 186.07931,
           "Y" : 163.06333}
    
    ans = round(sum([list(pmt.values())[list(pmt.keys()).index(elem)] for elem in s]),3)
    return ans
        
def shared_motif(s):
    splt_s=s.split("\n")
    maindict = {}
    k=0
    while (k<len(splt_s)):
        if (list(splt_s[k])[0]==">"):
            curr_val = ""
            q = k+1
            while (q < len(splt_s)):
                if (list(splt_s[q])[0]==">"):
                    break
                else:
                    curr_val+=splt_s[q]
                    q+=1
            maindict[splt_s[k]]=curr_val
            k+=1
        else:
            k+=1
    perdict = {}
    for elem in maindict:
        perdict[elem] = GCperc(maindict[elem])
    print(maindict)
    print((list(perdict.keys())[(list(perdict.values()).index(max(list(perdict.values()))))])+"\n"+ str(max(list(perdict.values()))))
    
def subsetgen(s):
    ind1 = 0
    mainlist = []
    while (ind1 < len(s)):
        ind2 = 0
        curr_list = []
        while (ind2 <= ind1):
            curr_list.append(list(s)[ind2])
            ind2+=1
        currsubstr = ''.join(curr_list)
        mainlist.append(currsubstr)
        ind1+=1
    return mainlist 

import random
        
def bpgen(n):
    newlst = []
    ind1=0
    while (ind1<n):
        x = random.randint(ind1,n)
        if (x%4==0):
            newlst.append('A')
            ind1+=1
        elif (x%4==1):
            newlst.append('C')
            ind1+=1
        elif (x%4==2):
            newlst.append('G')
            ind1+=1
        else:
            newlst.append('T')
            ind1+=1
    ansstr = ''.join(newlst)
    return ansstr 

def checker(a,alos):
    ans = True
    for elem in alos:
        if (a in elem):
            continue
        else: 
            ans = False
    return ans

def sharmo(n,alos):
    alon = [elem+1 for elem in range(n)]
    startind = 0
    endind = 0
    finalans = []
    while (endind < len(alon)):
        if (startind == endind):
            curr_ans = ""
            finalans.append(curr_ans)
        else:
            cll = alon[startind:endind]
            curr_ans = "".join(map(str,cll))
            finalans.append(curr_ans)
        if (checker(curr_ans, alos)):
            endind+=1
            finalans.append(curr_ans)
        else:
            startind += 1 
            endind += 1
            alon = alon[startind:]
            cll = alon[startind:endind]
            curr_ans = "".join(map(str,cll))
    return finalans
            
        