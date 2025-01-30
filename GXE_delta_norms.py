import os
os.environ["PATH"] += os.pathsep + '/Library/TeX/texbin'
import matplotlib.pyplot as plt # plot stuff
plt.rc('text', usetex=True)




import numpy as np

import scipy as sp
import csv



## This code was used to generate all reaction norm plots used in 
# 'How GXE interactions maintain variation in mutualisms'

def import_Johnson_data():
    # .csvfile is included in supplemental materials. 
    with open('/Users/chriscarlson/Documents/GXE/Data/Johnson_data.csv', newline='') as csvfile:
        data = list(csv.reader(csvfile))


    vals=['A','K','C']
    from itertools import product
    vals = ['A', 'K', 'C']
    keys=[]
    for comb in product(vals, repeat=len(vals)):
        keys.append(comb)

    myco_dict={}

    for key in keys: 
        myco_dict[key]=[[],[]]
        for item in data: 
            if key[0]==item[0] and key[1]==item[1] and key[2]==item[2]:
                myco_dict[key][0].append(float(item[3]))
                myco_dict[key][1].append(float(item[4]))
    
    return myco_dict

def Johnson_plot(myco_dict):
    
    ################### plant parameters: #####################

    # order is: soil - plant - fungus - hyphae - plant growth
    Hm1=np.mean(myco_dict[('K','K','K')][1])
    sehm1=sp.stats.sem(myco_dict[('K','K','K')][1])
    distHm1=myco_dict[('K','K','K')][1]
    print('Hm1',Hm1)

    # order is: soil - plant - fungus - hyphae - plant growth
    Hm2=np.mean(myco_dict[('C','C','C')][1])
    sehm2=sp.stats.sem(myco_dict[('C','C','C')][1])
    distHm2=myco_dict[('C','C','C')][1]
    print('Hm2',Hm2)

    # order is: soil - plant - fungus - hyphae - plant growth
    Hh1=np.mean(myco_dict[('K','K','C')][1])
    sehh1=sp.stats.sem(myco_dict[('K','K','C')][1])
    distHh1=myco_dict[('K','K','C')][1]
    print('Hh1',Hh1)

    # order is: soil - plant - fungus - hyphae - plant growth
    Hh2=np.mean(myco_dict[('C','C','K')][1])
    sehh2=sp.stats.sem(myco_dict[('C','C','K')][1])
    distHh2=myco_dict[('C','C','K')][1]
    print('Hh2',Hh2)

    # order is: soil - plant - fungus - hyphae - plant growth
    Hg1=np.mean(myco_dict[('C','K','K')][1])
    sehg1=sp.stats.sem(myco_dict[('C','K','K')][1])
    distHg1=myco_dict[('C','K','K')][1]
    print('Hg1',Hg1)

    # order is: soil - plant - fungus - hyphae - plant growth
    Hg2=np.mean(myco_dict[('K','C','C')][1])
    sehg2=sp.stats.sem(myco_dict[('K','C','C')][1])
    distHg2=myco_dict[('K','C','C')][1]
    print('Hg2',Hg2)
    # order is: soil - plant - fungus - hyphae - plant growth
    Hs1=np.mean(myco_dict[('C','K','C')][1])
    sehs1=sp.stats.sem(myco_dict[('C','K','C')][1])
    distHs1=myco_dict[('C','K','C')][1]
    print('Hs1',Hs1)

    # order is: soil - plant - fungus - hyphae - plant growth
    Hs2=np.mean(myco_dict[('K','C','K')][1])
    sehs2= sp.stats.sem(myco_dict[('K','C','K')][1])
    distHs2=myco_dict[('K','C','K')][1]
    print('Hs2',Hs2)

    ################### fungus  parameters: #####################
    
    # order is: soil - plant - fungus - hyphae - plant growth
    Sm1=np.mean(myco_dict[('K','K','K')][0])
    print('Sm1',Sm1)
    sesm1=np.std(myco_dict[('K','K','K')][0])**2/len(myco_dict[('K','K','K')][0])
    distSm1=myco_dict[('C','C','C')][0]

    Sm2=np.mean(myco_dict[('C','C','C')][0])
    print('Sm2',Sm2)
    sesm2=np.std(myco_dict[('C','C','C')][0])**2/len(myco_dict[('C','C','C')][0])
    distSm2=myco_dict[('C','C','C')][0]

    Ss1=np.mean(myco_dict[('K','C','K')][0])
    print('Ss1',Ss1)
    sess1=np.std(myco_dict[('K','C','K')][0])**2/len(myco_dict[('K','C','K')][0])
    distSs1=myco_dict[('K','C','K')][0]

    Ss2=np.mean(myco_dict[('C','K','C')][0])
    print('Ss2',Ss2)
    sess2=np.std(myco_dict[('C','K','C')][0])**2/len(myco_dict[('C','K','C')][0])
    distSs2=myco_dict[('C','K','C')][0]

    Sg1=np.mean(myco_dict[('C','K','K')][0])
    print('Sg1',Sg1)
    sesg1=np.std(myco_dict[('C','K','K')][0])**2/len(myco_dict[('C','K','K')][0])
    distSg1=myco_dict[('C','K','K')][0]

    Sg2=np.mean(myco_dict[('K','C','C')][0])
    print('Sg2',Sg2)
    sesg2=np.std(myco_dict[('K','C','C')][0])**2/len(myco_dict[('K','C','C')][0])
    distSg2=myco_dict[('K','C','C')][0]

    Sh1=np.mean(myco_dict[('C','C','K')][0])
    print('Sh1',Sh1)
    sesh1=np.std(myco_dict[('C','C','K')][0])**2/len(myco_dict[('C','C','K')][0])
    distSh1=myco_dict[('C','C','K')][0]

    Sh2=np.mean(myco_dict[('K','K','C')][0])
    print('Sh2',Sh2)
    sesh2= np.std(myco_dict[('K','K','C')][0])**2/len(myco_dict[('K','K','C')][0])
    distSh2=myco_dict[('K','K','C')][0]
    annotationsgxes2=[['$S_{s1}$','$S_{h1}$'],['$S_{g2}$','$S_{m2}$']]
    

    # x axis labels
    gxe=['Env. K', 'Env. C']
    
    # host K relative fitnesses: 
    gxeHK=[[Hm1-Hs2,Hg1-Hh2],[Hh1-Hg2,Hs1-Hm2],['Partner K','Partner C']]# fixed for symbiont 1
    #print('DeltawKK:',Hm1-Hs2,'DeltawCK: ',Hg1-Hh2,'DeltaxwKC: ',Hh1-Hg2,'DeltawCC: ',Hs1-Hm2)
    gxeHKerr=[[(sehm1+sehs2)**0.5,(sehg1+sehh2)**0.5],[(sehh1+sehg2)**0.5,(sehs1+sehm2)**0.5]]

    ### check that host differences are significant: 

    print('host testing: ')
    test1=[distHm1,distHh1,distHs2,distHg2]
    for i in range(4):
        print(sp.stats.shapiro(test1[i]))

    print(sp.stats.bartlett(test1[0],test1[1],test1[2],test1[3]))
    print(sp.stats.ttest_ind(test1[0],test1[2]))
    print(sp.stats.ttest_ind(test1[1],test1[3]))
    # no host differences are significant


    # symbiont K relative fitnesses: 
    gxeSK=[[Sm1-Sh2,Sg1-Ss2],[Ss1-Sg2,Sh1-Sm2],['Partner K','Partner C']]# fixed for symbiont 1
    #print('DeltavKK:',Sm1-Sh2,'DeltavCK: ',Sg1-Ss2,'DeltavKC: ',Ss1-Sg2,'DeltavCC: ',Sh1-Sm2)
    gxeSKerr=[[(sesm1+sesh2)**0.5,(sesg1+sess2)**0.5],[(sesh1+sesg2)**0.5,(sesh1+sesm2)**0.5]]


    ## check that differences are significant: 

    # symbiont K relative fitnesses: 
    gxeSK=[[Sm1-Sh2,Sg1-Ss2],[Ss1-Sg2,Sh1-Sm2],['Partner K','Partner C']]# fixed for symbiont 1
    #print('DeltavKK:',Sm1-Sh2,'DeltavCK: ',Sg1-Ss2,'DeltavKC: ',Ss1-Sg2,'DeltavCC: ',Sh1-Sm2)
    gxeSKerr=[[(sesm1+sesh2)**0.5,(sesg1+sess2)**0.5],[(sesh1+sesg2)**0.5,(sesh1+sesm2)**0.5]]

    print('symbiont testing: ')
    test1=[distSm1,distSs1,distSh2,distSg2]
    for i in range(4):
        print(sp.stats.shapiro(test1[i]))
        print('length',len(test1[i]))
    
    print(test1)
    ## Sm1 vs Sh2d
    print(sp.stats.bartlett(test1[0],test1[2]))
    ## Ss1 vs Ssh2
    print(sp.stats.bartlett(test1[1],test1[3]))
    ## Sm1 vs Sh2
    print(sp.stats.ttest_ind(test1[0],test1[2]))
    ## Ss1 vs Sg2
    print(sp.stats.ttest_ind(test1[1],test1[3],equal_var=True))

    # both differences are significant. 
    fig, (ax1, ax2) = plt.subplots(1,2 ,figsize=(12, 6))
    ax1.annotate("a)", xy=(-0.05, 1.04), xycoords="axes fraction",fontsize=24)
    ax2.annotate("b)", xy=(-0.05, 1.04), xycoords="axes fraction",fontsize=24)


    fig.subplots_adjust(left=0.1,
                        bottom=0.2, 
                        right=0.9, 
                        top=0.9, 
                        wspace=0.3, 
                        hspace=0.35)


    cols=['#c80404','#0874c4']
    jit=[-0.03,0.03]
    ax1.set_ylim(-1.15,1.15)
    xmax=1.45
    xmin=-0.45
    ax1.set_xlim(xmin,xmax)

    for i,line in enumerate(gxeHK[:2]): 
        xs=[0+jit[i],1+jit[i]]
        ax1.errorbar(xs,line,gxeHKerr[i],marker='o',linestyle='-', color=cols[i])
    
    
    ax1.annotate( '$ \Delta w_{CK}$' , (0.64,-0.35),color=cols[0])
    ax1.annotate( '$ \Delta w_{CC}$' , (1.06,-0.35),color=cols[1])

    ax1.annotate( '$ \Delta w_{KK}$' , (-0.4,0.4),color=cols[0])
    ax1.annotate( '$ \Delta w_{KC}$' , (0.06,0.42),color=cols[1])

    ax1.set_xticks([0,1])
    ax1.set_xticklabels(gxe)
    
    ax1.hlines([np.average(gxeHK[0]),np.average(gxeHK[1])],xmin=xmin,xmax=xmax,colors=cols,linestyles='dashed')
    
    ax1.hlines(0,xmin=xmin,xmax=xmax,colors='k',linestyles='solid',alpha = 0.5)
   
    ax1.set_title('Host K Relative Fitness')

    ax1.set_ylabel("Difference in shoot dry mass (g)")

    ymax=55
    ymin=-55
    ax2.set_ylim(ymin,ymax)
    ax2.set_xlim(-0.15,1.15)

    for i,line in enumerate(gxeSK[:2]): 
        xs=[0+jit[i],1+jit[i]]
        ax2.errorbar(xs,line,gxeSKerr[i],marker='o',linestyle='-', color=cols[i],label=gxeSK[2][i])

    ax2.annotate( '$ \Delta v_{CK}$' , (0.88,-7),color=cols[0])
    ax2.annotate( '$ \Delta v_{CC}$' , (0.91,-34),color=cols[1])

    ax2.annotate( '$ \Delta v_{KK}$' , (0.01,45),color=cols[0])
    ax2.annotate( '$ \Delta v_{KC}$' , (-0.1,23),color=cols[1])
    ax2.set_xticks([0,1])
    ax2.set_xticklabels(gxe)
  
    
    ax2.hlines([np.average(gxeSK[0]),np.average(gxeSK[1])],xmin=-0.25,xmax=1.25,colors=cols,linestyles='dashed')
    
    ax2.hlines(0,xmin=-0.25,xmax=1.25,colors='k',linestyles='solid',alpha = 0.5)
   
    ax2.set_title('Symbiont K Relative Fitness')

    ax2.set_ylabel("Difference in hyphae density (m/g soil)")

    
    plt.figlegend(ncol=2,loc='lower center')


    plt.savefig('/Users/chriscarlson/Documents/GXE/Johnson_Delta.png',dpi=600)

    plt.show()

def Ricks_plants():


    # means and standard errors are taken directly from 
    # ricks et al. 2023. 
    # 1 = nonsaline
    #2= saline
    Hm1=122.62
    sehm1=11.63

    Hm2=30.44
    sehm2=6.539220711039238

    Hh1=135.75
    sehh1=11.177333267555014

    Hh2=67.35
    sehh2=9.311800224202571

    Hg1=26.39
    sehg1=6.793028911897245

    Hg2=111.17
    sehg2=9.95032686547924 

    Hs1=14.72
    sehs1=19.209433729306237

    Hs2=110.55
    sehs2= 10.0887671586744



    gxe=['Nonsaline Env.', 'Saline Env.']


    # translating into model's fitness advantages: 
    gxeN=[[Hm1-Hs2,Hg1-Hh2],[Hh1-Hg2,Hs1-Hm2],['Nonsaline Microbiome','Saline Microbiome']]# fixed for symbiont 1
    gxeHKerr=[[(sehm1+sehs2)**0.5,(sehg1+sehh2)**0.5],[(sehh1+sehg2)**0.5,(sehs1+sehm2)**0.5]]
    
    plt.figure(figsize=(6, 6))

    ax1 = plt.subplot()

    
    jit=[-0.03,0.03]

    cols=['#c80404','#0874c4']
    ax1.set_ylim(-50,50)
    ax1.set_xlim(-0.25,1.25)

    for i,line in enumerate(gxeN[:2]): 
        xs=[0+jit[i],1+jit[i]]
        ax1.errorbar(xs,line,gxeHKerr[i],marker='o',linestyle='-', color=cols[i],label=gxeN[2][i])

        
    ax1.annotate( '$ \Delta w_{SN}$' , (0.65,-42),color=cols[0])
    ax1.annotate( '$ \Delta w_{SS}$' , (0.77,-25),color=cols[1])

    ax1.annotate( '$ \Delta w_{NN}$' , (-0.25,18),color=cols[0])
    ax1.annotate( '$ \Delta w_{NS}$' , (0.03,30),color=cols[1])

    ax1.set_xticks([0,1])
    ax1.set_xticklabels(gxe)



    ax1.hlines([np.average(gxeN[0]),np.average(gxeN[1])],xmin=-0.25,xmax=1.25,colors=cols,linestyles='dashed')

    ax1.hlines(0,xmin=-0.25,xmax=1.25,colors='k',linestyles='solid',alpha = 0.5)

    ax1.set_title('Nonsaline Plant Relative Fitness')
    ax1.set_ylabel("Difference in plant fitness")
    plt.subplots_adjust(left=0.2,
                        bottom=0.22, 
                        right=0.9, 
                        top=0.9, 
                        wspace=0.3, 
                        hspace=0.35)
    plt.figlegend(loc='lower center')

    plt.savefig('/Users/chriscarlson/Documents/GXE/ricks_norms_Delta.png',dpi=600)
    plt.show()

def Henry_flies():

    # 1 = control
    #2= high sugar

    # CCC
    Hm1=2.6199666
    sehm1=0.1476962

    # HS on HS with HS
    Hm2=2.2886833
    sehm2=0.1278010

    # C on C with HS
    Hh1=2.6538430
    sehh1=0.1430854

    # HS on HS with C
    Hh2=1.7764400
    sehh2=0.1446986

    # C with C on HS
    Hg1=1.9381047
    sehg1=0.1260228

    # HS with HS on C
    Hg2=2.7289512
    sehg2=0.1440706

    # C with HS on HS
    Hs1=2.0133790
    sehs1=0.1226357

    # HS with C on C
    Hs2=2.6389263
    sehs2=0.1444308



    gxe=['Control Food', 'High Sugar Food']


    # Translating into fitness advantages: 
    gxeN=[[Hm1-Hs2,Hg1-Hh2],[Hh1-Hg2,Hs1-Hm2],['Control Microbiome','High Sugar Microbiome']]# fixed for symbiont 1
    gxeHKerr=[[(sehm1+sehs2)**0.5,(sehg1+sehh2)**0.5],[(sehh1+sehg2)**0.5,(sehs1+sehm2)**0.5]]
    
    
    jit=[-0.03,0.03]
    
    cols=['#c80404','#0874c4']
    plt.figure(figsize=(6, 6))
    ax1 = plt.subplot()
    
    ax1.set_ylim(-0.8,0.8)
    ax1.set_xlim(-0.25,1.25)

    for i,line in enumerate(gxeN[:2]): 
        xs=[0+jit[i],1+jit[i]]
        ax1.errorbar(xs,line,gxeHKerr[i],marker='o',linestyle='-', color=cols[i],label=gxeN[2][i])


   
    ax1.annotate( '$ \Delta w_{CH}$' , (0.65,0.22),color=cols[0])
    ax1.annotate( '$ \Delta w_{HH}$' , (0.68,-0.4),color=cols[1])

    ax1.annotate( '$ \Delta w_{CC}$' , (0.06,0.15),color=cols[0])
    ax1.annotate( '$ \Delta w_{HC}$' , (0.1,-0.3),color=cols[1])


    ax1.set_xticks([0,1])
    ax1.set_xticklabels(gxe)

    ax1.hlines([np.average(gxeN[0]),np.average(gxeN[1])],xmin=-0.25,xmax=1.25,colors=cols,linestyles='dashed')

    ax1.hlines(0,xmin=-0.25,xmax=1.25,colors='k',linestyles='solid',alpha = 0.5)

    ax1.set_title('Control Fly Relative Fitness')

    plt.subplots_adjust(left=0.2,
                        bottom=0.22, 
                        right=0.9, 
                        top=0.9, 
                        wspace=0.3, 
                        hspace=0.35)
    plt.figlegend(loc='lower center')
    ax1.set_ylabel("Difference in fecundity (egg number)")

    plt.savefig('/Users/chriscarlson/Documents/GXE/henry_norms_Delta.png',dpi=600)
    plt.show()

if __name__ == "__main__":
    myco_dict = import_Johnson_data()
    Johnson_plot(myco_dict)

    Ricks_plants()
    Henry_flies()