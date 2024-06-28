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

def Johnson_myccorhizae(myco_dict):
    print('Konza vs. Cedar Creek, Fungi')


    # order is: soil - plant - inoculum - hyphae - plant growth
    Sm1=np.mean(myco_dict[('K','K','K')][0])
    print('Sm1',Sm1)
    sesm1=sp.stats.sem(myco_dict[('A','K','K')][0])
    distSm1=myco_dict[('K','K','K')][0]

    Sm2=np.mean(myco_dict[('C','C','C')][0])
    print('Sm2',Sm2)
    sesm2=sp.stats.sem(myco_dict[('C','C','C')][0])
    distSm2=myco_dict[('C','C','C')][0]

    Ss1=np.mean(myco_dict[('K','C','K')][0])
    print('Ss1',Ss1)
    sess1=sp.stats.sem(myco_dict[('K','C','K')][0])
    distSs1=myco_dict[('K','C','K')][0]

    Ss2=np.mean(myco_dict[('C','K','C')][0])
    print('Ss2',Ss2)
    sess2=sp.stats.sem(myco_dict[('C','K','C')][0])
    distSs2=myco_dict[('C','K','C')][0]

    Sg1=np.mean(myco_dict[('C','K','K')][0])
    print('Sg1',Sg1)
    sesg1=sp.stats.sem(myco_dict[('C','K','K')][0])
    distSg1=myco_dict[('C','K','K')][0]

    Sg2=np.mean(myco_dict[('K','C','C')][0])
    print('Sg2',Sg2)
    sesg2=sp.stats.sem(myco_dict[('K','C','C')][0])
    distSg2=myco_dict[('K','C','C')][0]

    Sh1=np.mean(myco_dict[('C','C','K')][0])
    print('Sh1',Sh1)
    sesh1=sp.stats.sem(myco_dict[('C','C','K')][0])
    distSh1=myco_dict[('C','C','K')][0]

    Sh2=np.mean(myco_dict[('K','K','C')][0])
    print('Sh2',Sh2)
    sesh2= sp.stats.sem(myco_dict[('K','K','C')][0])
    distSh2=myco_dict[('K','K','C')][0]
    annotationsgxes2=[['$S_{s1}$','$S_{h1}$'],['$S_{g2}$','$S_{m2}$']]
    gxg=['Host K', 'Host C']

    gxe=['Env. K', 'Env. C']

    symbs=['symb. K', 'symb. C']



    gxge1=[[Sm1,Ss1],[Sh2,Sg2]] # fixed in environment 1
    gxge1err=[[sesm1,sess1],[sesh2,sesg2]] # fixed in environment 1
    # host 1 env 1 symbiont 1 vs. host 2 env 1 symbiont 1
    # difference two is 0.012
    gxge2=[[Sg1,Sh1],[Ss2,Sm2]] # fixed in environment 2
    gxge2err=[[sesg1,sesh1],[sess2,sesm2]]


    gxes1=[[Sm1,Sg1],[Sh2,Ss2]]# fixed for symbiont 1
    gxes1err=[[sesm1,sesg1],[sesh2,sess2]]

    gxes2=[[Ss1,Sh1],[Sg2,Sm2]]# fixed for symbiont 2
    gxes2err=[[sess1,sesh1],[sesg2,sesm2]]

    cols=['#c80404','#0874c4']

    print('panel A testing: ')
    test1=[distSm1,distSs1,distSh2,distSg2]
    for i in range(4):
        print(sp.stats.shapiro(test1[i]))
        print('length',len(test1[i]))
    print(test1)
    print(sp.stats.bartlett(test1[0],test1[2]))
    print(sp.stats.bartlett(test1[1],test1[3]))
    print(sp.stats.ttest_ind(test1[0],test1[2]))
    print(sp.stats.ttest_ind(test1[1],test1[3],equal_var=True))

    # SMF and SsK are significantly different
    # SMF and SsK are significantly different


    print('panel C testing: ')
    test1=[distSg1,distSh1,distSs2,distSm2]
    for i in range(4):
        print(sp.stats.shapiro(test1[i]))
        print('length',len(test1[i]))

    print(sp.stats.bartlett(test1[0],test1[2]))
    print(sp.stats.bartlett(test1[1],test1[3]))
    print(sp.stats.ttest_ind(test1[0],test1[2]))
    print(sp.stats.ttest_ind(test1[1],test1[3],equal_var=True))

    # Sg1 and Ss2 are significantly different as are sh1 nd sm2





    print('panel B testing: ')
    test1=[distSm1,distSg1,distSh2,distSs2]
    for i in range(4):
        print(sp.stats.shapiro(test1[i]))
        print('length',len(test1[i]))

    print(sp.stats.levene(test1[0],test1[2]))
    print(sp.stats.bartlett(test1[1],test1[3]))
    print(sp.stats.mannwhitneyu(test1[0],test1[2]))
    print(sp.stats.ttest_ind(test1[1],test1[3],equal_var=True))


    # both are significant


    print('panel D testing: ')
    test1=[distSs1,distSh1,distSg2,distSm2]
    for i in range(4):
        print(sp.stats.shapiro(test1[i]))
        print('length',len(test1[i]))

    print(sp.stats.bartlett(test1[0],test1[2]))
    print(sp.stats.bartlett(test1[1],test1[3]))
    print(sp.stats.ttest_ind(test1[0],test1[2]))
    print(sp.stats.ttest_ind(test1[1],test1[3],equal_var=True))

    # both significant


    fig, ((ax1, ax2), (ax3, ax4)) = plt.subplots(2, 2,figsize=(10, 12))

    ax1.annotate("a)", xy=(-0.05, 1.04), xycoords="axes fraction",fontsize=24)
    ax2.annotate("b)", xy=(-0.05, 1.04), xycoords="axes fraction",fontsize=24)
    ax3.annotate("c)", xy=(-0.05, 1.04), xycoords="axes fraction",fontsize=24)
    ax4.annotate("d)", xy=(-0.05, 1.04), xycoords="axes fraction",fontsize=24)


    fig.subplots_adjust(left=0.1,
                        bottom=0.1, 
                        right=0.9, 
                        top=0.95, 
                        wspace=0.3, 
                        hspace=0.35)

    ax1.annotate("a)", xy=(-0.05, 1.04), xycoords="axes fraction",fontsize=24)
    ax2.annotate("b)", xy=(-0.05, 1.04), xycoords="axes fraction",fontsize=24)
    ax3.annotate("c)", xy=(-0.05, 1.04), xycoords="axes fraction",fontsize=24)
    ax4.annotate("d)", xy=(-0.05, 1.04), xycoords="axes fraction",fontsize=24)


    annotationsgxge1=[['$S_{mK}$','$S_{sK}$'],['$S_{hC}$','$S_{gC}$']]
    ax1.set_ylim(30,100)
    ax1.set_xlim(-0.25,1.25)
    jit=[-0.03,0.03]



    for i,line in enumerate(gxge1): 
        xs=[0+jit[i],1+jit[i]]
        ax1.errorbar(xs,line,gxge1err[i],marker='o',linestyle='-', color=cols[i], label=symbs[i])
        for j,label in enumerate(annotationsgxge1[i]):
            for j,label in enumerate(annotationsgxes2[i]):
                if i==0 and j ==0:
                    ax1.annotate(label,(xs[j]-0.2,line[j]+0.1))
                if i==0 and j ==1:
                    ax1.annotate(label,(xs[j]+0.035,line[j]+0.05))
                if i==1 and j ==0:
                    ax1.annotate(label,(xs[j]-0.2,line[j]+0.1))
                if i==1 and j ==1:
                    ax1.annotate(label,(xs[j]+0.035,line[j]+0.05))

    #ax1.annotate(r'$\underline{n.s.}$',[0,129])
    # p = 0.233
    ax1.annotate(r'$\underline{*}$',[1.01,144])
    ax1.set_xticks([0,1])
    ax1.set_xticklabels(gxg)
    # p = 0.012
    ax1.set_xlabel('Environment K')
    #ax1.set_ylabel('Fitness')
    #ax1.legend()
    #ax.legend()
    ax1.set_title('GXG Environment K')

    annotationsgxge2=[['$S_{gK}$','$S_{hK}$'],['$S_{sC}$','$S_{mC}$']]
    ax3.set_ylim(30,100)
    ax3.set_xlim(-0.25,1.25)


    for i,line in enumerate(gxge2): 
        xs=[0+jit[i],1+jit[i]]
        ax3.errorbar(xs,line,gxge2err[i],marker='o',linestyle='-', color=cols[i])
        for j,label in enumerate(annotationsgxge2[i]):
            if i==0 and j ==0:
                ax3.annotate(label,(xs[j]-0.2,line[j]+0.1))
            if i==0 and j ==1:
                ax3.annotate(label,(xs[j]+0.035,line[j]+0.05))
            if i==1 and j ==0:
                ax3.annotate(label,(xs[j]-0.2,line[j]+0.1))
            if i==1 and j ==1:
                ax3.annotate(label,(xs[j]+0.035,line[j]+0.05))

    #ax3.annotate(r'$\underline{*}$',[0,77])
    # p < 0.001
    #ax3.annotate(r'$\underline{*}$',[1,40])
    # p = 0.006
    ax3.set_xticks([0,1])
    ax3.set_xticklabels(gxg)

    ax3.set_xlabel('Environment C')
    #ax3.set_ylabel('Fitness')
    #ax3.legend(loc='lower left')
    ax3.set_title('GXG Environment C')





    annotationsgxes1=[['$S_{mK}$','$S_{gK}$'],['$S_{hC}$','$S_{sC}$']]
    ax2.set_ylim(30,100)
    ax2.set_xlim(-0.25,1.25)

    for i,line in enumerate(gxes1): 
        xs=[0+jit[i],1+jit[i]]
        ax2.errorbar(xs,line,gxes1err[i],marker='o',linestyle='-', color=cols[i])
        for j,label in enumerate(annotationsgxes1[i]):
            if i==0 and j ==0:
                ax2.annotate(label,(xs[j]-0.21,line[j]+2))
            if i==0 and j ==1:
                ax2.annotate(label,(xs[j]+0.037,line[j]+0.05))
            if i==1 and j ==0:
                ax2.annotate(label,(xs[j]-0.2,line[j]+1.5))
            if i==1 and j ==1:
                ax2.annotate(label,(xs[j]+0.037,line[j]+0.05))

    #ax2.annotate(r'$\underline{*}$',[1,84])
    # p = 0.001
    ax2.set_xticks([0,1])
    ax2.set_xticklabels(gxe)

    ax2.set_xlabel('Host K')
    #ax2.set_ylabel('Fitness')
    #ax2.legend()
    ax2.set_title('GXE Host K')

    annotationsgxes2=[['$S_{sK}$','$S_{hK}$'],['$S_{gC}$','$S_{mC}$']]
    ax4.set_ylim(30,100)
    ax4.set_xlim(-0.25,1.25)

    for i,line in enumerate(gxes2): 
        xs=[0+jit[i],1+jit[i]]
        ax4.errorbar(xs,line,gxes2err[i],marker='o',linestyle='-', color=cols[i])
        for j,label in enumerate(annotationsgxes2[i]):
            if i==0 and j ==0:
                ax4.annotate(label,(xs[j]-0.2,line[j]+0.1))
            if i==0 and j ==1:
                ax4.annotate(label,(xs[j]+0.037,line[j]+0.05))
            if i==1 and j ==0:
                ax4.annotate(label,(xs[j]-0.2,line[j]+0.1))
            if i==1 and j ==1:
                ax4.annotate(label,(xs[j]+0.037,line[j]+0.05))


    # p = 0.012
    ax4.set_xticks([0,1])
    ax4.set_xticklabels(gxe)
    ax4.set_xlabel('Host C')
    #ax4.set_ylabel('Fitness')
    #ax4.legend(loc='lower right')
    ax4.set_title('GXE Host C')

    print('DeltaSmh 1',Sm1-Sh2)
    print('DeltaSmh 2',Sm2-Sh1)

    print('DeltaSgs 1',Sg1-Ss2)
    print('DeltaSgs 2',Sg2-Ss1)

    print('total',Sm1-Sh2+Sm2-Sh1+Sg1-Ss2+Sg2-Ss1)

    print('h1',Sm1-Sh2+Sg1-Ss2)
    print('h2',Sm2-Sh1+Sg2-Ss1)

    fig.supylabel("Hyphae Density (m/g soil)")
    plt.figlegend(ncol=2,loc='center')
    plt.tight_layout(h_pad=2.2)

    plt.savefig('/Users/chriscarlson/Documents/GXE/Johnson_fungi_norms.png',dpi=600)

    plt.show()

def Johnson_plants(myco_dict):

    print('Konza vs. Cedar Creek, Plants')

    Hm1=np.mean(myco_dict[('K','K','K')][1])
    sehm1=sp.stats.sem(myco_dict[('K','K','K')][1])
    distHm1=myco_dict[('K','K','K')][1]
    print('Hm1',Hm1)

    Hm2=np.mean(myco_dict[('C','C','C')][1])
    sehm2=sp.stats.sem(myco_dict[('C','C','C')][1])
    distHm2=myco_dict[('C','C','C')][1]
    print('Hm2',Hm2)

    Hh1=np.mean(myco_dict[('K','K','C')][1])
    sehh1=sp.stats.sem(myco_dict[('K','K','C')][1])
    distHh1=myco_dict[('K','K','C')][1]
    print('Hh1',Hh1)

    Hh2=np.mean(myco_dict[('C','C','K')][1])
    sehh2=sp.stats.sem(myco_dict[('C','C','K')][1])
    distHh2=myco_dict[('C','C','K')][1]
    print('Hh2',Hh2)

    Hg1=np.mean(myco_dict[('C','K','K')][1])
    sehg1=sp.stats.sem(myco_dict[('C','K','K')][1])
    distHg1=myco_dict[('C','K','K')][1]
    print('Hg1',Hg1)

    Hg2=np.mean(myco_dict[('K','C','C')][1])
    sehg2=sp.stats.sem(myco_dict[('K','C','C')][1])
    distHg2=myco_dict[('K','C','C')][1]
    print('Hg2',Hg2)

    Hs1=np.mean(myco_dict[('C','K','C')][1])
    sehs1=sp.stats.sem(myco_dict[('C','K','C')][1])
    distHs1=myco_dict[('C','K','C')][1]
    print('Hs1',Hs1)


    Hs2=np.mean(myco_dict[('K','C','K')][1])
    sehs2= sp.stats.sem(myco_dict[('K','C','K')][1])
    distHs2=myco_dict[('K','C','K')][1]
    print('Hs2',Hs2)

    gxg=['Symbiont K', 'Symbiont C']

    gxe=['Env. K', 'Env. C']

    hosts=['host K', 'host C']




    gxge1=[[Hm1,Hh1],[Hs2,Hg2]] # fixed in environment 1
    gxge1err=[[sehm1,sehh1],[sehs2,sehg2]] # fixed in environment 1
    # host 1 env 1 symbiont 1 vs. host 2 env 1 symbiont 1
    # difference two is 0.012
    gxge2=[[Hg1,Hs1],[Hh2,Hm2]] # fixed in environment 2
    gxge2err=[[sehg1,sehs1],[sehh2,sehm2]]


    gxes1=[[Hm1,Hg1],[Hs2,Hh2]]# fixed for symbiont 1
    gxes1err=[[sehm1,sehg1],[sehs2,sehh2]]

    gxes2=[[Hh1,Hs1],[Hg2,Hm2]]# fixed for symbiont 2
    gxes2err=[[sehh1,sehs1],[sehg2,sehm2]]

    cols=['#c80404','#0874c4']

    print('panel A testing: ')
    test1=[distHm1,distHh1,distHs2,distHg2]
    for i in range(4):
        print(sp.stats.shapiro(test1[i]))

    print(sp.stats.bartlett(test1[0],test1[1],test1[2],test1[3]))
    print(sp.stats.ttest_ind(test1[0],test1[2]))
    print(sp.stats.ttest_ind(test1[1],test1[3]))



    # no significant differences

    print('panel C testing: ')
    test1=[distHg1,distHs1,distHh2,distHm2]
    for i in range(4):
        print(sp.stats.shapiro(test1[i]))

    print(sp.stats.bartlett(test1[0],test1[1],test1[2],test1[3]))
    print(sp.stats.ttest_ind(test1[0],test1[2],equal_var=False))
    print(sp.stats.ttest_ind(test1[1],test1[3],equal_var=False))

    # No significant differences


    print('panel B testing: ')
    test1=[distHm1,distHg1,distHs2,distHh2]
    for i in range(4):
        print(sp.stats.shapiro(test1[i]))

    print(sp.stats.bartlett(test1[0],test1[1],test1[2],test1[3]))
    print(sp.stats.ttest_ind(test1[0],test1[2],equal_var=False))
    print(sp.stats.ttest_ind(test1[1],test1[3],equal_var=False))

    # No significant differences


    print('panel D testing: ')
    test1=[distHh1,distHs1,distHg2,distHm2]
    for i in range(4):
        print(sp.stats.shapiro(test1[i]))

    print(sp.stats.bartlett(test1[0],test1[1],test1[2],test1[3]))
    print(sp.stats.ttest_ind(test1[0],test1[2],equal_var=False))
    print(sp.stats.ttest_ind(test1[1],test1[3],equal_var=False))

    # both diffreences are insignificant. 

    fig, ((ax1, ax2), (ax3, ax4)) = plt.subplots(2, 2,figsize=(10, 10))

    fig.subplots_adjust(left=0.1,
                        bottom=0.1, 
                        right=0.9, 
                        top=0.95, 
                        wspace=0.3, 
                        hspace=0.35)

    ax1.annotate("a)", xy=(-0.05, 1.04), xycoords="axes fraction",fontsize=24)
    ax2.annotate("b)", xy=(-0.05, 1.04), xycoords="axes fraction",fontsize=24)
    ax3.annotate("c)", xy=(-0.05, 1.04), xycoords="axes fraction",fontsize=24)
    ax4.annotate("d)", xy=(-0.05, 1.04), xycoords="axes fraction",fontsize=24)

    annotationsgxge1=[['$H_{mK}$','$H_{hK}$'],['$H_{sC}$','$H_{gC}$']]
    ax1.set_ylim(2,7)
    ax1.set_xlim(-0.25,1.25)
    jit=[-0.03,0.03]

    for i,line in enumerate(gxge1): 
        xs=[0+jit[i],1+jit[i]]
        ax1.errorbar(xs,line,gxge1err[i],marker='o',linestyle='-', color=cols[i], label=hosts[i])
        for j,label in enumerate(annotationsgxge1[i]):
            if i==0 and j ==0:
                ax1.annotate(label,(xs[j]-0.22,line[j]+0.35))
            if i==0 and j ==1:
                ax1.annotate(label,(xs[j]+0.035,line[j]+0.3))
            if i==1 and j ==0:
                ax1.annotate(label,(xs[j]-0.22,line[j]-0.3))
            if i==1 and j ==1:
                ax1.annotate(label,(xs[j]+0.035,line[j]+0.05))

    #ax1.annotate(r'$\underline{n.s.}$',[0,129])
    # p = 0.233
    ax1.annotate(r'$\underline{*}$',[1.01,144])
    ax1.set_xticks([0,1])
    ax1.set_xticklabels(gxg)
    # p = 0.012
    ax1.set_xlabel('Env. K')
    #ax1.set_ylabel('Fitness')
    ax1.legend(loc='upper left')
    #ax.legend()
    ax1.set_title('GXG Environment K')


    #plt.show()

    #fig, ax = plt.subplots()
    annotationsgxge2=[['$H_{gK}$','$H_{sK}$'],['$H_{hC}$','$H_{mC}$']]
    ax3.set_ylim(1.5,4.5)
    ax3.set_xlim(-0.25,1.25)


    for i,line in enumerate(gxge2): 
        xs=[0+jit[i],1+jit[i]]
        ax3.errorbar(xs,line,gxge2err[i],marker='o',linestyle='-', color=cols[i], label=hosts[i])
        for j,label in enumerate(annotationsgxge2[i]):
            if i==0 and j ==0:
                ax3.annotate(label,(xs[j]-0.22,line[j]+0.05))
            if i==0 and j ==1:
                ax3.annotate(label,(xs[j]+0.075,line[j]-0.05))
            if i==1 and j ==0:
                ax3.annotate(label,(xs[j]-0.22,line[j]+0.1))
            if i==1 and j ==1:
                ax3.annotate(label,(xs[j]+0.0125,line[j]+0.1))

    ax3.annotate(r'$\underline{*}$',[0,77])
    # p < 0.001
    ax3.annotate(r'$\underline{*}$',[1,40])
    # p = 0.006
    ax3.set_xticks([0,1])
    ax3.set_xticklabels(gxg)

    ax3.set_xlabel('Environment C')
    #ax3.set_ylabel('Fitness')
    ax3.legend()
    ax3.set_title('GXG Environment C')



    annotationsgxes1=[['$H_{mK}$','$H_{gK}$'],['$H_{sC}$','$H_{hC}$']]
    ax2.set_ylim(2,6)
    ax2.set_xlim(-0.25,1.25)

    for i,line in enumerate(gxes1): 
        xs=[0+jit[i],1+jit[i]]
        ax2.errorbar(xs,line,gxes1err[i],marker='o',linestyle='-', color=cols[i], label=hosts[i])
        for j,label in enumerate(annotationsgxes1[i]):
            if i==0 and j ==0:
                ax2.annotate(label,(xs[j]-0.22,line[j]+0.35))
            if i==0 and j ==1:
                ax2.annotate(label,(xs[j]+0.08,line[j]))
            if i==1 and j ==0:
                ax2.annotate(label,(xs[j]-0.22,line[j]-0.3))
            if i==1 and j ==1:
                ax2.annotate(label,(xs[j]+0.035,line[j]+0.15))

    ax2.annotate(r'$\underline{*}$',[1,84])
    # p = 0.001
    ax2.set_xticks([0,1])
    ax2.set_xticklabels(gxe)

    ax2.set_xlabel('Symbiont K')
    #ax2.set_ylabel('Fitness')
    ax2.legend()
    ax2.set_title('GXE Symbiont K')





    annotationsgxes2=[['$H_{hK}$','$H_{sK}$'],['$H_{gC}$','$H_{mC}$']]
    ax4.set_ylim(2,5)
    ax4.set_xlim(-0.25,1.25)

    for i,line in enumerate(gxes2): 
        xs=[0+jit[i],1+jit[i]]
        ax4.errorbar(xs,line,gxes2err[i],marker='o',linestyle='-', color=cols[i], label=hosts[i])
        for j,label in enumerate(annotationsgxes2[i]):
            if i==0 and j ==0:
                ax4.annotate(label,(xs[j]-0.22,line[j]+0.35))
            if i==0 and j ==1:
                ax4.annotate(label,(xs[j]+0.08,line[j]))
            if i==1 and j ==0:
                ax4.annotate(label,(xs[j]-0.26,line[j]-0.3))
            if i==1 and j ==1:
                ax4.annotate(label,(xs[j]+0.01,line[j]+0.15))

    ax4.annotate(r'$\underline{*}$',[0,148])
    # p = 0.012

    ax4.annotate(r'$\underline{*}$',[1,48])
    # p = 0.012
    ax4.set_xticks([0,1])
    ax4.set_xticklabels(gxe)
    ax4.set_xlabel('Symbiont C')
    #ax4.set_ylabel('Fitness')
    ax4.legend()
    ax4.set_title('GXE Symbiont C')

    fig.supylabel("Shoot Dry Mass (g)")
    plt.savefig('/Users/chriscarlson/Documents/GXE/Johnson_plant_norms.png',dpi=600)
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

    gxg=['Symbiont 1', 'Symbiont 2']

    gxe=['Env. 1', 'Env. 2']

    hosts=['Host 1', 'Host 2']



    gxge1=[[Hm1,Hh1],[Hs2,Hg2]] # fixed in environment 1
    gxge1err=[[sehm1,sehh1],[sehs2,sehg2]] # fixed in environment 1
    # host 1 env 1 symbiont 1 vs. host 2 env 1 symbiont 1
    # difference two is 0.012
    gxge2=[[Hg1,Hs1],[Hh2,Hm2]] # fixed in environment 2
    gxge2err=[[sehg1,sehs1],[sehh2,sehm2]]


    gxes1=[[Hm1,Hg1],[Hs2,Hh2]]# fixed for symbiont 1
    gxes1err=[[sehm1,sehg1],[sehs2,sehh2]]

    gxes2=[[Hh1,Hs1],[Hg2,Hm2]]# fixed for symbiont 2
    gxes2err=[[sehh1,sehs1],[sehg2,sehm2]]

    cols=['#c80404','#0874c4']
    fig, ((ax1, ax2), (ax3, ax4)) = plt.subplots(2, 2,figsize=(10, 10))

    fig.subplots_adjust(left=0.1,
                        bottom=0.1, 
                        right=0.9, 
                        top=0.95, 
                        wspace=0.3, 
                        hspace=0.35)

    annotationsgxge1=[['$H_{m1}$','$H_{h1}$'],['$H_{s1}$','$H_{g1}$']]
    ax1.set_ylim(100,max(max(gxge1))+20)
    ax1.set_xlim(-0.25,1.25)
    jit=[-0.03,0.03]

    for i,line in enumerate(gxge1): 
        xs=[0+jit[i],1+jit[i]]
        ax1.errorbar(xs,line,gxge1err[i],marker='o',linestyle='-', color=cols[i], label=hosts[i])
        for j,label in enumerate(annotationsgxge1[i]):
            ax1.annotate(label,(xs[j]+0.025,line[j]+3.5))

    #ax.annotate(r'$\underline{n.s.}$',[0,129])
    # p = 0.233
    ax1.annotate('**',[1.01,144])
    ax1.set_xticks([0,1])
    ax1.set_xticklabels(gxg)
    # p = 0.012
    ax1.set_xlabel('Environment 1')
    #ax1.set_ylabel('Fitness')
    ax1.legend(loc='upper left')
    #ax.legend()
    ax1.set_title('GXG Environment 1',fontweight="bold")

    ax1.annotate("a)", xy=(-0.05, 1.04), xycoords="axes fraction",fontsize=26)


    #plt.show()

    #fig, ax = plt.subplots()
    annotationsgxge2=[['$H_{g1}$','$H_{s1}$'],['$H_{h2}$','$H_{m2}$']]
    ax3.set_ylim(14.72-20,max(max(gxge2))+14)
    ax3.set_xlim(-0.25,1.25)


    for i,line in enumerate(gxge2): 
        xs=[0+jit[i],1+jit[i]]
        ax3.errorbar(xs,line,gxge2err[i],marker='o',linestyle='-', color=cols[i], label=hosts[i])
        for j,label in enumerate(annotationsgxge2[i]):
            ax3.annotate(label,(xs[j]+0.015,line[j]+2))

    ax3.annotate('**',[0,77])
    # p < 0.001
    ax3.annotate('*',[1,40])
    # p = 0.006
    ax3.set_xticks([0,1])
    ax3.set_xticklabels(gxg)

    ax3.set_xlabel('Environment 2')
    #ax3.set_ylabel('Fitness')
    ax3.legend()
    ax3.set_title('GXG Environment 2',fontweight="bold")

    ax3.annotate("c)", xy=(-0.05, 1.04), xycoords="axes fraction",fontsize=26)



    annotationsgxes1=[['$H_{m1}$','$H_{g1}$'],['$H_{s2}$','$H_{h2}$']]
    ax2.set_ylim(26.39-10,max(max(gxes1))+14)
    ax2.set_xlim(-0.25,1.25)

    for i,line in enumerate(gxes1): 
        xs=[0+jit[i],1+jit[i]]
        ax2.errorbar(xs,line,gxes1err[i],marker='o',linestyle='-', color=cols[i], label=hosts[i])
        for j,label in enumerate(annotationsgxes1[i]):
            ax2.annotate(label,(xs[j]+0.025,line[j]+3))

    ax2.annotate('**',[1,84])
    # p = 0.001
    ax2.set_xticks([0,1])
    ax2.set_xticklabels(gxe)

    ax2.set_xlabel('Symbiont 1')
    #ax2.set_ylabel('Fitness')
    ax2.legend()
    ax2.set_title('GXE Symbiont 1',fontweight="bold")
    ax2.annotate("b)", xy=(-0.05, 1.04), xycoords="axes fraction",fontsize=26)




    annotationsgxes2=[['$H_{h1}$','$H_{s1}$'],['$H_{g2}$','$H_{m2}$']]
    ax4.set_ylim(14.6-20,max(max(gxes2))+20)
    ax4.set_xlim(-0.25,1.25)

    for i,line in enumerate(gxes2): 
        xs=[0+jit[i],1+jit[i]]
        ax4.errorbar(xs,line,gxes2err[i],marker='o',linestyle='-', color=cols[i], label=hosts[i])
        for j,label in enumerate(annotationsgxes2[i]):
            if label=='$H_{g2}$':
                #print('hooray')
                ax4.annotate(label,(xs[j]-0.2,line[j]+2))
            else: 
                ax4.annotate(label,(xs[j]+0.015,line[j]+2))

    ax4.annotate('**',[0,148])
    # p = 0.012

    ax4.annotate('*',[1,48])
    # p = 0.012
    ax4.set_xticks([0,1])
    ax4.set_xticklabels(gxe)
    ax4.set_xlabel('Symbiont 2')
    #ax4.set_ylabel('Fitness')
    ax4.legend()
    ax4.set_title('GXE Symbiont 2',fontweight="bold")
    ax4.annotate("d)", xy=(-0.05, 1.04), xycoords="axes fraction",fontsize=26)

    fig.supylabel("Plant Fitness")




    plt.savefig('/Users/chriscarlson/Documents/GXE/ricks_norms.png',dpi=600)
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

    gxg=['Symbiont 1', 'Symbiont 2']

    gxe=['Env. 1', 'Env. 2']

    hosts=['Host 1', 'Host 2']

    # Significant differences (as calculated by Henry et al. ): 
    # Hg1 - Hs1
    # Hg1 - hm2
    # hs1 - hm2
    # hh2 - hm2


    gxge1=[[Hm1,Hh1],[Hs2,Hg2]] # fixed in environment 1
    gxge1err=[[sehm1,sehh1],[sehs2,sehg2]] # fixed in environment 1
    # host 1 env 1 symbiont 1 vs. host 2 env 1 symbiont 1
    # difference two is 0.012
    gxge2=[[Hg1,Hs1],[Hh2,Hm2]] # fixed in environment 2
    gxge2err=[[sehg1,sehs1],[sehh2,sehm2]]


    gxes1=[[Hm1,Hg1],[Hs2,Hh2]]# fixed for symbiont 1
    gxes1err=[[sehm1,sehg1],[sehs2,sehh2]]

    gxes2=[[Hh1,Hs1],[Hg2,Hm2]]# fixed for symbiont 2
    gxes2err=[[sehh1,sehs1],[sehg2,sehm2]]


    cols=['#c80404','#0874c4']
    fig, ((ax1, ax2), (ax3, ax4)) = plt.subplots(2, 2,figsize=(10, 10))

    fig.subplots_adjust(left=0.1,
                        bottom=0.1, 
                        right=0.9, 
                        top=0.95, 
                        wspace=0.3, 
                        hspace=0.35)

    annotationsgxge1=[['$H_{m1}$','$H_{h1}$'],['$H_{s2}$','$H_{g2}$']]
    ax1.set_ylim(2,3)
    ax1.set_xlim(-0.25,1.25)
    jit=[-0.03,0.03]
    count=0
    for i,line in enumerate(gxge1): 
        xs=[0+jit[i],1+jit[i]]
        ax1.errorbar(xs,line,gxge1err[i],marker='o',linestyle='-', color=cols[i], label=hosts[i])
        for j,label in enumerate(annotationsgxge1[i]):
            if i==0 and j ==0:
                ax1.annotate(label,(xs[j]-0.2,line[j]-0.075))
            if i==0 and j ==1:
                ax1.annotate(label,(xs[j]+0.075,line[j]))
            if i==1 and j ==0:
                ax1.annotate(label,(xs[j]+0.05,line[j]+0.06))
            if i==1 and j ==1:
                ax1.annotate(label,(xs[j]+0.05,line[j]+0.01))
            
            
            


    #ax.annotate(r'$\underline{n.s.}$',[0,129])
    # p = 0.233
    ax1.annotate('**',[1.01,144])
    ax1.set_xticks([0,1])
    ax1.set_xticklabels(gxg)
    # p = 0.012
    ax1.set_xlabel('Environment 1')
    #ax1.set_ylabel('Fitness')
    ax1.legend(loc='lower left')
    #ax.legend()
    ax1.set_title('GXG Environment 1',fontweight="bold")

    ax1.annotate("a)", xy=(-0.05, 1.04), xycoords="axes fraction",fontsize=26)


    #plt.show()

    #fig, ax = plt.subplots()
    annotationsgxge2=[['$H_{g1}$','$H_{s1}$'],['$H_{h2}$','$H_{m2}$']]
    ax3.set_ylim(1,3)
    ax3.set_xlim(-0.25,1.25)


    for i,line in enumerate(gxge2): 
        xs=[0+jit[i],1+jit[i]]
        ax3.errorbar(xs,line,gxge2err[i],marker='o',linestyle='-', color=cols[i], label=hosts[i])
        for j,label in enumerate(annotationsgxge2[i]):
            if i==0 and j ==0:
                ax3.annotate(label,(xs[j]-0.05,line[j]+0.15))
            if i==0 and j ==1:
                ax3.annotate(label,(xs[j]+0.035,line[j]+0.015))
            if i==1 and j ==0:
                ax3.annotate(label,(xs[j]-0.1,line[j]-0.25))
            if i==1 and j ==1:
                ax3.annotate(label,(xs[j]+0.025,line[j]+0.02))


    # p < 0.001
    ax3.annotate('*',[1,2.5])
    # p = 0.006
    ax3.set_xticks([0,1])
    ax3.set_xticklabels(gxg)

    ax3.set_xlabel('Environment 2')
    #ax3.set_ylabel('Fitness')
    ax3.legend(loc='upper left')
    ax3.set_title('GXG Environment 2',fontweight="bold")

    ax3.annotate("c)", xy=(-0.05, 1.04), xycoords="axes fraction",fontsize=26)



    annotationsgxes1=[['$H_{m1}$','$H_{g1}$'],['$H_{s2}$','$H_{h2}$']]
    ax2.set_ylim(1,3)
    ax2.set_xlim(-0.25,1.25)

    for i,line in enumerate(gxes1): 
        xs=[0+jit[i],1+jit[i]]
        ax2.errorbar(xs,line,gxes1err[i],marker='o',linestyle='-', color=cols[i], label=hosts[i])
        for j,label in enumerate(annotationsgxes1[i]):
            if i==0 and j ==0:
                ax2.annotate(label,(xs[j]-0.2,line[j]+0.075))
            if i==0 and j ==1:
                ax2.annotate(label,(xs[j]+0.025,line[j]+0.05))
            if i==1 and j ==0:
                ax2.annotate(label,(xs[j]+0.025,line[j]+0.075))
            if i==1 and j ==1:
                ax2.annotate(label,(xs[j]+0.025,line[j]+0.05))



    # p = 0.001
    ax2.set_xticks([0,1])
    ax2.set_xticklabels(gxe)

    ax2.set_xlabel('Symbiont 1')
    #ax2.set_ylabel('Fitness')
    ax2.legend()
    ax2.set_title('GXE Symbiont 1',fontweight="bold")
    ax2.annotate("b)", xy=(-0.05, 1.04), xycoords="axes fraction",fontsize=26)


    annotationsgxes2=[['$H_{h1}$','$H_{s1}$'],['$H_{g2}$','$H_{m2}$']]
    ax4.set_ylim(1,3)
    ax4.set_xlim(-0.25,1.25)

    for i,line in enumerate(gxes2): 
        xs=[0+jit[i],1+jit[i]]
        ax4.errorbar(xs,line,gxes2err[i],marker='o',linestyle='-', color=cols[i], label=hosts[i])
        for j,label in enumerate(annotationsgxes2[i]):
            if i==0 and j ==0:
                ax4.annotate(label,(xs[j]-0.2,line[j]+0.075))
            if i==0 and j ==1:
                ax4.annotate(label,(xs[j]+0.025,line[j]+0.05))
            if i==1 and j ==0:
                ax4.annotate(label,(xs[j]+0.025,line[j]+0.075))
            if i==1 and j ==1:
                ax4.annotate(label,(xs[j]+0.025,line[j]+0.05))


    # p = 0.012

    ax4.annotate('*',[1,2.6])
    # p = 0.012 
    ax4.set_xticks([0,1])
    ax4.set_xticklabels(gxe)
    ax4.set_xlabel('Symbiont 2')
    #ax4.set_ylabel('Fitness')
    ax4.legend()
    ax4.set_title('GXE Symbiont 2',fontweight="bold")
    ax4.annotate("d)", xy=(-0.05, 1.04), xycoords="axes fraction",fontsize=26)

    fig.supylabel("Fly Fecundity")

    plt.savefig('/Users/chriscarlson/Documents/GXE/henry_norms.png',dpi=600)
    plt.show()

if __name__ == "__main__":
    myco_dict = import_Johnson_data()
    #Johnson_myccorhizae(myco_dict)
    Johnson_plants(myco_dict)
    #Ricks_plants()
    #Henry_flies()

