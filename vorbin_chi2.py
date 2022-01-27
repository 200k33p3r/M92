#This code requires two inputs: mc_num and s/n
#This code outputs a list of chi2 for each age, dm and reddening combination
import numpy as np
import pandas as pd
import sys
import os

class chi2:

    def read_input(self):
        #read vorbin file
        dp = pd.read_csv("/home/mying/Desktop/vorbin_chi2/sn={}".format(self.sn))
        self.xBar = dp['xBar'].values
        self.yBar = dp['yBar'].values
        #read observed bin counts
        dp = pd.read_csv("/home/mying/Desktop/vorbin_chi2/bin_count_std_sn={}".format(self.sn))
        self.bin_count_std = dp['bin_count_std'].values
        self.obs_size = np.sum(self.bin_count_std)

    #search for bin number for each data point
    def search_point_location_bc(self,x, y, xBar, yBar):
        lenx = len(x)
        delta_x = xBar - x.reshape(lenx,1)
        delta_y = yBar - y.reshape(lenx,1)
        distance = np.square(delta_x) + np.square(delta_y)
        bin_num = np.argmin(distance, axis = 1)
        return bin_num
    def search_vorbin(self):
        #go through the search process
        Tb_size = 1000
        dms = np.linspace(14.62,14.82,21)
        reds = np.linspace(0.0,0.05,6)
        ages = ['09000','10000','11000','12000','13000','14000','15000','16000']
        chi2 = []
        for age in ages:
            #read iso files
            dp = pd.read_csv("/data/outcmd/mc{}.a{}".format(self.mc_num,age),sep='\s+',names=['vi','v','i'],skiprows=3)
            for dm in dms:
                for red in reds:
                    bin_count = np.zeros(len(self.xBar))
                    total_pt = len(dp)
                    n_div = total_pt // Tb_size
                    for i in range(n_div):
                        bin_num = self.search_point_location_bc((dp['vi'].values[i*Tb_size:(i+1)*Tb_size] + red)*12.5, dp['v'].values[i*Tb_size:(i+1)*Tb_size] + dm, self.xBar, self.yBar)
                        for j in range(Tb_size):
                            bin_count[bin_num[j]] += 1
        #do the last bit
                    bin_num = self.search_point_location_bc((dp['vi'].values[n_div*Tb_size:] + red)*12.5, dp['v'].values[n_div*Tb_size:] + dm, self.xBar, self.yBar)
                    for j in range(total_pt - n_div*Tb_size):
                        bin_count[bin_num[j]] += 1
                    #to avoid divde by 0
                    for i in range(len(bin_count)):
                        if bin_count[i] == 0:
                            bin_count[i] += 1
                    chi2.append([age, dm, red, np.inner(np.divide(self.bin_count_std,bin_count/(total_pt/self.obs_size)) - 1, self.bin_count_std - bin_count/(total_pt/self.obs_size))])
        self.chi2 = chi2
    def writeout(self):
        #write chi2 to csv file
        dp = pd.DataFrame(data=self.chi2,columns=['age','dm','red','chi2'])
        path = "/data/outchi2/chi2_mc{}_sn{}".format(self.mc_num,self.sn)
        dp.to_csv(path)

    def __init__(self, mc_num, sn):
        self.mc_num = mc_num
        self.sn = sn
        self.read_input()
        self.search_vorbin()        
        self.writeout()
        print("done mc{}sn{}".format(self.mc_num,self.sn))
