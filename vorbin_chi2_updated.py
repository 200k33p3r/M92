#This code requires two inputs: mc_num and s/n
#This code outputs a list of chi2 for each age, dm and reddening combination
import numpy as np
import pandas as pd
import sys
import os
import voronoi_2d_binning

class chi2:

    def read_input(self):
        #read M92 observed data
        self.obs_data = pd.read_csv('/home/mying/Desktop/ObservationalData/M92_fitstars.dat',sep='\s+',names=['v','i','vi'],skiprows=3)
        #self.obs_size = len(self.obs_data)

    #search for bin number for each data point
    def search_point_location_bc(self,x, y, xBar, yBar):
        lenx = len(x)
        delta_x = xBar - x.reshape(lenx,1)
        delta_y = yBar - y.reshape(lenx,1)
        distance = np.square(delta_x) + np.square(delta_y)
        bin_num = np.argmin(distance, axis = 1)
        return bin_num
        
    def writevorbin(self, xBar, yBar, bin_count_std, age):
        #save the vorbin information
        bin_loc = np.vstack((xBar,yBar,bin_count_std)).T
        dp = pd.DataFrame(data=bin_loc, columns = ['xBar', 'yBar','bin_count_std'])
        path = "/data/outchi2_1/bin_mc{}_sn{}.age{}".format(self.mc_num,self.sn,age)
        dp.to_csv(path)
        
    def search_vorbin(self):
        #go through the search process
        Tb_size = 1000
        obs_vi_max = 0.721000016
        obs_vi_min = 0.351000011
        obs_v_max = 19.9249992
        obs_v_min = 15.9259996
        dm_max = 14.82
        dm_min = 14.62
        red_max = 0.05
        red_min = 0.0
        vorbin_pt = 180000
        dms = np.linspace(dm_min,dm_max,21)
        reds = np.linspace(red_min,red_max,6)
        ages = ['09000','10000','11000','12000','13000','14000','15000','16000']
        chi2 = []
        for age in ages:
            #read iso files
            dp = pd.read_csv("/data/outcmd/mc{}.a{}".format(self.mc_num,age),sep='\s+',names=['vi','v','i'],skiprows=3)
            #filter out data points that is out of boundary
            dp = dp[(dp['vi'] < (obs_vi_max - red_max)) & (dp['vi'] > (obs_vi_min - red_min))& (dp['v'] < (obs_v_max - dm_max)) & (dp['v'] > (obs_v_min - dm_min))]
            #generate vorbin use the first 180000 data points
            y = dp['v'].values[:vorbin_pt]
            x = dp['vi'].values[:vorbin_pt] * 12.5
            signal = np.array([1]*vorbin_pt)
            noise = np.array([1]*vorbin_pt)
            targetSN = self.sn
            binNum, xNode, yNode, xBar, yBar, sN, nPixels, scale = voronoi_2d_binning.voronoi_2d_binning(x, y, signal, noise, targetSN, plot=0, quiet=1, pixelsize=1)
            #find standard bin count by search through all the theoretical data points
            bin_count_std = np.zeros(len(xBar))
            total_pt = len(dp)
            n_div = total_pt // Tb_size
            for i in range(n_div):
                bin_num = self.search_point_location_bc((dp['vi'].values[i*Tb_size:(i+1)*Tb_size])*12.5, dp['v'].values[i*Tb_size:(i+1)*Tb_size], xBar, yBar)
                for j in range(Tb_size):
                    bin_count_std[bin_num[j]] += 1
            #do the last bit
            bin_num = self.search_point_location_bc((dp['vi'].values[n_div*Tb_size:])*12.5, dp['v'].values[n_div*Tb_size:], xBar, yBar)
            for j in range(total_pt - n_div*Tb_size):
                bin_count_std[bin_num[j]] += 1
            #to avoid divde by 0
            for i in range(len(bin_count_std)):
                if bin_count_std[i] == 0:
                    bin_count_std[i] += 1
            self.writevorbin(xBar, yBar, bin_count_std, age)
            #search through observed data
            for dm in dms:
                for red in reds:
                    bin_count = np.zeros(len(xBar))
                    dp = self.obs_data
                    dp = dp[(dp['vi'] - red < (obs_vi_max - red_max)) & (dp['vi'] - red > (obs_vi_min - red_min))& (dp['v'] - dm < (obs_v_max - dm_max)) & (dp['v'] - dm > (obs_v_min - dm_min))]
                    obs_size = len(dp)
                    bin_num = self.search_point_location_bc((dp['vi'].values - red)*12.5, dp['v'].values - dm, xBar, yBar)
                    for j in range(obs_size):
                        bin_count[bin_num[j]] += 1
                    #calculate chi2
                    chi2.append([age, dm, red, np.inner(np.divide(bin_count,bin_count_std/(total_pt/obs_size)) - 1, bin_count - bin_count_std/(total_pt/obs_size))])
        self.chi2 = chi2


    def writeout(self):
        #write chi2 to csv file
        dp = pd.DataFrame(data=self.chi2,columns=['age','dm','red','chi2'])
        path = "/data/outchi2_1/chi2_mc{}_sn{}".format(self.mc_num,self.sn)
        dp.to_csv(path)


    def __init__(self, mc_num, sn):
        self.mc_num = mc_num
        self.sn = sn
        self.read_input()
        self.search_vorbin()        
        self.writeout()
        print("done mc{}sn{}".format(self.mc_num,self.sn))
