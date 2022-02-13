#This code calculate chi square value from random sample from observation data
#Take start and end as input
import numpy as np
import pandas as pd
import sys
import os
import voronoi_2d_binning

class chi2:

	def read_input(self):
		#read M92 observed data
		self.obs_data = pd.read_csv('C:\\Users\\irisd\\Desktop\\share\\M92_fits_with_err')
		self.len_obs = len(self.obs_data)

	#search for bin number for each data point
	def search_point_location_bc(self, x, y, xBar, yBar):
		lenx = len(x)
		delta_x = xBar - x.reshape(lenx,1)
		delta_y = yBar - y.reshape(lenx,1)
		distance = np.square(delta_x) + np.square(delta_y)
		bin_num = np.argmin(distance, axis = 1)
		return bin_num
		
	def writevorbin(self, xBar, yBar, bin_count_std):
		#save the vorbin information
		bin_loc = np.vstack((xBar,yBar,bin_count_std)).T
		dp = pd.DataFrame(data=bin_loc, columns = ['xBar', 'yBar','bin_count_std'])
		path = "C:\\Users\\irisd\\Desktop\\share\\result\\vorbin_info"
		dp.to_csv(path)

	def resample(self, i):
		sample_list = np.random.randint(0,self.len_obs,size=self.sample_pt)
		Vvega = np.random.normal(self.obs_data['vv'].values[sample_list], self.obs_data['vv_err'].values[sample_list])
		#Ivega = np.random.normal(self.obs_data['ii'].values[sample_list], self.obs_data['ii_err'].values[sample_list])
		VIvega = np.random.normal(self.obs_data['vi'].values[sample_list], self.obs_data['vi_err'].values[sample_list])
		#Vvega = self.obs_data['vv']
		#Ivega = self.obs_data['ii']
		#VIvega = self.obs_data['vi']
		#data_resample = {'v':Vvega, 'i':Ivega, 'vi':VIvega}
		data_resample = {'v':Vvega, 'vi':VIvega}
		dp = pd.DataFrame(data=data_resample)
		df_resample = dp[(dp['vi'] < (self.obs_vi_max)) & (dp['vi'] > (self.obs_vi_min))& (dp['v'] < (self.obs_v_max)) & (dp['v'] > (self.obs_v_min))]
		self.total_pt = len(df_resample)
		path = "C:\\Users\\irisd\\Desktop\\share\\resample\\resample_{}".format(i)
		df_resample.to_csv(path,index=False)
		return df_resample
	
	def search_vorbin(self):
		Tb_size = 1000
		vorbin_sample_size = 180000
		chi2 = []
		for k in range(self.start, self.end):
			print("Starting {}th resample".format(k))
			dp = self.resample(k)
			#y = dp['v'].values[:self.len_obs]
			#x = dp['vi'].values[:self.len_obs] * 12.5
			#signal = np.array([1]*self.len_obs)
			#noise = np.array([1]*self.len_obs)
			#targetSN = 5
			y = dp['v'].values[:vorbin_sample_size]
			x = dp['vi'].values[:vorbin_sample_size] * 12.5
			signal = np.array([1]*vorbin_sample_size)
			noise = np.array([1]*vorbin_sample_size)
			targetSN = 15
			binNum, xNode, yNode, xBar, yBar, sN, nPixels, scale = voronoi_2d_binning.voronoi_2d_binning(x, y, signal, noise, targetSN, plot=0, quiet=1, pixelsize=1)
			#find standard bin count by search through all the theoretical data points
			bin_count_std = np.zeros(len(xBar))
			n_div = self.total_pt // Tb_size
			for i in range(n_div):
				bin_num = self.search_point_location_bc((dp['vi'].values[i*Tb_size:(i+1)*Tb_size])*12.5, dp['v'].values[i*Tb_size:(i+1)*Tb_size], xBar, yBar)
				for j in range(Tb_size):
					bin_count_std[bin_num[j]] += 1
			#do the last bit
			bin_num = self.search_point_location_bc((dp['vi'].values[n_div*Tb_size:])*12.5, dp['v'].values[n_div*Tb_size:], xBar, yBar)
			for j in range(self.total_pt - n_div*Tb_size):
				bin_count_std[bin_num[j]] += 1
			#to avoid divde by 0
			for i in range(len(bin_count_std)):
				if bin_count_std[i] == 0:
					bin_count_std[i] += 1
			#self.writevorbin(xBar, yBar, bin_count_std)
			bin_num = self.search_point_location_bc((self.obs_data['vi'].values)*12.5, self.obs_data['vv'].values, xBar, yBar)
			bin_count = np.zeros(len(xBar))
			for j in range(self.len_obs):
				bin_count[bin_num[j]] += 1
			#calculate chi2
			chi2.append([i, np.inner(np.divide(bin_count,bin_count_std/(self.total_pt/self.len_obs)) - 1, bin_count - bin_count_std/(self.total_pt/self.len_obs))])
		self.chi2 = chi2

	def writeout(self):
		#write chi2 to csv file
		dp = pd.DataFrame(data=self.chi2,columns=['i','chi2'])
		path = "C:\\Users\\irisd\\Desktop\\share\\result\\resample_chi2"
		dp.to_csv(path)


	def __init__(self, start, end):
		self.obs_vi_max = 0.721000016
		self.obs_vi_min = 0.351000011
		self.obs_v_max = 19.9249992
		self.obs_v_min = 15.9259996
		self.sample_pt = 2000000
		self.start = start
		self.end = end
		self.read_input()
		self.search_vorbin()        
		self.writeout()