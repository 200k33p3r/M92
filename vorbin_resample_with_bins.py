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
		self.obs_data = pd.read_csv('C:\\Users\\marti\\Desktop\\school work\\Dartmouth\\M92 project\\M92_fits_with_err_with_loc_with_bin')
		self.len_obs = len(self.obs_data)
		#read AS test error
		self.dps_Ierr = []
		self.dps_Verr = []
		for i in range(80):
			self.dps_Ierr.append(pd.read_csv("C:\\Users\\marti\\Desktop\\school work\\Dartmouth\\M92 project\\M92-20211208T215351Z-001\\M92\\SimulateCMD\\inputfiles\\Ierr{:02d}s.dat".format(i + 1),sep='\s+',skiprows=3,names=['Ierr']))
			self.dps_Verr.append(pd.read_csv("C:\\Users\\marti\\Desktop\\school work\\Dartmouth\\M92 project\\M92-20211208T215351Z-001\\M92\\SimulateCMD\\inputfiles\\Verr{:02d}s.dat".format(i + 1),sep='\s+',skiprows=3,names=['Verr']))

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
		path = "C:\\Users\\marti\\Desktop\\school work\\Dartmouth\\M92 project\\resample_with_bin\\result\\vorbin_info"
		dp.to_csv(path)

	def resample(self, j):
		sample_list = np.random.randint(0,self.len_obs,size=self.sample_pt)
		Ierr = np.zeros(self.sample_pt)
		Verr = np.zeros(self.sample_pt)
		for i in range(self.sample_pt):
			Ierr[i] = self.dps_Ierr[self.obs_data['Ibin'].values[sample_list[i]]]['Ierr'].values[np.random.randint(0, high=len(self.dps_Ierr[self.obs_data['Ibin'].values[sample_list[i]]]))]
			Verr[i] = self.dps_Verr[self.obs_data['Vbin'].values[sample_list[i]]]['Verr'].values[np.random.randint(0, high=len(self.dps_Verr[self.obs_data['Vbin'].values[sample_list[i]]]))]
		Vvega = self.obs_data['vv'].values[sample_list] - Verr
		Ivega = self.obs_data['ii'].values[sample_list] - Ierr
		VIvega = Vvega - Ivega
		#Vvega = self.obs_data['vv']
		#Ivega = self.obs_data['ii']
		#VIvega = self.obs_data['vi']
		#data_resample = {'v':Vvega, 'i':Ivega, 'vi':VIvega}
		data_resample = {'v':Vvega, 'vi':VIvega}
		dp = pd.DataFrame(data=data_resample)
		df_resample = dp[(dp['vi'] < (self.obs_vi_max)) & (dp['vi'] > (self.obs_vi_min))& (dp['v'] < (self.obs_v_max)) & (dp['v'] > (self.obs_v_min))]
		self.total_pt = len(df_resample)
		path = "C:\\Users\\marti\\Desktop\\school work\\Dartmouth\\M92 project\\resample_with_bin\\resample_{}".format(j)
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
			chi2.append([k, np.inner(np.divide(bin_count,bin_count_std/(self.total_pt/self.len_obs)) - 1, bin_count - bin_count_std/(self.total_pt/self.len_obs))])
		self.chi2 = chi2

	def writeout(self):
		#write chi2 to csv file
		dp = pd.DataFrame(data=self.chi2,columns=['i','chi2'])
		path = "C:\\Users\\marti\\Desktop\\school work\\Dartmouth\\M92 project\\resample_with_bin\\result\\resample_chi2_{}".format(self.start)
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