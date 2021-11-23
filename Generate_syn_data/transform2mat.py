import os
import numpy as np
import scipy.io as io

for root,dirs,files in os.walk(r"./input"):
	for file in files:
		print(file)
		data = np.loadtxt(os.path.join(root,file))
		#data = np.transpose(data)
		io.savemat("./input_mat/"+file+".mat",{file:data})
		print("save "+file+" success!")
