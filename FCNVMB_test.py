# -*- coding: utf-8 -*-
"""

"""
################################################
########        IMPORT LIBARIES         ########
################################################

from ParamConfig import *
from PathConfig import *
from LibConfig import *

################################################
########         LOAD    NETWORK        ########
################################################

# Here indicating the GPU you want to use. if you don't have GPU, just leave it.
cuda_available = torch.cuda.is_available()
device         = torch.device("cuda" if cuda_available else "cpu")
model_file     = models_dir+modelname+'_epoch'+str(Epochs)+'.pkl'
net            = UnetModel(n_classes=Nclasses,in_channels=Inchannels, \
                           is_deconv=True,is_batchnorm=True,data_dim=DataDim,label_dim=ModelDim) 
net.load_state_dict(torch.load(model_file))
if torch.cuda.is_available():
    net.cuda()

################################################
########    LOADING TESTING DATA       ########
################################################
print('***************** Loading Testing DataSet *****************')

test_set,label_set,data_dsp_dim,label_dsp_dim = DataLoad_Test(test_size=TestSize,test_data_dir=test_data_dir, \
                                                              data_dim=DataDim,in_channels=Inchannels, \
                                                              model_dim=ModelDim,data_dsp_blk=data_dsp_blk, \
                                                              label_dsp_blk=label_dsp_blk,start=1, \
                                                              datafilename=datafilename,dataname=dataname, \
                                                              truthfilename=truthfilename,truthname=truthname)

test        = data_utils.TensorDataset(torch.from_numpy(test_set),torch.from_numpy(label_set))
test_loader = data_utils.DataLoader(test,batch_size=TestBatchSize,shuffle=False)


################################################
########            TESTING             ########
################################################

print() 
print('*******************************************') 
print('*******************************************') 
print('            START TESTING                  ') 
print('*******************************************') 
print('*******************************************') 
print()

# Initialization
since      = time.time()
Prediction = np.zeros((TestSize,label_dsp_dim[0],label_dsp_dim[1]),dtype=float)
GT         = np.zeros((TestSize,label_dsp_dim[0],label_dsp_dim[1]),dtype=float)
total      = 0
for i, (images,labels) in enumerate(test_loader):        
    images = images.view(TestBatchSize,Inchannels,data_dsp_dim[0],data_dsp_dim[1])
    labels = labels.view(TestBatchSize,Nclasses,label_dsp_dim[0],label_dsp_dim[1])
    images = images.to(device)
    labels = labels.to(device)
    
    # Predictions
    net.eval() 
    outputs  = net(images,label_dsp_dim)
    outputs  = outputs.view(TestBatchSize,label_dsp_dim[0],label_dsp_dim[1])
    outputs  = outputs.data.cpu().numpy()
    gts      = labels.data.cpu().numpy()
    
    # Calculate the PSNR, SSIM
    for k in range(TestBatchSize):
        pd   = outputs[k,:,:].reshape(label_dsp_dim[0],label_dsp_dim[1])
        gt   = gts[k,:,:].reshape(label_dsp_dim[0],label_dsp_dim[1])
        pd   = turn(pd)
        gt   = turn(gt)
        Prediction[i*TestBatchSize+k,:,:] = pd
        GT[i*TestBatchSize+k,:,:] = gt
        total = total + 1

# Save Results
SaveTestResults(Prediction,GT,results_dir)
        
# Record the consuming time
time_elapsed = time.time() - since
print('Testing complete in  {:.0f}m {:.0f}s' .format(time_elapsed // 60, time_elapsed % 60))


