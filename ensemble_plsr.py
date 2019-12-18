############################## ALAB PLSR Adaptation ##########################
# The purpose of this code is to take a generically formatted csv database
# of foliar chemistries and spectra, and run PLSR calibrations on those data.  While multiple
# methods of running this analysis are viable using this code (through use of
# the associated settings file), this is adapted directly from Chadwick, K.D. and G.P. Asner. 
# "Organismic-scale remote sensing of canopy foliar traits in lowland tropical forests." 
# Remote Sensing 8.2 (2016):87. The original code used the autopls package, was slow to execute.
# This updated code adapts large portions of that work in a more clearly documented, easier to modify,
# and faster way.  The paper for the autopls package is Schmidtlein, S., H. Feilhauer, and H. Bruelheide.
# "Mapping plant strategy types using remote sensing." Journal of Vegetation Science 23.3 (2012):395-405.
# These adaptations were originally performed by Phil Brodrick under CAO research contract to The 
# Rainforest Trust (CIW #10719). Contact P. Brodrick (pbrodrick@ciw.edu) for code questions. Contact G.P.
# Asner (gpa@ciw.edu) or R.E. Martin (rmartin@ciw.edu) for applications and programmatic questions.
#############################################################################

import pandas as pd
import numpy as np
import numpy.matlib

from sklearn import linear_model
from scipy import stats
from sklearn import preprocessing
from sklearn import ensemble
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
from sklearn.cross_decomposition import PLSRegression
import read_settings_file
from sklearn.model_selection import KFold
import multiprocessing 
import time
from functools import partial
import subprocess
import matplotlib.gridspec as gridspec
import sys, os

# internal function for autoselect_plsr - runs for a particular number of principles/components
def internal_autoselect_iteration(X,Y,ncomp):
    comp = ncomp
    kf = KFold(n_splits=10,shuffle=True,random_state=10)
    kf.get_n_splits(X)
    loc_score = []
    loc_slope = []
    loc_inter = []
    for train_index, test_index in kf.split(X):
      model = PLSRegression(n_components=ncomp,scale=False,tol=1e-9)
      model.fit(X[train_index,:],Y[train_index])
      pred_y = model.predict(X[test_index,:])
      loc_score.append(np.sqrt(np.nanmean(np.power(pred_y-Y[test_index],2))))

    return np.nanmean(loc_score)

# make the autopls coefficient selection without attempting to drop any additional bands out (similar to A1 in autopls)
def autoselect_plsr(X,Y,max_comp):
  complist = np.arange(2,max_comp)[::-1]
  scores = map(partial(internal_autoselect_iteration,X,Y), complist.tolist())

  ncomp =  complist[np.argmin(scores)]

  model = PLSRegression(n_components=ncomp,scale=False,tol=1e-9)
  model.fit(X,Y)
  loc_coeff = np.zeros(X.shape[1]+1)
  loc_coeff[1:] = np.squeeze(model.coef_)
  loc_coeff[0] = model.y_mean_ - np.dot(model.x_mean_ , model.coef_)
  return loc_coeff, list(scores)


# apply a forward or backward transformation, intended for chem data if needed
def apply_transform(Y,transform,invert=False):
  if (transform == 'None' or transform == None):
    return Y
  if (invert == False):
    if (transform == 'log'):
      return np.log(Y)
    if (transform == 'sqrt'):
      return np.sqrt(Y)
    if (transform == 'square'):
      return np.power(Y,2)
  else:
    if (transform == 'log'):
      return np.exp(Y)
    if (transform == 'sqrt'):
      return np.power(Y,2)
    if (transform == 'square'):
      return np.sqrt(Y)

# perform brightness normalization on either a matrix or a vector of reflectance coefficients
def brightness_normalize(X,return_norm=False):
  if (len(X.shape) == 1):
    norm = np.sqrt(np.sum(np.power(X,2)))
    res = X.astype(float) / norm
  else:
    norm = np.sqrt(np.sum(np.power(X,2),axis=1))
    res = X.astype(float) / np.transpose(np.matlib.repmat(norm,X.shape[1],1))
  if (return_norm == True):
    return res,norm
  else:
    return res

# print active performance results
def print_res(x,y,label,print_xy=False):
  print('######### {}  #########'.format(label))
  print('slope','intercept','r squared','rmse')
  slope, intercept, r_value, p_value, std_err = stats.linregress(np.squeeze(x),np.squeeze(y))
  print('{}, {}, {}, {}'.format(slope,intercept,r_value**2,np.sqrt(np.mean(np.power(x-y,2)))))
  if (print_xy == True):
    print(x)
    print(y)
 
# basic plsr prediction assuming that the coefficient vector has the intercept as the first element,
# and coefficients matching the reflectance matrix as the rest of the elements
def plsr_predict(in_coeff,x_mat,invert_scaling=False):
  pred = np.sum(np.matlib.repmat(in_coeff[1:],x_mat.shape[0],1) * x_mat,axis=1) + in_coeff[0]
  if (invert_scaling == True):
    if (sf.get_setting('scale response')):
      pred = y_scaler.inverse_transform(pred)
    pred = apply_transform(pred,sf.get_setting('chem transforms')[_chem],invert=True)
  return pred


# cleanup and perform specified pruning of the input data.
# listin is a list of all of the objects to get filtered.  The first element of the list should
# be the spectra ONLY.
def cleanup_spectra(listin):
  x = listin[0]
  x, norm = brightness_normalize(x,return_norm=True)
  good_row = np.ones(norm.shape).astype(bool)
  good_row[np.sum(np.isnan(listin[0]),axis=1) > 0] = False
  good_row[np.isnan(listin[1]) ] = False
  if (sf.get_setting('brightness minimum') != -1):
    good_row[norm < sf.get_setting('brightness minimum')] = False
  if (sf.get_setting('brightness maximum') != -1):
    good_row[norm > sf.get_setting('brightness maximum')] = False
  if (sf.get_setting('brightness normalize') == True):
    listin[0] = x

  # ndvi bands (46 - 34) / (46 + 34) [1 - based]
  if (sf.get_setting('ndvi minimum') != -1 or sf.get_setting('ndvi maximum') != -1):
    bad_bands = sf.get_setting('bad bands')
    ndvi_bands = sf.get_setting('ndvi bands')

    if (ndvi_bands[0] in bad_bands or ndvi_bands[1] in bad_bands):
      print('cannot use NDVI min without bands {}'.format(ndvi_bands))
    ndvi_bands = ndvi_bands -1
    ndvi = (x[:,ndvi_bands[1]] - x[:,ndvi_bands[0]]) / (x[:,ndvi_bands[1]] + x[:,ndvi_bands[0]])
    if (sf.get_setting('ndvi minimum') != -1):
      good_row[ndvi < sf.get_setting('ndvi minimum')] = False
    if (sf.get_setting('ndvi maximum') != -1):
      good_row[ndvi > sf.get_setting('ndvi maximum')] = False
  
 ### What does this do? 
  for n in range(0,len(listin)):
    if (len(listin[n].shape) == 1):
      listin[n] = listin[n][good_row]
    elif (len(listin[n].shape) == 2):
      listin[n] = listin[n][good_row,:]
    else:
      print('bad listin length',n)
      quit()

  return listin
    
  
########################################### Begin Main Code ########################################


np.random.seed(13)
sf = read_settings_file.settings(sys.argv[1])
sf.clean_and_check()
df = pd.read_csv(sf.get_setting('csv file'),sep=',')

for col in sf.get_setting('ignore columns'):
    df.pop(col)

header =  list(df)
df = np.array(df)

if (os.path.isdir('figs') == False):
    subprocess.call('mkdir figs',shell=True)
if (os.path.isdir('stats') == False):
    subprocess.call('mkdir stats',shell=True)
if (os.path.isdir('plot_points') == False):
    subprocess.call('mkdir plot_points',shell=True)
if (os.path.isdir('coeff') == False):
    subprocess.call('mkdir coeff',shell=True)
if (os.path.isdir('lvs') == False):
    subprocess.call('mkdir lvs',shell=True)


################## Data frame prep - get it cleaned up before looping through chems #####
# get the chem list / chem transforms of interest
chemlist = sf.get_setting('chems')
chemtransform_list = sf.get_setting('chem transforms')

# use the settings file and input csv to find the columns of interest 
crown_col = -1
x_col = []
for i in range(0,len(header)):
  if (header[i] == sf.get_setting('crown col')):
    crown_col = i
  if (header[i] == sf.get_setting('test training col')):
    test_col = i
  if (sf.get_setting('band preface') in header[i]):
    lbn = int(header[i].replace(sf.get_setting('band preface'),''))
    if (lbn > 0 and lbn <= sf.get_setting('max band')):
      if (np.sum(sf.get_setting('bad bands') == lbn) == 0):
        x_col.append(i)
x_col = np.array(x_col)
if (crown_col == -1):
  print('no crown index \'',sf.get_setting('crown col'),'\' found, terminating.')
      
# remove cronws where there aren't at least min pixels per crown
un_crown = np.unique(df[:,crown_col])
good_rows = np.ones(len(df)).astype(bool)
for n in range(0,len(un_crown)):
  if (np.sum(df[:,crown_col] == un_crown[n]) < sf.get_setting('min pixel per crown')):
    good_rows[df[:,crown_col] == un_crown[n]] = False

prev_size = df.shape[0]
df = df[good_rows,:]
print('Number of pixels trimmed by min pixel per crown:',prev_size - df.shape[0])

     
# convert the crown columns to a list of unique IDs that are numerical (for convenience)
# in numpy datatype specific arrays
crown_names = df[:,crown_col].copy() 
un_crown_names = np.unique(df[:,crown_col]).copy()
for n in range(0,len(un_crown_names)):
  df[df[:,crown_col] == un_crown_names[n],crown_col] = n 
crown_ids = np.array(df[:,crown_col]).copy()
un_crown_ids = np.unique(df[:,crown_col]).copy()


# Define test set
test_rows = df[:,test_col] == sf.get_setting('test set value') 

#Remove columns that have been filtered on already
df = np.delete(df,test_col,axis=1)
header.pop(test_col)
if (test_col < np.min(x_col)):
    x_col -= 1

# convert the dataframe to an array of floats
df = np.float64(df)

# separate out training and testing sets
test_df = df[test_rows]
test_crown_ids = crown_ids[test_rows].copy()
un_test_crown_ids = np.unique(test_crown_ids)
test_crown_names = crown_names[test_rows].copy()
un_test_crown_names = np.unique(test_crown_names)

df = df[np.logical_not(test_rows)]
complete_crown_names = crown_names[np.logical_not(test_rows)].copy()

#complete_crown_names = crown_names[test_rows].copy()

lv_file = 'lvs/' + sf.get_setting('version name') + '.txt'
outstr = 'Chem,# LV Selected'
for n in range(2,sf.get_setting('max components')):
  outstr += ',LV-' + str(n) + 'RMSE'
subprocess.call('echo \"' + outstr + '\" > ' + lv_file,shell=True)


stat_file = 'stats/' + sf.get_setting('version name') + '.txt'
outstr = 'Chem, Val-Slope, Val-Intercept, Val-R2, Val-RMSE, Val-%RMSE, Cal-Slope, Cal-Intercept, Cal-R2, Cal-RMSE, Cal-%RMSE'
subprocess.call('echo \"' + outstr + '\" > ' + stat_file,shell=True)

print('preparing data...')

################# step through each chem, run all PLSR code.  Each chem is independent ############
for _chem in range(0,len(chemlist)):
  print(chemlist[_chem])
 
  # find the column of the dataframe corresponding to the chem of interest
  resp_col = -1
  for i in range(0,len(header)):
    if (header[i] == chemlist[_chem]):
      resp_col = i
      break

  if (resp_col == -1):
    print('no response column\'',chemlist[_chem],'\' found')
    quit()
  
  # convert the dataframe to the components we'll use for this chem.  Apply any normalizations
  # and scalings specified in the settings file
  crown = df[:,crown_col]
  X = df[:,x_col]
  Y = df[:,resp_col]
  
  # eliminate and data with bad response columns
  good_rows = np.logical_and(np.isnan(Y) == False,np.isinf(Y) == False)
  good_rows[np.all(X == 0,axis=1)] = False
  good_rows[Y == -9999] = False
  X = X[good_rows,:]
  Y = Y[good_rows]
  crown = crown[good_rows]
  crown_names = complete_crown_names[good_rows]

  Y = apply_transform(Y,sf.get_setting('chem transforms')[_chem])
  [X,Y,crown,crown_names] = cleanup_spectra([X,Y,crown,crown_names])

  # do a check on X and convert any nans or infs to 0s, which the PLSR will ignore
  X[np.isnan(X)] = 0
  X[np.isinf(X)] = 0

  # if requested in the settings file, do a check to see which bands are all 0.  The program
  # terminates after doing this check, assuming the user intends to put these back into
  # the settings file
  if (sf.get_setting('find bad bands')):
    print('bad band guesses: ')
    lh = np.array(header)[x_col]
    for t in range(0,X.shape[1]):
      if np.all(X[:,t] == 0):
        print(lh[t].split('B')[1])
    print('')
    quit() 
  
  # set up a 'global holdout' - this is a true test set, not a validation set.
  # perform the holdout on a crown-level basis, and set up the crown averaging ahead of time
  global_holdout = np.ones(len(X))
  un_crown = np.unique(crown)
  un_crown = un_crown[np.random.permutation(len(un_crown))]
  ind = 0
  for m in range(0,int(len(un_crown)*sf.get_setting('test set holdout fraction'))):
    global_holdout[crown == un_crown[m]] = 0
  
  # define test sets
  global_X_test = test_df[:,x_col]
  global_Y_test = test_df[:,resp_col]
  global_crown_test = test_df[:,crown_col]
  global_crown_name_test = test_crown_names.copy()

  good_row = np.logical_not(np.all(global_X_test == 0,axis=1))
  good_row[global_Y_test == -9999] = False

  global_X_test = global_X_test[good_row]
  global_Y_test = global_Y_test[good_row]
  global_crown_test = global_crown_test[good_row]
  global_crown_name_test = global_crown_name_test[good_row]

  [global_X_test,global_Y_test,global_crown_test,global_crown_name_test] = cleanup_spectra([global_X_test,global_Y_test,global_crown_test,global_crown_name_test])
  
  global_crown_test_un,_gctu = np.unique(global_crown_test,return_index=True)
  global_crown_Y_test = []
  for m in range(0,len(global_crown_test_un)):
    global_crown_Y_test.append(np.mean(global_Y_test[global_crown_test == global_crown_test_un[m]]))
  global_crown_Y_test = np.array(global_crown_Y_test)
  global_crown_name_test = global_crown_name_test[_gctu]

  # filter down our X,Y, and crown sets to NOT include the test set data
  X = X[global_holdout != 0,:]
  Y = Y[global_holdout != 0].reshape(-1,1)
  crown = crown[global_holdout != 0]
  crown_names = crown_names[global_holdout != 0]
  

  # fit the x and y scalers (application will come later if it's supposed to)
  x_scaler = preprocessing.StandardScaler()
  x_scaler.fit(X)
  y_scaler = preprocessing.StandardScaler()
  y_scaler.fit(Y.reshape(-1,1))

  # scale the X if necessary (we won't scale the response, because we'll back-transform the predicted
  # y's in place
  if (sf.get_setting('scale features')):
    global_X_test_scaled = x_scaler.transform(global_X_test)
  else:
    global_X_test_scaled = global_X_test
  

  # write training and testing data to an output file for external use (r call, or verification)
  train_file = 'munged/train_dat_chem' + chemlist[_chem] + '.csv'
  test_file = 'munged/test_dat_chem' + chemlist[_chem] + '.csv'

    
  
  ##################### start the iteration loop ########################
  # at each iteration, we will pull out a random sample of pixels from different crowns,
  # and then fit the PLSR to that data set.  As we go, we'll update the overall coefficients
  # and check performance of the local PLSR, along with the cumulative, developing model
  all_coeff = np.zeros((sf.get_setting('iterations'),X.shape[1]+1))
  coeff_perf = np.ones(sf.get_setting('iterations'))*-9999
  culm_lv = np.ones(sf.get_setting('iterations'))*-9999
  coeff = np.zeros(X.shape[1]+1)
  for iter in range(0,sf.get_setting('iterations')):

    un_crown = np.unique(crown)
    un_crown = un_crown[np.random.permutation(len(un_crown))]
    holdout = np.ones(len(crown))
    for n in range(0,int(len(un_crown)*sf.get_setting('iteration holdout fraction'))):
      holdout[crown == un_crown[n]] = 0
    
    
    train = holdout != 0
    validation = holdout == 0
    

    # subsample a specified number of samples per crown for the training data
    train_crown_un = np.unique(crown[train])
    train_rows = np.zeros(len(X)).astype(bool)
    if (sf.get_setting('samples per crown') == -1):
      for m in range(0,len(train_crown_un)):
        train_rows[crown == train_crown_un[m]] = True
    else:
     for m in range(0,len(train_crown_un)):
        train_rows[np.random.choice(np.arange(len(X))[crown == train_crown_un[m]],\
                                    size=sf.get_setting('samples per crown'),\
                                    replace=True)] = True

    X_train = X[train_rows,:]
    Y_train = Y[train_rows]
    

    ###### scale and/or transform as necessary ############
    if (sf.get_setting('scale features')):
      X_train = x_scaler.transform(X_train)
    if (sf.get_setting('scale response')):
      Y_train = y_scaler.transform(Y_train)


    if (np.sum(validation) > 0):
      validation_crown_un = np.unique(crown[validation])
      X_validation = X[validation,:]
      Y_validation = Y[validation]

      if (sf.get_setting('scale features')):
        X_validation = x_scaler.transform(X_validation)
      if (sf.get_setting('scale response')):
        Y_validation = y_scaler.transform(Y_validation)
    
    
    ################# initialize and train model
    iter_coeff, score = autoselect_plsr(X_train,Y_train,sf.get_setting('max components'))
    outstr = chemlist[_chem] + ',' + str(np.argmin(score)+1)
    for nnn in range(0,len(score)):
      outstr += ',' + str(round(score[nnn],4))
    subprocess.call('echo \"' + outstr + '\" >> ' + lv_file,shell=True)

   ############# make local predictions #####################
    if (np.sum(validation) > 0):
      loc_pred_y = np.squeeze(plsr_predict(iter_coeff,X_validation,invert_scaling=True))
      pred_y = loc_pred_y.reshape(-1,1)
      act_y = Y[validation]
   
      ##### find crown averages
      crown_pred_y = []
      crown_y = []
      for m in range(0,len(validation_crown_un)):
        loc_crown = crown[validation] == validation_crown_un[m]
        crown_pred_y.append(np.mean(loc_pred_y[loc_crown]))
        crown_y.append(np.mean(Y[validation][loc_crown]))
   
      crown_pred_y = np.array(crown_pred_y)
      crown_y = np.array(crown_y)
      
    ################ update the cumulative coefficients (and stored references)
    culm_lv[iter] = np.argmin(score)+1
    if (np.sum(validation) > 0):
      coeff_perf[iter] = np.sqrt(np.nanmean(np.power(pred_y-act_y,2)))
      all_coeff[iter,:] = iter_coeff.copy()

      use_coeff = np.zeros(len(coeff_perf)).astype(bool)
      use_coeff[coeff_perf <= np.percentile(coeff_perf[coeff_perf >= 0], (sf.get_setting('iteration fraction used'))*100)] = True
      use_coeff[coeff_perf == -9999] = False

      lv = np.nanmean(culm_lv[use_coeff])
      coeff = np.nanmean(all_coeff[use_coeff,:],axis=0)
    else:
      coeff_perf[iter] = 1
      all_coeff[iter,:] = iter_coeff.copy()
      use_coeff = np.zeros(len(coeff_perf)).astype(bool)
      use_coeff[coeff_perf >= 0] = True
      lv = np.nanmean(culm_lv[use_coeff])
      coeff = np.nanmean(all_coeff[use_coeff,:],axis=0)
 

    ################# make global test set predictions for current iteration ##############
    global_pred_y = plsr_predict(iter_coeff,global_X_test_scaled,invert_scaling=True) 

    ##### find crown averages
    global_crown_pred_y = []
    for m in range(0,len(global_crown_test_un)):
      global_crown_pred_y.append(np.mean(global_pred_y[global_crown_test == global_crown_test_un[m]]))
      
    global_crown_pred_y = np.array(global_crown_pred_y)
  
    print_res(global_crown_pred_y,global_crown_Y_test,'global holdout per crown pred_y-y (iteration)')
  
    ################# make global test set predictions for cumulative model ##############
    global_pred_y = plsr_predict(coeff,global_X_test_scaled,invert_scaling=True) 

    global_crown_pred_y = []
    for m in range(0,len(global_crown_test_un)):
      global_crown_pred_y.append(np.mean(global_pred_y[global_crown_test == global_crown_test_un[m]]))
    global_crown_pred_y = np.array(global_crown_pred_y)
  
    print_res(global_crown_pred_y,global_crown_Y_test,'global holdout per crown pred_y-y (cumulative model)')
    print_res(global_Y_test,global_pred_y,'global holdout per pix y-pred_y (cumulative model)')

    ################# make training set predictions #################
    loc_X = X.copy()
    loc_Y = Y.copy()
    # undo any transform done on the Y space
    loc_Y = apply_transform(loc_Y,sf.get_setting('chem transforms')[_chem],invert=True)
    if (sf.get_setting('scale features')):
      loc_X = x_scaler.transform(loc_X)

    tpy = plsr_predict(coeff,loc_X,invert_scaling=True) 
    train_crown_pred_y = []
    train_crown_y = []
    un_crown,_un_crown_ind = np.unique(crown,return_index=True)
    un_crown_names = crown_names[_un_crown_ind]
    for m in range(0,len(un_crown)):
      train_crown_pred_y.append(np.mean(tpy[crown == un_crown[m]]))
      train_crown_y.append(np.mean(loc_Y[crown == un_crown[m]]))
    train_crown_pred_y = np.array(train_crown_pred_y)
    train_crown_y = np.array(train_crown_y)
  
    print_res(train_crown_y,train_crown_pred_y,'train set y-pred_y (cumulative model)')

  
    ################### plot current progress (for cumulative model) ############### 
    ub = np.max(np.append(global_Y_test,global_pred_y))
    lb = np.min(np.append(global_Y_test,global_pred_y))

    fig = plt.figure(facecolor='white',figsize=(12,6))
    gs1 = gridspec.GridSpec(1, 2)
    plt.subplot(gs1[0])
    plt.scatter(train_crown_y,train_crown_pred_y,color='grey',s=3)
    plt.scatter(global_crown_Y_test,global_crown_pred_y,color='red',s=5)
    plt.plot([lb,ub],[lb,ub],color='black',ls='--')
    plt.xlim([lb,ub])
    plt.ylim([lb,ub])
    plt.xlabel('Observed Data')
    plt.ylabel('Model Prediction')

    plt.subplot(gs1[1])
    plt.scatter(train_crown_pred_y,train_crown_y,color='grey',s=3)
    plt.scatter(global_crown_pred_y,global_crown_Y_test,color='red',s=5)
    plt.plot([lb,ub],[lb,ub],color='black',ls='--')
    plt.xlim([lb,ub])
    plt.ylim([lb,ub])
    plt.xlabel('Model Prediction')
    plt.ylabel('Observed Data')

    fig.savefig('figs/' + sf.get_setting('version name') + '_validation_' + chemlist[_chem] + '.png',dpi=100,format='png',facecolor='white',bbox_inches='tight')
    plt.close()

    print_cyt = np.append(global_crown_Y_test,train_crown_y)
    print_cytp = np.append(global_crown_pred_y,train_crown_pred_y)
    print_cn = np.append(global_crown_name_test,un_crown_names)
    print_cn_tt = np.append(np.array(['Validation' for x in range(0,len(global_crown_Y_test))]),np.array(['Calibration' for x in range(0,len(train_crown_y))]))
    out_df = pd.DataFrame(data=np.transpose(np.vstack([print_cn,print_cyt,print_cytp,print_cn_tt])))
    out_df.to_csv('plot_points/' + sf.get_setting('version name') + '_' + chemlist[_chem] + '.csv',sep=',',index=False)


    ################### plot current coefficients (for cumulative model) ############### 
    fig = plt.figure(facecolor='white',figsize=(10,6))
    wavelengths = np.arange(0,X.shape[1])*sf.get_setting('wavelength interval') + sf.get_setting('lowest wavelength')
    if (np.sum(coeff_perf >= 0) > 1):
      plt.fill_between(wavelengths,coeff[1:]-np.std(all_coeff[use_coeff,1:],axis=0),coeff[1:]+np.std(all_coeff[use_coeff,1:],axis=0),facecolor='blue',alpha=0.6)
    plt.plot(wavelengths,coeff[1:],color='black')
    plt.xlabel('Wavelength [nm]')
    plt.ylabel('PLSR Coeff')
    fig.savefig('figs/' + sf.get_setting('version name') + '_coeff_' + chemlist[_chem] + '.png',dpi=100,format='png',facecolor='white',bbox_inches='tight')
    plt.close()


  ###############################################################################################
  ############################## all iterations finished for single chem ########################
  ###############################################################################################

  ############ Write coefficient file and std file #############
  coeff_file = 'coeff/coeff_values_' + sf.get_setting('version name') + '.csv'
  coeff_std_file = 'coeff/coeff_std_' + sf.get_setting('version name') + '.csv'
  coeff_stds=np.std(all_coeff[use_coeff,:],axis=0)
  if (_chem == 0): # make header on first chem, overwrite previous version
    outstr = 'Chem,Intercept'
    for n in range(1,215):
      outstr += ',B' + str(n)
    subprocess.call('echo \"' + outstr + '\" > ' + coeff_file,shell=True)
    subprocess.call('echo \"' + outstr + '\" > ' + coeff_std_file,shell=True)

  outstr = chemlist[_chem] + ',' + str(coeff[0])
  outstr_std = chemlist[_chem] + ',' + str(coeff_stds[0])
  it = 1
  for n in range(1,1+sf.get_setting('max band')):
    if (n in sf.get_setting('bad bands')):
      outstr += ',0'
      outstr_std += ',0'
    else:
      outstr += ',' + str(coeff[it])
      outstr_std += ',' + str(coeff_stds[it])
      it += 1
  subprocess.call('echo \"' + outstr + '\" >> ' + coeff_file,shell=True)
  subprocess.call('echo \"' + outstr_std + '\" >> ' + coeff_std_file,shell=True)

  ############ Update coefficient plot to include spacing of null bands #############
  fig = plt.figure(facecolor='white',figsize=(10,6))
  wavelengths = np.arange(0,sf.get_setting('max band'))*sf.get_setting('lowest wavelength') + sf.get_setting('wavelength interval')
  plot_coeff = []
  plot_coeff_std = []
  it = 1
  for n in range(1,1+sf.get_setting('max band')):
    if (n in sf.get_setting('bad bands')):
      plot_coeff.append(0)
      plot_coeff_std.append(0)
    else:
      plot_coeff.append(coeff[it])
      plot_coeff_std.append(coeff_stds[it])
      it+=1
  plot_coeff = np.array(plot_coeff)    
  plot_coeff_std = np.array(plot_coeff_std)    

  plt.fill_between(wavelengths,plot_coeff-plot_coeff_std,plot_coeff+plot_coeff_std,facecolor='blue',alpha=0.6)
  plt.plot(wavelengths,plot_coeff,color='black')
  plt.xlabel('Wavelength [nm]')
  plt.ylabel('PLSR Coeff')
  fig.savefig('figs/' + sf.get_setting('version name') + '_coeff_' + chemlist[_chem] + '.png',dpi=100,format='png',facecolor='white',bbox_inches='tight')
  plt.close()

  train_slope, train_intercept, train_r_value, p_value, std_err = stats.linregress(np.squeeze(train_crown_y),np.squeeze(train_crown_pred_y))
  train_rmse = np.sqrt(np.mean(np.power(train_crown_pred_y-train_crown_y,2)))
  train_rel_rmse = np.sqrt(np.mean(np.power((train_crown_pred_y-train_crown_y)/train_crown_y,2)))

  test_slope, test_intercept, test_r_value, p_value, std_err = stats.linregress(np.squeeze(global_crown_Y_test),np.squeeze(global_crown_pred_y))
  test_rmse = np.sqrt(np.mean(np.power(global_crown_pred_y-global_crown_Y_test,2)))
  test_rel_rmse = np.sqrt(np.mean(np.power((global_crown_pred_y-global_crown_Y_test)/global_crown_Y_test,2)))
  outstr = str(chemlist[_chem])+', '+str(test_slope)+', '+str(test_intercept)+', '+str(test_r_value**2)+', '+str(test_rmse)+', ' + str(test_rel_rmse) + ',' + str(train_slope) + ', ' + str(train_intercept) + ', ' + str(train_r_value**2) + ', ' + str(train_rmse) + ',' + str(train_rel_rmse)
  subprocess.call('echo \"' + outstr + '\" >> ' + stat_file,shell=True)






