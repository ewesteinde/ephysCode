import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from sklearn.linear_model import LinearRegression
bData = pd.read_csv(r'Z:\Dropbox (HMS)\Wilson_Lab_Data\ephys\identified_PFL2\newPFL2\20221129\fly_1\cell_1\jumpCL_trial_6\processedData\bData_trial_1.csv')
tData = pd.read_csv(r'Z:\Dropbox (HMS)\Wilson_Lab_Data\ephys\identified_PFL2\newPFL2\20221129\fly_1\cell_1\jumpCL_trial_6\processedData\tData_trial_1.csv')
# print(bData.head(10))
# print(bData.shape)

# print(tData.head(10))
# print(tData.shape)

activity = tData.loc[:,'smoothVm'].to_numpy()
print(activity[0:10])
behaviour = bData.loc[:,'vel_for'].to_numpy()
print(behaviour[0:10])

# plt.scatter(behaviour,activity, s = 0.2)
# plt.xlabel('behaviour')
# plt.ylabel('activity')

count = 0
lags = np.arange(-500, 501, 10)
rsq_values = [0] * lags.shape[0]
lag_values = [0] * lags.shape[0]

for lag in lags:
    print(lag)
    if lag > 0 :
        b = behaviour[abs(lag):]
        a = activity[:-lag]
    elif lag < 0:
        b = behaviour[:lag]
        a = activity[abs(lag):]
    else:
        b = behaviour
        a = activity

    bnan = np.isnan(b)
    numbnan = sum(bnan)
   # print(numbnan)
    tnan = np.isnan(a)
    numtnan = sum(tnan)
    #print(numtnan)

    if numbnan != 0 or numtnan != 0:
        bidx = ~bnan
        tidx = ~tnan       
        b = b[bidx & tidx]
        a = a[bidx & tidx]

    # fig, ax = plt.subplots(figsize = [9,7])
    # ax.plot(a, color = 'b')
    # ax2 = ax.twinx()
    # ax2.plot(b,color = 'r')
    # plt.show()


    X = b.reshape(-1, 1)
    #print(X.shape)
    y = a
    #print(y.shape)

    model = LinearRegression().fit(X, y)

    r_sq = model.score(X, y)

    # print(f"slope: {model.coef_}")
    # print(f"rsquared: {r_sq}")
    # print(count)
    
    rsq_values[count] = r_sq
    lag_values[count] = lag
    
    count += 1

plt.plot(lag_values,rsq_values)
plt.xlabel('lag (ms)')
plt.ylabel('r2')
plt.show