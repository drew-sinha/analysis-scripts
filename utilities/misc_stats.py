import sklearn
import numpy as np
import scipy.stats as scipy_stats

def breusch_pagan(X,Y):
    X, Y = np.array(X), np.array(Y)

    # Fit data with linear fit; TODO implement a moving average.
    linear_fit = sklearn.linear_model.LinearRegression()
    linear_fit.fit(X, Y)
    resid = Y - linear_fit.predict(X)

    resid_var = np.mean(resid)
    linear_resid_fit.fit(X, (resid_var**-2)*resid)
    resid_R2 = linear_resid_fit.score(X ,(resid_var**-2)*resid)
    return 1-scipy.stats.chi2.cdf(0.5*resid_R2,1)
