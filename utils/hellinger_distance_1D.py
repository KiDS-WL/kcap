"""
calculates sample-based tension between two 1D distributions
"""

import scipy.stats
from scipy.spatial.distance import euclidean
import numpy as np


### code

# analytic Jeffreys divergence for 2 1D Gaussians with free mean and std
def kl_symm_gauss(mean1,std1,mean2,std2):
    s1sq = np.square(std1)
    s2sq = np.square(std2)
    msq = np.square(mean1) - np.square(mean2)
    return((np.square(s1sq - s2sq) + (s1sq+s2sq)*np.square(msq))/(2.*s1sq*s2sq))


# numerical Jeffreys divergence based on two samples - NOT STABLE
def kl_symm_sample(sample1,sample2):
    thresh = 1.e-111  # min. histogram value to avoid 0
    bin_low = np.floor(np.minimum(min(sample1),min(sample2)))
    bin_upp = np.ceil(np.maximum(max(sample1),max(sample2)))
    bin_n = int(np.ceil(np.sqrt(np.maximum(len(sample1),len(sample2)))))
    #print(bin_low,bin_upp,bin_n)
    bins = np.linspace(bin_low,bin_upp,bin_n)
    histogram1, bins = np.histogram(sample1, bins=bins, density=True)
    histogram2, bins = np.histogram(sample2, bins=bins, density=True)
    histogram1 = np.where(histogram1 < thresh, thresh, histogram1)
    histogram2 = np.where(histogram2 < thresh, thresh, histogram2)
    #print(bins,histogram1,histogram2)   
    kl12 = scipy.stats.entropy(histogram1,histogram2)
    kl21 = scipy.stats.entropy(histogram2,histogram1)
    return(kl12+kl21)


# analytic Helinger distance for 2 1D Gaussians with free mean and std
def hellinger_gauss(mean_diff,std1,std2):
    s1sq = np.square(std1)
    s2sq = np.square(std2)
    expo = -0.25 * np.square(mean_diff) / (s1sq+s2sq)
    root = 2*std1*std2 / (s1sq+s2sq)
    return(np.sqrt(1. - np.sqrt(root) * np.exp(expo)))


# inverse function for given std deviations
def inverse_hellinger_gauss(H,std1,std2):
    s1sq = np.square(std1)
    s2sq = np.square(std2)
    root = np.sqrt(2*std1*std2 / (s1sq+s2sq))
    return( np.sqrt( (-4.)*(s1sq+s2sq) * np.log((1-np.square(H))/root) ) )


# numerical Hellinger distance based on two samples (histogram)
def hellinger_sample(sample1,weight1,sample2,weight2):
    bin_low = np.floor(np.minimum(min(sample1),min(sample2)))
    bin_upp = np.ceil(np.maximum(max(sample1),max(sample2)))
    bin_n = int(np.ceil(np.sqrt(np.maximum(len(sample1),len(sample2)))))
    bins = np.linspace(bin_low,bin_upp,bin_n)
    
    histogram1, bins = np.histogram(sample1, bins=bins, weights=weight1, density=True)
    histogram2, bins = np.histogram(sample2, bins=bins, weights=weight2, density=True)
    histogram1 = histogram1 / np.sum(histogram1)  #normalise
    histogram2 = histogram2 / np.sum(histogram2)  #normalise
    return(euclidean(np.sqrt(histogram1), np.sqrt(histogram2)) / np.sqrt(2.))

# numerical Hellinger distance based on two samples (KDE)
def hellinger_sample_kde(sample1,weight1,sample2,weight2):
    bin_low = np.floor(np.minimum(min(sample1),min(sample2)))
    bin_upp = np.ceil(np.maximum(max(sample1),max(sample2)))
    bin_n = int(np.ceil(np.sqrt(np.maximum(len(sample1),len(sample2)))))
    bins = np.linspace(bin_low,bin_upp,bin_n)

    kde1 = scipy.stats.gaussian_kde(sample1, weights=weight1)
    kde2 = scipy.stats.gaussian_kde(sample2, weights=weight2)
    pdf1 = scipy.stats.gaussian_kde.pdf(kde1,bins)
    pdf2 = scipy.stats.gaussian_kde.pdf(kde2,bins)
    pdf1 = pdf1 / np.sum(pdf1)  # normalise
    pdf2 = pdf2 / np.sum(pdf2)  # normalise
    return(euclidean(np.sqrt(pdf1), np.sqrt(pdf2)) / np.sqrt(2.))


# standard tension estimate for two samples
def classic_tension(sample1,sample2):
    mean1 = np.mean(sample1)
    mean2 = np.mean(sample2)
    var1 = np.var(sample1, ddof=1)
    var2 = np.var(sample2, ddof=1)
    return(np.abs(mean1-mean2)/np.sqrt(var1+var2))


# Hellinger distance-based tension estimate for two samples, each with associated weights
def hellinger_tension(sample1, sample2, weight1=None, weight2=None, htype=3):
    if weight1 is None:
        weight1 = np.ones_like(sample1)
    if weight2 is None:
        weight2 = np.ones_like(sample2)
    
    if htype == 1:
        d_hell = hellinger_sample(sample1, weight1, sample2, weight2)
    elif htype == 2:
        d_hell = hellinger_sample_kde(sample1, weight1, sample2, weight2)
    elif htype == 3:
        d_hell = 0.5*(hellinger_sample(sample1,weight1,sample2,weight2)+hellinger_sample_kde(sample1,weight1,sample2,weight2))
    else:
        raise ValueError("Incorrect value for htype.")
    
    var1 = np.cov(sample1, aweights=weight1, ddof=1)
    var2 = np.cov(sample2, aweights=weight2, ddof=1)

    # find mean_diff for which hellinger_gauss(mean_diff,std1,std2)=d_hell
    mean_diff = inverse_hellinger_gauss(d_hell, np.sqrt(var1), np.sqrt(var2))
    return d_hell, mean_diff/np.sqrt(var1+var2)


if __name__ == "__main__":
    ### settings
    nsample=10000
    m1 = 1.0
    s1 = 0.5
    m2 = 3.5
    s2 = 1.0

    ### main
    np.random.seed(1)

    # create two sampled distributions
    d1 = np.random.normal(loc=m1,scale=s1,size=nsample)
    d2 = np.random.normal(loc=m2,scale=s2,size=nsample)
    w = np.full(nsample,0.5)

    #dist_exp = hellinger_gauss(m1-m2,s1,s2)
    #dist_obs = hellinger_sample(d1,d2)
    #dist_obs2 = hellinger_sample_kde(d1,d2)

    #print(dist_exp,dist_obs,dist_obs2)

    # calculate tension between samples
    tension_c = classic_tension(d1,d2)
    #tension_h1 = hellinger_tension(d1,d2,1)
    #tension_h2 = hellinger_tension(d1,d2,2)
    tension_h = hellinger_tension(d1,w,d2,w)

    print(tension_c,tension_h)
