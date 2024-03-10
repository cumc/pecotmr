#!/usr/bin/env python

"""
Markov Chain Monte Carlo (MCMC) sampler for polygenic prediction with continuous shrinkage (CS) priors.

"""


import scipy as sp
from scipy import linalg 
from numpy import random


import math
from numpy import random
import numpy as np


def fpsi(x, alpha, lam):
    f = -alpha*(math.cosh(x)-1.0)-lam*(math.exp(x)-x-1.0)
    return f


def dfpsi(x, alpha, lam):
    f = -alpha*math.sinh(x)-lam*(math.exp(x)-1.0)
    return f


def g(x, sd, td, f1, f2):
    if (x >= -sd) and (x <= td):
        f = 1.0
    elif x > td:
        f = f1
    elif x < -sd:
        f = f2

    return f


def gigrnd(p, a, b):
    # setup -- sample from the two-parameter version gig(lam,omega)
    p = float(p); a = float(a); b = float(b)
    lam = p
    omega = math.sqrt(a*b)

    if lam < 0:
        lam = -lam
        swap = True
    else:
        swap = False

    alpha = math.sqrt(math.pow(omega,2)+math.pow(lam,2))-lam

    # find t
    x = -fpsi(1.0, alpha, lam)
    if (x >= 0.5) and (x <= 2.0):
        t = 1.0
    elif x > 2.0:
        if (alpha == 0) and (lam == 0):
            t = 1.0
        else:
            t = math.sqrt(2.0/(alpha+lam))
    elif x < 0.5:
        if (alpha == 0) and (lam == 0):
            t = 1.0
        else:
            t = math.log(4.0/(alpha+2.0*lam))

    # find s
    x = -fpsi(-1.0, alpha, lam)
    if (x >= 0.5) and (x <= 2.0):
        s = 1.0
    elif x > 2.0:
        if (alpha == 0) and (lam == 0):
            s = 1.0
        else:
            s = math.sqrt(4.0/(alpha*math.cosh(1)+lam))
    elif x < 0.5:
        if (alpha == 0) and (lam == 0):
            s = 1.0
        elif alpha == 0:
            s = 1.0/lam
        elif lam == 0:
            s = math.log(1.0+1.0/alpha+math.sqrt(1.0/math.pow(alpha,2)+2.0/alpha))
        else:
            s = min(1.0/lam, math.log(1.0+1.0/alpha+math.sqrt(1.0/math.pow(alpha,2)+2.0/alpha)))

    # find auxiliary parameters
    eta = -fpsi(t, alpha, lam)
    zeta = -dfpsi(t, alpha, lam)
    theta = -fpsi(-s, alpha, lam)
    xi = dfpsi(-s, alpha, lam)

    p = 1.0/xi
    r = 1.0/zeta

    td = t-r*eta
    sd = s-p*theta
    q = td+sd

    # random variate generation
    while True:
        U = random.random()
        V = random.random()
        W = random.random()
        if U < q/(p+q+r):
            rnd = -sd+q*V
        elif U < (q+r)/(p+q+r):
            rnd = td-r*math.log(V)
        else:
            rnd = -sd+p*math.log(V)

        f1 = math.exp(-eta-zeta*(rnd-t))
        f2 = math.exp(-theta+xi*(rnd+s))
        if W*g(rnd, sd, td, f1, f2) <= math.exp(fpsi(rnd, alpha, lam)):
            break

    # transform back to the three-parameter version gig(p,a,b)
    rnd = math.exp(rnd)*(lam/omega+math.sqrt(1.0+math.pow(lam,2)/math.pow(omega,2)))
    if swap:
        rnd = 1.0/rnd

    rnd = rnd/math.sqrt(a/b)
    return rnd


def mcmc(a, b, phi, sst_dict, n, ld_blk, blk_size, n_iter, n_burnin, thin, chrom, out_dir, beta_std, write_psi, seed):
    print('... MCMC ...')

    # seed
    if seed != None:
        random.seed(seed)

    # derived stats
    beta_mrg = np.array(sst_dict['BETA'], ndmin=2).T
    maf = np.array(sst_dict['MAF'], ndmin=2).T
    n_pst = (n_iter-n_burnin)/thin
    p = len(sst_dict['MAF'])
    n_blk = len(ld_blk)

    # initialization
    beta = np.zeros((p,1))
    psi = np.ones((p,1))
    sigma = 1.0
    if phi == None:
        phi = 1.0; phi_updt = True
    else:
        phi_updt = False

    beta_est = np.zeros((p,1))
    psi_est = np.zeros((p,1))
    sigma_est = 0.0
    phi_est = 0.0

    # MCMC
    for itr in range(1,n_iter+1):
        if itr % 100 == 0:
            print('--- iter-' + str(itr) + ' ---')

        mm = 0; quad = 0.0
        for kk in range(n_blk):
            if blk_size[kk] == 0:
                continue
            else:
                idx_blk = range(mm,mm+blk_size[kk])
                #print(np.diag(1.0/psi[idx_blk].T))
                dinvt = ld_blk[kk]+np.diag(1.0/psi[idx_blk].T[0])
                dinvt_chol = linalg.cholesky(dinvt)
                beta_tmp = linalg.solve_triangular(dinvt_chol, beta_mrg[idx_blk], trans='T') + np.sqrt(sigma/n)*random.randn(len(idx_blk),1)
                beta[idx_blk] = linalg.solve_triangular(dinvt_chol, beta_tmp, trans='N')
                quad += np.dot(np.dot(beta[idx_blk].T, dinvt), beta[idx_blk])
                mm += blk_size[kk]
                #print("dinvt", dinvt)
                #print("dinvt_chol", dinvt_chol)
                #print("beta_tmp", beta_tmp)
                #print("quad", quad)
                #print("idx_blk", idx_blk)
                #print("psi", psi)
                #print("LD", ld_blk[kk])
        #print("beta", beta)
        err = max(n/2.0*(1.0-2.0*sum(beta*beta_mrg)+quad), n/2.0*sum(beta**2/psi))
        sigma = 1.0/random.gamma((n+p)/2.0, 1.0/err)
        #print(sigma)
        #print(beta)

        delta = random.gamma(a+b, 1.0/(psi+phi))
        #print("delta", delta)

        for jj in range(p):
            psi[jj] = gigrnd(a-0.5, 2.0*delta[jj], n*beta[jj]**2/sigma)
            #print(psi[jj])
        psi[psi>1] = 1.0

        if phi_updt == True:
            w = random.gamma(1.0, 1.0/(phi+1.0))
            phi = random.gamma(p*b+0.5, 1.0/(sum(delta)+w))

        # posterior
        if (itr>n_burnin) and (itr % thin == 0):
            beta_est = beta_est + beta/n_pst
            psi_est = psi_est + psi/n_pst
            sigma_est = sigma_est + sigma/n_pst
            phi_est = phi_est + phi/n_pst

    # convert standardized beta to per-allele beta
    if beta_std == 'False':
        beta_est /= sp.sqrt(2.0*maf*(1.0-maf))

    # write posterior effect sizes
    if phi_updt == True:
        eff_file = out_dir + '_pst_eff_a%d_b%.1f_phiauto_chr%d.txt' % (a, b, chrom)
    else:
        eff_file = out_dir + '_pst_eff_a%d_b%.1f_phi%1.0e_chr%d.txt' % (a, b, phi, chrom)


    # write posterior estimates of psi
    if write_psi == 'True':
        if phi_updt == True:
            psi_file = out_dir + '_pst_psi_a%d_b%.1f_phiauto_chr%d.txt' % (a, b, chrom)
        else:
            psi_file = out_dir + '_pst_psi_a%d_b%.1f_phi%1.0e_chr%d.txt' % (a, b, phi, chrom)


    # print estimated phi
    if phi_updt == True:
        print('... Estimated global shrinkage parameter: %1.2e ...' % phi_est )

    print('... Done ...')
    print(beta_est)



# Load sumstats from the text file
sumstats_data = np.loadtxt("sumstats.txt", delimiter="\t", skiprows=1)
sst_dict = {
    'BETA': sumstats_data[:, 0],
    'MAF': sumstats_data[:, 1]
}

# Load LD from the text matrix
ld_matrix = np.loadtxt("LD.txt", delimiter="\t")
ld_blk = {
    0 : ld_matrix
}

blk_size = [ld_blk[x].shape[0] for x in ld_blk]
print(blk_size)

res = mcmc(1,0.5,None, sst_dict,350, ld_blk, blk_size,
            1000, 500,5, 1, "./",False, False, None)


#python run_prs_cs_test.py 
#... MCMC ...
#/home/gw/GIT/software/pecotmr/inst/prototype/run_prs_cs_test.py:42: DeprecationWarning: Conversion of an array with ndim > 0 to a scalar is deprecated, and will error in future. Ensure you extract a single element from your array before performing this operation. (Deprecated NumPy 1.25.)
#  p = float(p); a = float(a); b = float(b)
#--- iter-100 ---
#--- iter-200 ---
#--- iter-300 ---
#--- iter-400 ---
#--- iter-500 ---
#--- iter-600 ---
#--- iter-700 ---
#--- iter-800 ---
#--- iter-900 ---
#--- iter-1000 ---
#/home/gw/GIT/software/pecotmr/inst/prototype/run_prs_cs_test.py:226: DeprecationWarning: Conversion of an array with ndim > 0 to a scalar is deprecated, and will error in future. Ensure you extract a single element from your array before performing this operation. (Deprecated NumPy 1.25.)
#  print('... Estimated global shrinkage parameter: %1.2e ...' % phi_est )
#... Estimated global shrinkage parameter: 4.17e-01 ...
#... Done ...
#[[ 0.08410789]
# [ 0.74511205]
# [-1.18767133]
# [ 0.77778277]
# [-0.52743742]
# [-0.55978771]
# [-1.10077025]
# [ 1.26027262]
# [ 1.43920873]
# [-1.26155271]
# [-0.51717922]
# [ 4.3908392 ]
# [-0.32931302]
# [-1.39188668]
# [ 0.53230795]
# [-0.13282161]]

