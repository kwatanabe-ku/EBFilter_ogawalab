#! /usr/bin/env python

import beta_binomial
import utils
import math

def get_eb_score(vardepth_target, vardepth_controls, epsilon, whole=False):

    """calculate the EBCall score from the numbers of variant reads and depths of tumor and control samples"""

    if whole:
        depthCounts_control_p = [ x[7] for x in vardepth_controls ]
        depthCounts_control_n = [ x[9] for x in vardepth_controls ]
        varCounts_control_p   = [ x[6] for x in vardepth_controls ]
        varCounts_control_n   = [ x[8] for x in vardepth_controls ]
        
        depthCounts_target_p = vardepth_target[7]
        depthCounts_target_n = vardepth_target[9]
        varCounts_target_p   = vardepth_target[6]
        varCounts_target_n   = vardepth_target[8]
    else:
        depthCounts_control_p = [ x[1] for x in vardepth_controls ]
        depthCounts_control_n = [ x[3] for x in vardepth_controls ]
        varCounts_control_p   = [ x[0] for x in vardepth_controls ]
        varCounts_control_n   = [ x[2] for x in vardepth_controls ]
        
        depthCounts_target_p = vardepth_target[1]
        depthCounts_target_n = vardepth_target[3]
        varCounts_target_p   = vardepth_target[0]
        varCounts_target_n   = vardepth_target[2]

    # estimate the beta-binomial parameters for positive and negative strands
    alpha_p, beta_p = beta_binomial.fit_beta_binomial(depthCounts_control_p, varCounts_control_p, epsilon=epsilon)
    alpha_n, beta_n = beta_binomial.fit_beta_binomial(depthCounts_control_n, varCounts_control_n, epsilon=epsilon)
    
    # evaluate the p-values of target mismatch numbers for positive and negative strands
    pvalue_p = beta_binomial.beta_binom_pvalue([alpha_p, beta_p], depthCounts_target_p, varCounts_target_p)
    pvalue_n = beta_binomial.beta_binom_pvalue([alpha_n, beta_n], depthCounts_target_n, varCounts_target_n)

    # perform Fisher's combination methods for integrating two p-values of positive and negative strands
    EB_pvalue = utils.fisher_combination([pvalue_p, pvalue_n])
    EB_score = 0
    if EB_pvalue < 1e-60:
        EB_score = 60
    elif EB_pvalue > 1.0 - 1e-10:
        EB_score = 0
    else:
        EB_score = - round(math.log10(EB_pvalue), 3)

    return EB_score

