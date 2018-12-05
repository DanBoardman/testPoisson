"""

Project description
Motivation: You're looking at foot traffic counts from two available storefronts. One is higher, but the samples were taken over just one week and you're not sure how confident to be in the apparent difference. The storefront that more people passed also costs more to rent and if you're wrong it could make it harder to be profitable and cost your job at the company. If only you knew the odds!

Inspiration: Shameless MATLAB Port of Krishnamoorthy's "Two Poisson Means" fortran code. https://github.com/leiferlab/testPoisson

Original inspiration:  5. Two Poisson Means, 6. Power: Two Poisson Means from http://www.ucs.louisiana.edu/~kxk4695/StatCalc.htm

Theory: Krishnamoorthy, K and Thomson, J. (2004) A more powerful test for comparing two Poisson means. Journal of Statistical Planning and Inference, 119, 249-267
"""



# define Poisson E-test function
def Etest(k1, k2, n1 = 1, n2 = 1, d = 0, iside = 2):
    """
    Computes the p-value for the unconditional 'E-test' of difference between two Poisson means.
    k1: sample 1 count
    k2: sample 2 count
    n1: sample 1 duration
    n2: sample 2 duration
    d: difference of means under null hypothesis (population 1 mean - d = population 2 mean)
    iside: 1 for right tail-test or 2 for two-tail test (default)

    Source: Krishnamoorthy, K and Thomson, J. (2004) A more powerful test for comparing two Poisson means. Journal of Statistical Planning and Inference, 119, 249-267
    Matlab source: https://github.com/leiferlab/testPoisson
    Fortran source (original): http://www.ucs.louisiana.edu/~kxk4695/statcalc/pois2pval.for
    """
    # Initialize stats under null hypothesis
    from scipy.stats import poisson
    import math
    import numpy

    lambdaHatK = (k1+k2) / (n1+n2) - d * n1 / (n1+n2)
    varHat = ( k1 / (n1**2) + k2 / (n2**2) )
    t_k1k2 = (k1 / n1 - k2 / n2 - d) / (varHat**(1/2))

    lambdaHat1 = n1*(lambdaHatK + d)
    lambdaHat2 = n2*lambdaHatK

    i1mode = math.floor(lambdaHat1)
    i2mode = math.floor(lambdaHat2)

    # initialize the probability at the i1mode
    pi1mode = poisson.pmf(i1mode, lambdaHat1, 0)
    pi1 = pi1mode

    # initialize the probability at the i2mode
    pi2mode = poisson.pmf(i2mode, lambdaHat2, 0)

    #define sub-function
    def sumi2(iside, n1, n2, lambdaHat2, t_k1k2, i1, pi1, i2mode, pi2mode, d, pValue):
        """
        Carries out the pValue sum over i2
        """

        pi2 = pi2mode

        # first recursive series
        for i2 in range(i2mode, 1000, 1):
            if pi2 < 1*10**-7:
                break
            lambdaHat_i1 = i1/n1
            lambdaHat_i2 = i2/n2
            diffi = lambdaHat_i1 - lambdaHat_i2 - d
            vari = (lambdaHat_i1/n1 + lambdaHat_i2/n2)
            if iside == 1:
                if i1/n1 - i2/n2 <= d:
                    t_i1i2 = 0.0000
                else:
                    t_i1i2 = diffi/(vari**(1/2))
                if t_i1i2 >= t_k1k2:
                    pValue = pValue + pi1*pi2
            elif iside == 2:
                if abs(i1/n1 - i2/n2) <= d:
                    t_i1i2 = 0
                else:
                    t_i1i2 = diffi/(vari**(1/2))
                if abs(t_i1i2) >= abs(t_k1k2):
                    #print(type(pi2))
                    pValue = numpy.float64(pValue) + float(pi1) * float(pi2)
            pi2 = lambdaHat2 * pi2 / (i2 + 1)

        i2 = i2mode - 1
        pi2 = pi2mode
        pi2 = i2mode * pi2 / lambdaHat2

        # second recursive series
        for i2 in range(i2mode - 1, 0, -1):
            #print('2nd loop, pi2 = ' + str(pi2))
            if pi2 <  1*10**-7:
                break
            lambdaHat_i1 = i1/n1
            lambdaHat_i2 = i2/n2
            diffi = lambdaHat_i1 - lambdaHat_i2 - d
            vari = (lambdaHat_i1/n1 + lambdaHat_i2/n2)
            if iside == 1:
                if i1/n1 - i2/n2 <= d:
                    t_i1i2 = 0
                else:
                    t_i1i2 = diffi/(vari**(1/2))
                if t_i1i2 >= t_k1k2:
                    pValue = pValue + pi1 * pi2
            elif iside == 2:
                if abs(i1/n1 - i2/n2) <= d:
                    t_i1i2 = 0
                else:
                    t_i1i2 = diffi/(vari**(1/2))
                if abs(t_i1i2) >= abs(t_k1k2):
                    pValue = pValue + pi1 * pi2
            pi2 = i2 * pi2 / lambdaHat2

        return pValue


    # initialize pValue
    pValue = 0.0000
    # first outer recursive series
    for i1 in range(i1mode, 1000, 1):
        if pi1 < 1*10**-7:
            break
        pValue = sumi2(iside, n1, n2, lambdaHat2, t_k1k2, i1, pi1, i2mode, pi2mode, d,pValue)
        pi1 = lambdaHat1 * pi1 / (i1 + 1)

    i1 = i1mode - 1
    pi1 = pi1mode
    pi1 = i1mode * pi1 / lambdaHat1

    # second outer recursive series
    for i1 in range(i1mode-1,0,-1):
        if pi1 < 1*10**-7:
            break
        pValue = sumi2(iside, n1, n2, lambdaHat2, t_k1k2, i1, pi1, i2mode, pi2mode, d,pValue)
        pi1 = i1 * pi1 / lambdaHat1

    return pValue
