#simulation testing
"""
Simulation to estimate coverage properties of the E-test
"""
# false positive rate (alpha)

# import tests and supporting functions
import poissonEtest as pet
from scipy.stats import poisson
import pandas as pd
import matplotlib.pyplot as plt


# initialize w/ proposed samples/timings
n1 = 1
n2 = 1

# set trials per experiment
iters = 100

# hack warning -- set range and step for expected mean arrival rates
# challenge: change integer to continuous testing
mu_lower = 1
mu_upper = 5
step = 1

# set significance level
alpha = 0.1

# initialize results lists
means = []
fp_rates= []

for mu in range(mu_lower, mu_upper + 1, step):
    mu1 = mu
    mu2 = mu
    positives = 0
    rate = 0.0
    tests = 0
    for i in range(1, iters):
        k1 = sum(poisson.rvs(mu1, size = n1))
        k2 = sum(poisson.rvs(mu2, size = n2))
        if pet.Etest(k1,k2,n1,n2, iside = 2) < alpha:
            positives = positives + 1
        tests = tests + 1
        rate = float(positives) / float(tests)
    means.append(mu)
    fp_rates.append(rate)

#bind in dataframe, just in case
false_pos = pd.DataFrame(
{'means': means,
 'fp_rates': fp_rates
})

# print simulation result
plt.scatter(x = means, y = fp_rates, data = false_pos, alpha = 0.5)
plt.axhline(y = alpha)
plt.show()
print('Overall false positive rate = ' + str(sum(fp_rates)/len(fp_rates)))
print('Experiments = ' + str(len(range(mu_lower, mu_upper + 1, step))))
print('Total trials = ' + str(iters * len(range(mu_lower, mu_upper + 1, step))))
