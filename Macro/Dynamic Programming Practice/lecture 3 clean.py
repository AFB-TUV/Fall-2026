import numpy as np
import numba
from scipy import optimize
import matplotlib.pyplot as plt

# some useful plot defaults
plt.rcParams.update({'font.size' : 20, 'lines.linewidth' : 3.5, 'figure.figsize' : (13,7)})

def discretize_assets(amin, amax, n_a):
    # find maximum ubar of uniform grid corresponding to desired maximum amax of asset grid
    ubar = np.log(1 + np.log(1 + amax - amin))
    
    # make uniform grid
    u_grid = np.linspace(0, ubar, n_a)
    
    # double-exponentiate uniform grid and add amin to get grid from amin to amax
    return amin + np.exp(np.exp(u_grid) - 1) - 1

def rouwenhorst_Pi(N, p):
    # base case Pi_2
    Pi = np.array([[p, 1 - p],
                   [1 - p, p]])
    
    # recursion to build up from Pi_2 to Pi_N
    for n in range(3, N + 1):
        Pi_old = Pi
        Pi = np.zeros((n, n))
        
        Pi[:-1, :-1] += p * Pi_old
        Pi[:-1, 1:] += (1 - p) * Pi_old
        Pi[1:, :-1] += (1 - p) * Pi_old
        Pi[1:, 1:] += p * Pi_old
        Pi[1:-1, :] /= 2
        
    return Pi

def stationary_markov(Pi, tol=1E-14):
    # start with uniform distribution over all states
    n = Pi.shape[0]
    pi = np.full(n, 1/n)
    
    # update distribution using Pi until successive iterations differ by less than tol
    for _ in range(10_000):
        pi_new = Pi.T @ pi
        if np.max(np.abs(pi_new - pi)) < tol:
            return pi_new
        pi = pi_new

def discretize_income(rho, sigma, n_e):
    # choose inner-switching probability p to match persistence rho
    p = (1+rho)/2
    
    # start with states from 0 to n_e-1, scale by alpha to match standard deviation sigma
    e = np.arange(n_e)
    alpha = 2*sigma/np.sqrt(n_e-1)
    e = alpha*e
    
    # obtain Markov transition matrix Pi and its stationary distribution
    Pi = rouwenhorst_Pi(n_e, p)
    pi = stationary_markov(Pi)
    
    # e is log income, get income y and scale so that mean is 1
    y = np.exp(e)
    y /= np.vdot(pi, y)
    
    return y, pi, Pi

def backward_iteration(Va, Pi, a_grid, y, r, beta, eis):
    # step 1: discounting and expectations
    Wa = (beta * Pi) @ Va
    
    # step 2: solving for asset policy using the first-order condition
    c_endog = Wa**(-eis)
    coh = y[:, np.newaxis] + (1+r)*a_grid
    a = np.empty_like(coh)
    for e in range(len(y)):
        a[e, :] = np.interp(coh[e, :], c_endog[e, :] + a_grid, a_grid)
        
    # step 3: enforcing the borrowing constraint and backing out consumption
    a = np.maximum(a, a_grid[0])
    c = coh - a
    
    # step 4: using the envelope condition to recover the derivative of the value function
    Va = (1+r) * c**(-1/eis)
    return Va, a, c

def policy_ss(Pi, a_grid, y, r, beta, eis, tol=1E-5):
    # initial guess for Va: assume consumption 5% of cash-on-hand, then get Va from envelope condition
    coh = y[:, np.newaxis] + (1+r)*a_grid
    c = 0.05*coh
    Va = (1+r) * c**(-1/eis)
    
    # iterate until maximum distance between two iterations falls below tol, fail-safe max of 10,000 iterations
    for it in range(10_000):
        Va, a, c = backward_iteration(Va, Pi, a_grid, y, r, beta, eis)
        #if it < 200:
        #    print("iteration ", it, ":   ", a.sum())

        # after iteration 0, can compare new policy function to old one
        if it > 0 and np.max(np.abs(a - a_old)) < tol:
            return Va, a, c
        
        a_old = a


a_grid = discretize_assets(0, 10_000, 500)
y, pi, Pi = discretize_income(0.975, 0.7, 7)

r = 0.01/4
beta = 1-0.08/4
eis = 1

Va, a, c = policy_ss(Pi, a_grid, y, r, beta, eis)
#coh = y[:,np.newaxis] + (1+r)*a_grid
#c = 0.05*coh
#Va = (1+r) * c**(-1/eis)
#Wa = (beta * Pi) @ Va
#c_endog = Wa**(-eis)
#coh = y[:, np.newaxis] + (1+r)*a_grid
#print(c_endog.shape)

for e, ye in enumerate(y):
    plt.plot(a_grid, c[e, :], label=f'y={ye:.2f}')
plt.legend()
plt.xlabel('assets')
plt.ylabel('consumption');

for e, ye in enumerate(y):
    plt.plot(a_grid[:120], c[e, :120], label=f'y={ye:.2f}')
plt.legend()
plt.xlabel('assets')
plt.ylabel('consumption');

for e, ye in enumerate(y[:4]):
    plt.plot(a_grid[:120], c[e, :120], label=f'y={ye:.2f}')
plt.legend()
plt.xlabel('assets')
plt.ylabel('consumption');

for e, ye in enumerate(y):
    plt.plot(a_grid[:300], a[e, :300] - a_grid[:300], label=f'y={ye:.2f}')
plt.axhline(0, linestyle='--', color='k')
plt.legend()
plt.xlabel('assets')
plt.ylabel('savings');

mpcs = np.empty_like(c)
    
# symmetric differences away from boundaries
mpcs[:, 1:-1] = (c[:, 2:] - c[:, 0:-2]) / (a_grid[2:] - a_grid[:-2]) / (1+r)

# asymmetric first differences at boundaries
mpcs[:, 0]  = (c[:, 1] - c[:, 0]) / (a_grid[1] - a_grid[0]) / (1+r)
mpcs[:, -1] = (c[:, -1] - c[:, -2]) / (a_grid[-1] - a_grid[-2]) / (1+r)

# special case of constrained
mpcs[a == a_grid[0]] = 1

for e in range(7):
    plt.plot(a_grid[:50], mpcs[e, :50], linewidth=3, label=f'y = {y[e]:.2f}')
plt.xlabel('Assets')
plt.ylabel('MPC')
plt.title('Quarterly marginal propensities to consume by income state y(s)')
plt.legend();

def get_lottery(a, a_grid):
    # step 1: find the i such that a' lies between gridpoints a_i and a_(i+1)
    a_i = np.searchsorted(a_grid, a) - 1
    
    # step 2: obtain lottery probabilities pi
    a_pi = (a_grid[a_i+1] - a)/(a_grid[a_i+1] - a_grid[a_i])
    
    return a_i, a_pi

a_i, a_pi = get_lottery(a, a_grid)

@numba.njit
def forward_policy(D, a_i, a_pi):
    Dend = np.zeros_like(D)
    for e in range(a_i.shape[0]):
        for a in range(a_i.shape[1]):
            # send pi(e,a) of the mass to gridpoint i(e,a)
            Dend[e, a_i[e,a]] += a_pi[e,a]*D[e,a]
            
            # send 1-pi(e,a) of the mass to gridpoint i(e,a)+1
            Dend[e, a_i[e,a]+1] += (1-a_pi[e,a])*D[e,a]
            
    return Dend

def forward_iteration(D, Pi, a_i, a_pi):
    Dend = forward_policy(D, a_i, a_pi)    
    return Pi.T @ Dend

def distribution_ss(Pi, a, a_grid, tol=1E-10):
    a_i, a_pi = get_lottery(a, a_grid)
    
    # as initial D, use stationary distribution for e, plus uniform over a
    pi = stationary_markov(Pi)
    D = pi[:, np.newaxis] * np.ones_like(a_grid) / len(a_grid)
    
    # now iterate until convergence to acceptable threshold
    for _ in range(10_000):
        D_new = forward_iteration(D, Pi, a_i, a_pi)
        if np.max(np.abs(D_new - D)) < tol:
            return D_new
        D = D_new

D = distribution_ss(Pi, a, a_grid)
#print(D[:,20])

i = np.argmax(a_grid > 20) # first gridpoint above 20
plt.plot(a_grid[:i], D.sum(axis=0)[:i].cumsum())
plt.hlines(1,0, a_grid[i-1], colors='k', linestyles='dotted');

i = np.argmax(a_grid > 2) # first gridpoint above 2
plt.plot(a_grid[:150], D.sum(axis=0)[:150].cumsum())
plt.hlines(1,0, a_grid[149], colors='k', linestyles='dotted');

for e, ye in enumerate(y):
    plt.plot(a_grid[:150], D[e][:150].cumsum()/pi[e], label=f'y={ye:.2f}')
plt.legend();

def steady_state(Pi, a_grid, y, r, beta, eis):
    Va, a, c = policy_ss(Pi, a_grid, y, r, beta, eis)
    D = distribution_ss(Pi, a, a_grid)
    
    return dict(D=D, Va=Va, 
                a=a, c=c,
                A=np.vdot(a, D), C=np.vdot(c, D), # aggregation
                Pi=Pi, a_grid=a_grid, y=y, r=r, beta=beta, eis=eis)

rs = r + np.linspace(-0.02, 0.015, 15)
As = [steady_state(Pi, a_grid, y, r, beta, eis)['A'] for r in rs]
plt.plot(rs, As)
plt.xlabel('Real interest rate')
plt.ylabel('Aggregate assets over quarterly income');

print(As)
'''
sigmas = np.linspace(0.3, 1.2, 8) # our benchmark was sigma=0.7
As = []
for sigma in sigmas:
    y_new, pi_new, Pi_new = discretize_income(0.975, sigma, 7)
    As.append(steady_state(Pi_new, a_grid, y_new, r, beta, eis)['A'])
plt.plot(sigmas, As)
plt.xlabel('Cross-sectional standard deviation of income')
plt.ylabel('Assets over quarterly income');

eis_vec = np.linspace(0.4, 2, 10) # our benchmark was eis=1
As = [steady_state(Pi, a_grid, y, r, beta, eis)['A'] for eis in eis_vec]
plt.plot(eis_vec, As)
plt.xlabel('Elasticity of intertemporal substitution')
plt.ylabel('Assets over quarterly income');

beta_calib = optimize.brentq(
                lambda beta: steady_state(Pi, a_grid, y, r, beta, eis)['A'] - 5.6,
                0.98, 0.995)
beta_calib

ss_calib = steady_state(Pi, a_grid, y, r, beta_calib, eis)
ss_calib['A'], ss_calib['C']

ss_calib['C'] - (1 + ss_calib['r']*ss_calib['A'])

B = 5.6   # annual bonds/GDP is 140%, so quarterly is 560%, and quarterly GDP is 1
tau = r*B # labor tax needed to balance steady-state government budget
e = y     # use our previous y, which had mean 1, for the labor endowment process

beta_ge = optimize.brentq(
                lambda beta: steady_state(Pi, a_grid, (1-tau)*e, r, beta, eis)['A'] - B,
                0.98, 0.995)
beta_ge

ss_ge = steady_state(Pi, a_grid, (1-tau)*e, r, beta_ge, eis)
ss_ge['A'] - B, ss_ge['C'] - 1 # asset mkt clearing, goods mkt clearing

sigmas = np.linspace(0.3, 1.2, 8) # our benchmark was sigma=0.7
rs = []
for sigma in sigmas:
    e_new, pi_new, Pi_new = discretize_income(0.975, sigma, 7)
    rs.append(optimize.brentq(
                lambda r: steady_state(Pi_new, a_grid, (1-r*B)*e_new, r, beta_ge, eis)['A'] - B,
                -0.02, 0.015))
plt.plot(sigmas, 4*np.array(rs))
plt.xlabel('Cross-sectional standard deviation of income')
plt.ylabel('Equilibrium real interest rate (annualized)');
'''