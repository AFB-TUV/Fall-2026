import sim_steady_state_fast as sim
ss = sim.steady_state(**sim.example_calibration())

import numpy as np
import numba
from scipy import optimize
import matplotlib.pyplot as plt

# some useful plot defaults
plt.rcParams.update({'font.size' : 20, 'lines.linewidth' : 3.5, 'figure.figsize' : (13,7)})

T = 11 # we'll calculate household responses up to T
ys = np.outer(1 + 0.01*(np.arange(T) == 5), ss['y'])

Va = ss['Va'] # initialize with steady-state marginal value function
a, c = [np.empty((T,) + ss['a'].shape) for _ in range(2)] # empty arrays to store a and c paths
ss_inputs = {k: ss[k] for k in ('Pi', 'a_grid', 'r', 'beta', 'eis')} # dict of non-changing ss inputs
for t in reversed(range(T)):
    Va, a[t], c[t] = sim.backward_iteration(Va=Va, y=ys[t], **ss_inputs)

c_impulse = 100*(c - ss['c']) / ss['c'] # percent change in consumption policy function
i_ave = np.argmax(ss['a_grid'] > ss['A']) # first index of assets above average
for e in (0, 3, 6):
    plt.plot(c_impulse[:, e, i_ave], label=f"y={ss['y'][e]:.2f}")
plt.legend()
plt.title(f"% consumption policy response to +1% y increase,\n at t=4, given average assets");

D = np.empty_like(a)
D[0] = ss['D']

for t in range(T-1):
    # use special function from sim module to get lottery for time-t policy
    a_i, a_pi = sim.interpolate_lottery_loop(a[t], ss['a_grid'])
    
    D[t+1] = sim.forward_iteration(D[t], ss['Pi'], a_i, a_pi)

A, C = np.empty(T), np.empty(T)
for t in range(T):
    A[t] = np.vdot(D[t], a[t])
    C[t] = np.vdot(D[t], c[t])

C_impulse = 100*(C - ss['C'])/ss['C']
A_impulse = 100*(A - ss['A'])/ss['C']
plt.plot(C_impulse, label='C')
plt.plot(A_impulse, '--', label='A')
plt.legend()
plt.title('% change in response to +1% y increase at t=4');

def policy_impulse(ss, shocks, T):
    # check that all values in "shocks" have first dimension T
    assert all(x.shape[0] == T for x in shocks.values())
    
    # extract inputs to backward_iteration function from ss
    inputs = {k: ss[k] for k in ('Va', 'Pi', 'a_grid', 'y', 'r', 'beta', 'eis')}
    
    # create a T*nS*nA array to store each outputs of backward iteration
    Va, a, c  = (np.empty((T,) + ss['Va'].shape) for _ in range(3))

    for t in reversed(range(T)):
        # add this period's perturbation to parameters that are shocked
        for k in ('y', 'r', 'beta'):
            if k in shocks:
                inputs[k] = ss[k] + shocks[k][t]
                
        Va[t], a[t], c[t] = sim.backward_iteration(**inputs)
        
        # use this Va for next iteration
        inputs['Va'] = Va[t]
        
    return Va, a, c

_, a_alt, c_alt = policy_impulse(ss, {'y': ys - ss['y']}, T)
assert np.array_equal(a, a_alt) and np.array_equal(c, c_alt)

def distribution_impulse(ss, a, T):
    assert a.shape[0] == T
    D = np.empty_like(a)
    D[0] = ss['D']
    
    for t in range(T-1):
        a_i, a_pi = sim.interpolate_lottery_loop(a[t], ss['a_grid'])
        D[t+1] = sim.forward_iteration(D[t], ss['Pi'], a_i, a_pi)
        
    return D



def household_impulse(ss, shocks, T):
    Va, a, c = policy_impulse(ss, shocks, T)
    D = distribution_impulse(ss, a, T)
    return dict(D=D, Va=Va, a=a, c=c, # time-varying stuff
                A=np.sum(a*D, axis=(1,2)), C=np.sum(c*D, axis=(1,2))) # aggregate everything else quickly

impulse = household_impulse(ss, {'y': ys - ss['y']}, T)
assert np.max(np.abs(impulse['A'] - A)) < 2E-15 and np.max(np.abs(impulse['C'] - C)) < 2E-15

for t in (0, 5, 10, 15):
    plt.plot(100*(household_impulse(ss, {'r': 0.01*(np.arange(20) == t)}, 20)['C'] - ss['C'])/ss['C'],
                 label=f'shock at {t}')
plt.legend()
plt.title('% agg C response to 1pp r shocks at different dates');

B = 5.6
tau = ss['r']*B
e = ss['y']

beta_ge = optimize.brentq(lambda beta: sim.steady_state(ss['Pi'], ss['a_grid'],
                                    (1-tau)*e, ss['r'], beta, ss['eis'])['A'] - B, 0.98, 0.995)
ss_ge = sim.steady_state(ss['Pi'], ss['a_grid'], (1-tau)*e, ss['r'], beta_ge, ss['eis'])
ss_ge['A'] - B # check asset market clearing

T = 300
Xs = 1 + 0.01*0.95**np.arange(T)

pi = ss_ge['D'].sum(1) # steady-state distribution of e (unchanging)

def impulse_map(rs, Xs):
    # government budget balance at each date: tau_t = r_t*B (given average e=1) 
    taus = rs * B
    ys = np.outer(Xs - taus, e)
    
    # calculate household impulse
    impulse = household_impulse(ss_ge, {'r': rs - ss_ge['r'], 'y': ys - ss_ge['y']}, len(rs))
    
    return impulse['A'] - B, impulse

asset_mkt_error, _  = impulse_map(np.full(T, ss_ge['r']), Xs)
plt.plot(asset_mkt_error)
plt.title('Excess of asset demand over supply');

rs = np.full(T, ss_ge['r']) # initial guess constant
for it in range(400):
    asset_mkt_error, impulse = impulse_map(rs, Xs)
    if it % 10 == 0 and it <= 100:
        plt.plot(asset_mkt_error, label=f'iteration {it}')
    if np.max(np.abs(asset_mkt_error)) < 5E-4:
        print(f'Converged to 0.0005 accuracy after {it} iterations')
        break
    rs[1:] -= 0.002*asset_mkt_error[:-1] # remember: can't change rs[0], which was fixed here before date 0
plt.legend();

rs = np.full(T, ss_ge['r']) # initial guess constant
for it in range(10):
    asset_mkt_error, impulse = impulse_map(rs, Xs)
    plt.plot(asset_mkt_error, label=f'iteration {it}')
    rs[1:] -= 0.003*asset_mkt_error[:-1]
plt.legend();

%%time
J = np.empty((T, T))
h = 1E-4 # small shock to r is 1E-4
no_shock = impulse_map(np.full(T, ss_ge['r']), 1)[0] # "ghost run" with no shock to subtract off
for tshock in range(T):
    J[:, tshock] = (impulse_map(ss_ge['r'] + h*(np.arange(T) == tshock), 1)[0] - no_shock) / h

plt.plot(J[:60, [0, 5, 10, 15, 20]]);

rs = np.full(T, ss_ge['r']) # initial guess constant
Jbar = J[:-1, 1:] # Jbar for our guesses eliminates last row and first column
errs = []
for it in range(30):
    asset_mkt_error, impulse = impulse_map(rs, Xs)
    plt.plot(asset_mkt_error, label=f'iteration {it}')
    err = np.max(np.abs(asset_mkt_error[:-1])) # solve for asset market clearing at 0, 1, ... , T-2
    errs.append(err)
    if err < 1E-10:
        print(f'Asset market clearing up to 12 digits after {it} iterations')
        break
    rs[1:] -= np.linalg.solve(Jbar, asset_mkt_error[:-1]) # adjust r_1, ..., r_(T-1)
plt.legend();

plt.plot(np.log10(np.array(errs)))

plt.plot(rs[1:101] - ss_ge['r'])
plt.title('Impulse response of $r_t^{ante}$ to AR(1) $X_t$ shock');

impulse_s4 = (impulse_map(ss_ge['r'] + h*(np.arange(T) == 4), 1)[1]['a'] - ss_ge['a'])/h
impulse_s5 = (impulse_map(ss_ge['r'] + h*(np.arange(T) == 5), 1)[1]['a'] - ss_ge['a'])/h
impulse_s6 = (impulse_map(ss_ge['r'] + h*(np.arange(T) == 6), 1)[1]['a'] - ss_ge['a'])/h

plt.plot(range(4, -5, -1), impulse_s4[:9, 3, 50], label='r increase at s=4')
plt.plot(range(5, -5, -1), impulse_s5[:10, 3, 50], '--', label='r increase at s=5')
plt.plot(range(6, -5, -1), impulse_s6[:11, 3, 50], '-.', label='r increase at s=6')
plt.xlabel(r'Number of periods until increase in $r$')
plt.ylabel('Asset response vs. steady state')
plt.legend();

%%time
da = impulse_map(ss_ge['r'] + h*(np.arange(T) == T-1), 1)[1]['a'] - ss_ge['a']

%%time
J_alt = np.empty((T, T))
a_ss = np.broadcast_to(ss_ge['a'], da.shape)
for s in range(T):            # for each dr_s from s=0,...,T-1
    a = a_ss.copy()           # start with steady-state asset policy
    a[:s+1] += da[T-1-s:]     # shock at s impacts 'a' up through s, use response to shock at T-1 starting at T-1-s
    D = distribution_impulse(ss_ge, a, T)
    A = np.sum(a*D, axis=(1,2))
    J_alt[:, s] = (A - B - no_shock) / h

np.max(np.abs(J_alt - J))

plt.plot(range(1, 100), np.diag(J)[1:100], label='$J_{t,t}$')
plt.plot(range(1, 100), np.diag(J)[:99], '--', label='$J_{t-1,t-1}$')
plt.plot(range(1, 100), np.diag(J)[1:100] - np.diag(J)[:99], '-.', label='$J_{t,t} - J_{t-1,t-1}$')
plt.legend();

plt.plot(J[:140, [0, 20, 40, 60, 80, 100]]);

a = a_ss.copy()
a[0] += da[T-1-10] # only shock is to asset policy at date 0 from anticipating s=10 then, nothing else
D = distribution_impulse(ss_ge, a, T)
A = np.sum(a*D, axis=(1,2))
anticipation_effect = (A - B - no_shock) / h

plt.plot(anticipation_effect[:100], label='Impulse from anticipating date-10 shock at date 0')
F = J.copy()
F[1:, 1:] -= J[:-1, :-1]
plt.plot(F[:100, 10], '--', label='$J_{t, 10} - J_{t-1, 9}$')
plt.legend();

%%time
F_alt = np.empty((T,T))
for s in range(T):
    # calculate impulse response to anticipation of date-s shock at date 0
    a = a_ss.copy()
    a[0] += da[T-1-s] 
    D = distribution_impulse(ss_ge, a, T)
    A = np.sum(a*D, axis=(1,2))
    F_alt[:, s] = (A - B - no_shock) / h
np.max(np.abs(F_alt - F))

def J_from_F(F):
    J = F.copy()
    for t in range(1, F.shape[0]):
        J[1:, t] += J[:-1, t-1]
    return J

np.array_equal(J_from_F(F), J)

sim.expectation_functions(ss_ge['a'], ss_ge['Pi'], ss_ge['a_i'], ss_ge['a_pi'], T)

%time curlyE = sim.expectation_functions(ss_ge['a'], ss_ge['Pi'], ss_ge['a_i'], ss_ge['a_pi'], T)

plt.plot(curlyE[:130, 3, 0], label='middle income, low assets today')
plt.plot(curlyE[:130, 3, 150], label='middle income, medium assets today')
plt.plot(curlyE[:130, 3, 250], label='middle income, high assets today')
plt.title('Expected assets by horizon for different states today')
plt.legend();

%%time
F = np.empty((T,T))
D1_noshock = sim.forward_iteration(ss_ge['D'], ss_ge['Pi'], ss_ge['a_i'], ss_ge['a_pi'])
for s in range(T):
    # F_(0,s): change in asset policy times steady-state incoming distribution
    F[0, s] = np.vdot(ss_ge['D'], da[T-1-s]) / h
    
    # change in D_1 from this change
    a_i_shocked, a_pi_shocked = sim.interpolate_lottery_loop(ss_ge['a'] + da[T-1-s], ss_ge['a_grid'])
    dD1 = sim.forward_iteration(ss_ge['D'], ss_ge['Pi'], a_i_shocked, a_pi_shocked) - D1_noshock
    
    # use expectation functions to project effect on aggregate
    F[1:, s] = (curlyE[:T-1].reshape(T-1, -1) @ dD1.ravel()) / h

J_alt = J_from_F(F)

np.max(np.abs(J_alt - J))

:

def step1_backward(ss, shock, T, h):
    # preliminaries: D_1 with no shock, ss inputs to backward_iteration
    D1_noshock = sim.forward_iteration(ss['D'], ss['Pi'], ss['a_i'], ss['a_pi'])
    ss_inputs = {k: ss[k] for k in ('Va', 'Pi', 'a_grid', 'y', 'r', 'beta', 'eis')}
    
    # allocate space for results
    curlyY = {'A': np.empty(T), 'C': np.empty(T)}
    curlyD = np.empty((T,) + ss['D'].shape)
    
    # backward iterate
    for s in range(T):
        if s == 0:
            # at horizon of s=0, 'shock' actually hits, override ss_inputs with shock
            shocked_inputs = {k: ss[k] + h*shock[k] for k in shock}
            Va, a, c = sim.backward_iteration(**{**ss_inputs, **shocked_inputs})
        else:
            # now the only effect is anticipation, so it's just Va being different
            Va, a, c = sim.backward_iteration(**{**ss_inputs, 'Va': Va})
        
        # aggregate effects on A and C
        curlyY['A'][s] = np.vdot(ss['D'], a - ss['a']) / h
        curlyY['C'][s] = np.vdot(ss['D'], c - ss['c']) / h
        
        # what is effect on one-period-ahead distribution?
        a_i_shocked, a_pi_shocked = sim.interpolate_lottery_loop(a, ss['a_grid'])
        curlyD[s] = (sim.forward_iteration(ss_ge['D'], ss_ge['Pi'], a_i_shocked, a_pi_shocked) - D1_noshock) / h
        
    return curlyY, curlyD

def jacobian(ss, shocks, T):
    # step 1 for all shocks i, allocate to curlyY[o][i] and curlyD[i]
    curlyY = {'A': {}, 'C': {}}
    curlyD = {}
    for i, shock in shocks.items():
        curlyYi, curlyD[i] = step1_backward(ss, shock, T, 1E-4)
        curlyY['A'][i], curlyY['C'][i] = curlyYi['A'], curlyYi['C']
    
    # step 2 for all outputs o of interest (here A and C)
    curlyE = {}
    for o in ('A', 'C'):
        curlyE[o] = sim.expectation_functions(ss[o.lower()], ss['Pi'], ss['a_i'], ss['a_pi'], T-1)
                                            
    # steps 3 and 4: build fake news matrices, convert to Jacobians
    Js = {'A': {}, 'C': {}}
    for o in Js:
        for i in shocks:
            F = np.empty((T, T))
            F[0, :] = curlyY[o][i]
            F[1:, :] = curlyE[o].reshape(T-1, -1) @ curlyD[i].reshape(T, -1).T
            Js[o][i] = J_from_F(F)
    
    return Js

J_alt2 = Js['A']['r_direct'] - B*Js['A']['X']
# compare to original brute-force Jacobian, fake-news Jacobian
np.max(np.abs(J_alt2 - J)), np.max(np.abs(J_alt2 - J_alt))

plt.plot(Js['A']['X'][:40, [0, 10, 20]]);

plt.plot(Js['C']['X'][:40, [0, 10, 20]]);

plt.plot(Js['C']['T'][:40, [0, 10, 20]]);

plt.plot(G[:50, [0, 10, 20]])
plt.legend(['s=0', 's=10', 's=20'])
plt.title('First-order response of $r^{ante}_t$ to $X_s$ shocks');

rhos = np.array([0.5, 0.8, 0.9, 0.95, 0.975])
dXs = rhos**np.arange(T-1)[:, np.newaxis] # each column is a dX impulse with different persistence
drs = G @ dXs # simple command obtains impulses to all these simultaneously!
plt.plot(drs[:25])
plt.legend([fr'$\rho={rho}$' for rho in rhos])
plt.title('First-order % response of $r^{ante}_t$ to different 1% AR(1) X shocks');

@numba.njit(parallel=True)
def simul_shock(a, epsilons):
    """Take in any impulse response a to epsilon shock, plus path of epsilons, and simulate"""    
    # if I have T_simul epsilons, can simulate length T_simul - T + 1 dY
    # by indexing as eps_(-T-1), ... , eps_(T_simul-T+1) and implementing formula
    T = len(a)
    T_simul = len(epsilons)
    dY = np.empty(T_simul - T + 1) 
    
    a_flipped = a[::-1].copy() # flip because a_s multiplies eps_(t-s)
    for t in numba.prange(T_simul - T + 1):
        dY[t] = np.vdot(a_flipped, epsilons[t:t + T]) # sum as single dot product

    return dY

%time epsilons = np.random.randn(10_000 + len(G) - 1)
simul_shock(0.01*drs[:, 1], epsilons) # burn-in Numba
%time dr_simulated = simul_shock(0.01*drs[:, 1], epsilons)

plt.plot(dr_simulated[:100]);

:

def autocovariances(a):
    T = len(a)
    return np.array(
        [sum(a[s]*a[s+dt] for s in range(T-dt)) for dt in range(T)])

gammas = autocovariances(0.01*drs[:, 1])
plt.plot(gammas[:50] / gammas[0])
plt.title('Autocorrelation of $r_t$ by lag');

T_sim = len(dr_simulated)
gammas_simulated = np.array([np.vdot(dr_simulated[dt:], dr_simulated[:T_sim-dt])/(T_sim-dt) for dt in range(len(drs))])

plt.plot(gammas[:50] / gammas[0], label='theoretical autocorrelation')
plt.plot(gammas_simulated[:50] / gammas_simulated[0], '--', label='simulated autocorrelation')
plt.legend();

def autocovariances_fft(a):
    T = len(a)
    dft = np.fft.rfftn(a, s=(2 * T - 1,), axes=(0,))
    total = dft.conjugate() @ dft.swapaxes(1, 2)
    return np.fft.irfftn(total, s=(2 * T - 1,), axes=(0,))[:T]

%time gammas_new = autocovariances_fft(0.01*drs[:, 1][:, np.newaxis, np.newaxis]).squeeze()
np.max(np.abs(gammas_new - gammas))

# make T * n_Y * n_X array of impulse responses of outcomes of interest vs. shocks
a = 0.01*np.swapaxes(np.stack([dXs, drs]), 0, 1)
a.shape

%time gammas = autocovariances_fft(a)

plt.plot(np.arange(50), gammas[:50, 0, 1] / np.sqrt(gammas[0, 0, 0]*gammas[0, 1, 1]))
plt.plot(-np.arange(50), gammas[:50, 1, 0] / np.sqrt(gammas[0, 0, 0]*gammas[0, 1, 1])) 
plt.title('Correlation between $dX_t$ and $dr_{t+s}$')
plt.xlabel('s');
