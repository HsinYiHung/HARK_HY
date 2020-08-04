import numpy as np
# %% Preferences

# Relative risk aversion
CRRA = 10
# Discount factor
DiscFac = 0.96

# Survival probabilities from the author's Fortran code
n = 80
#survprob = np.zeros(n+1)
survprob = 0.99845

# Fix indexing problem (fortran starts at 1, python at 0)
#survprob = np.delete(survprob, [0])
# Now we have 80 probabilities of death,
# for ages 20 to 99.

# Labor income

# They assume its a polinomial of age. Here are the coefficients
a=-2.170042+2.700381
b1=0.16818
b2=-0.0323371/10
b3=0.0019704/100

time_params = {'Age_born': 20, 'Age_retire': 65, 'Age_death': 100}
t_start = time_params['Age_born']
t_ret   = time_params['Age_retire'] # We are currently interpreting this as the last period of work
t_end   = time_params['Age_death']

# They assume retirement income is a fraction of labor income in the
# last working period
repl_fac = 0.68212

# Compute average income at each point in (working) life
f = np.arange(t_start, t_ret+1,1)
f = a + b1*f + b2*(f**2) + b3*(f**3)
det_work_inc = np.exp(f)

# Retirement income
det_ret_inc = repl_fac*det_work_inc[-1]*np.ones(t_end - t_ret)

# Get a full vector of the deterministic part of income
det_income = np.concatenate((det_work_inc, det_ret_inc))

# ln Gamma_t+1 = ln f_t+1 - ln f_t
gr_fac = np.exp(np.diff(np.log(det_income)))

# Now we have growth factors for T_end-1 periods.

# Finally define the normalization factor used by CGM, for plots.
# ### IMPORTANT ###
# We adjust this normalization factor for what we believe is a typo in the
# original article. See the REMARK jupyter notebook for details.
norm_factor = det_income * np.exp(0)

# %% Shocks

# Transitory and permanent shock variance from the paper
std_tran_shock = np.sqrt(0.0738)
std_perm_shock = np.sqrt(0.0106)

# Vectorize. (HARK turns off these shocks after T_retirement)
std_tran_vec = np.array([std_tran_shock]*(t_end-t_start))
std_perm_vec = np.array([std_perm_shock]*(t_end-t_start))

# %% Financial instruments

# Risk-free factor
Rfree = 1.02

# Creation of risky asset return distributions

Mu = 0.06 # Equity premium
Std = 0.157 # standard deviation of rate-of-return shocks

RiskyAvg = Mu + Rfree
RiskyStd = Std

# Make a dictionary to specify the rest of params
dict_portfolio = { 
                   # Usual params
                   'CRRA': CRRA,
                   'Rfree': Rfree,
                   'DiscFac': DiscFac,
                    
                   # Life cycle
                   'T_age' : t_end-t_start+1, # Time of death
                   'T_cycle' : t_end-t_start, # Number of non-terminal periods
                   'T_retire':t_ret-t_start+1,
                   'LivPrb': survprob,
                   'PermGroFac': gr_fac.tolist(),
                   'cycles': 1,
        
                   # Income shocks
                   'PermShkStd': std_perm_vec,
                   'PermShkCount': 3,
                   'TranShkStd': std_tran_vec,
                   'TranShkCount': 3,
                   'UnempPrb': 0,
                   'UnempPrbRet': 0,
                   'IncUnemp': 0,
                   'IncUnempRet': 0,
                   'BoroCnstArt': 0,
                   'tax_rate':0.0,
                   
                    # Portfolio related params
                   'RiskyAvg': RiskyAvg,
                   'RiskyStd': RiskyStd,
                   'RiskyCount': 3,
                   'RiskyShareCount': 30,
                  
                   # Grid 
                   'aXtraMin': 0.001,
                   'aXtraMax': 400,
                   'aXtraCount': 400,
                   'aXtraExtra': [None],
                   'aXtraNestFac': 3,
                   
                   # General
                   'vFuncBool': False,
                   'CubicBool': False,
                   
                   # Simulation params
                   'AgentCount': 10,
                   'pLvlInitMean' : np.log(det_income[0]), # Mean of log initial permanent income (only matters for simulation)
                   'pLvlInitStd' : std_perm_shock,  # Standard deviation of log initial permanent income (only matters for simulation)
                   'T_sim': (t_end - t_start+1)*50,
                   
                   # Unused params required for simulation
                   'PermGroFacAgg': 1,
                   'aNrmInitMean': -50.0, # Agents start with 0 assets (this is log-mean)
                   'aNrmInitStd' : 0.0
}

age_plot_params = [20, 30, 55, 75]
