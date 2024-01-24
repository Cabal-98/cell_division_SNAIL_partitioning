timescale = 6 # timescale = numero di timestep in 1 ora (es. timescale=6 significa 1t=10 min)

#COSTANTI SIMULAZIONE TEMPORALE

g_mu200 = 2.1e3
g_mZ = 11
g_Z = 0.1e3

k_mu200 = 0.05
k_mZ = 0.5
k_Z = 0.1

Z_0_mu200 = 220e3
S_0_mu200 = 180e3
Z_0_mZ = 25e3
S_0_mZ = 180e3
    
mu200_0 = 10e3

n_Z_mu200 = 3
n_S_mu200 = 2
n_Z_mZ = 2
n_S_mZ = 2

lam_Z_mu200 = 0.1
lam_S_mu200 = 0.1
lam_Z_mZ = 7.5
lam_S_mZ = 10

l     =    [1.0, 0.6, 0.3, 0.1, 0.05, 0.05, 0.05]
gamma_m  = [0.0, 0.04, 0.2, 1.0, 1.0, 1.0, 1.0]
gamma_mu = [0.0, 0.005, 0.05, 0.5, 0.5, 0.5, 0.5]

n_Z_circuit = 6

dt = 0.1

#COSTANTI DI STUDIO PARTIZIONE
#valori per test
p_vector = [0.40]
var_vector = [0.05]
#valori da usare
#p_vector = [0.25,0.30,0.35,0.40,0.45,0.50]
#var_vector = [0.005,0.01,0.02,0.05,0.1]

#COSTANTI DEL RUMORE
#Rumore duplicazione
eta_mean = 0
eta_var  = 1
eta2   = 0.2


state = 4
generation_mean=215e3

