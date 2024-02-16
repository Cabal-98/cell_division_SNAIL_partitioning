from cell_noise_LIBRERIE import *
from cell_noise_COSTANTI import *

#definisce la modalità di generazione delle cellule iniziali
def start_status(s): 
    generation_mean_0=0
    if s==0: # EPITELIALE
        generation_mean_0 = 150e3
    elif s==1: # IBRIDO
        generation_mean_0 = 215e3
    elif s==2: # MESENCHIMALE
        generation_mean_0 = 275e3
    elif s==3: # GENERALE
        generation_mean_0 = 215e3
    elif s==4: # SELEZIONARE LE FRAZIONI
        generation_mean_0 = 215e3
    else:
        print("Errore nell'inizializzazione dello stato")
    return generation_mean_0,s

def save_object(obj, filename):
    with open(filename, 'wb') as outp:  # Overwrites any existing file.
        pickle.dump(obj, outp, pickle.HIGHEST_PROTOCOL)

def read_object(filename):
    with open(filename, 'rb') as file:
        data = pickle.load(file)
    return data

#FUNZIONI DELL'EVOLUZIONE TEMPORALE
def H(a,lam,a0,n):
    output = 1 / (1 + np.power(a/a0,n))
    output = output + lam*(1 - output)
    return output

def M(mu,mu0,j,n_sum):
    output = np.power((mu/mu200_0),j) / np.power((1+(mu/mu0)),n_sum)
    return output

def Y(gamma,mu,mu0,n_sum):
    output=0
    i = np.array(range(n_sum+1))
    def comb(x):
        return scipy.special.binom(n_sum,x)
    output = np.sum(np.apply_along_axis(comb, 0, i)*gamma*M(mu,mu0,i,n_sum))
    return output

def Yi(gamma,mu,mu0,n_sum):
    output=0
    i = np.array(range(n_sum+1))
    def comb(x):
        return scipy.special.binom(n_sum,x)
    output = np.sum(np.apply_along_axis(comb, 0, i)*gamma*i*M(mu,mu0,i,n_sum))
    return output


#FUNZIONI DELL'EVOLUZIONE DELLA POPOLAZIONE
#Generazione Lognormale per generare un valore di SNAIL
def SNAILgen(generation_mean):
    #print(generation_mean)
    sigma = np.sqrt(np.log(2))
    mu = np.log(generation_mean) - 0.5 * sigma**2
    a = np.random.lognormal(mu,sigma,1000)
    a = a[(a>10e3).nonzero()]
    return a[0]

#Partizione di SNAIL tra le due cellule figlie
#QUESTA VERSIONE DELLA FUNZIONE, UTILIZZATA DALL'ARTICOLO, FUNZIONA IN MODO STUPIDO E DIFETTOSO
#def partition(a,a_noise,p,var):
#    r = np.abs(np.random.normal(p,var,1)) #Fluttuazione dovuta alla partizione 
#    if r[0] > 1:
#        r[0]=1
#    b = a + 0.5*a_noise + r[0]*(2*a+a_noise)
#    c = a + 0.5*a_noise - r[0]*(2*a+a_noise)
#    if b<=0 or c<=0:
#        if b<=0:
#            c=c-b-1
#            b=1
#        else:
#            b=b-c-1
#            c=1
#    #print("cell1 = %f, noise = %f, figlie = %f - %f" %(a,a_noise,b,c))
#    return b,c

def partition(a,a_noise,p,var,index):
    if index == 'SNAIL':
        min_cap = 10e3
    else:
        min_cap = 10
    r = np.abs(np.random.normal(p,var,1)) #Fluttuazione dovuta alla partizione 
    if r[0] > 1:
        r[0]=1
    a_doubled = 2*a + a_noise
    b = a_doubled*r[0]
    c = a_doubled*(1-r[0])
    if b<=min_cap or c<=min_cap:
        if b<=min_cap and c>min_cap:
            alpha = min_cap - b
            c=c+alpha
            b=min_cap
        elif c<=min_cap and b>min_cap:
            beta = min_cap - c
            b=b+beta
            c=min_cap
        else:
            b=min_cap
            c=min_cap
    #print("cell1 = %f, noise = %f, figlie = %f - %f" %(a,a_noise,b,c))
    return b,c
    
def partition_2(a,a_noise,p,var,index,r):
    if index == 'SNAIL':
        min_cap = 10e3
    else:
        min_cap = 10
    a_doubled = 2*a + a_noise
    b = a_doubled*r
    c = a_doubled*(1-r)
    if b<=min_cap or c<=min_cap:
        if b<=min_cap and c>min_cap:
            alpha = min_cap - b
            c=c+alpha
            b=min_cap
        elif c<=min_cap and b>min_cap:
            beta = min_cap - c
            b=b+beta
            c=min_cap
        else:
            b=min_cap
            c=min_cap
    #print("cell1 = %f, noise = %f, figlie = %f - %f" %(a,a_noise,b,c))
    return b,c

#Duplicazione SNAIL pre-divisione + rumore di partizione
def duplicate(a,p,var,index):
    b = 0
    c = 0
    if index == 'SNAIL':
        min_cap = 10e3
    else:
        min_cap = 10
    noise = np.random.normal(eta_mean,eta_var,1000) #Fluttuazione dovuta all'errore di duplicazione
    a_noise = noise*eta2*a
    #try:
    a_noise = a_noise[(2*a+a_noise>min_cap).nonzero()]
    if len(a_noise)>0:
        b,c = partition(a,a_noise[0],p,var,min_cap)
    else:
        print(f"Errore, il parametro {index} è inferiore ai limiti: a = {a}")        
        if a<min_cap:
            a = min_cap
        b=a
        c=a

    #except IndexError:
    #    b = 1
    #    c = 1
    #    print(f"errore, la cellula è vuota: a = {a}")
    #print(b,c)
    return b,c

def duplicate_2(a,p,var,index,r):
    b = 0
    c = 0
    if index == 'SNAIL':
        min_cap = 10e3
    else:
        min_cap = 10
    noise = np.random.normal(eta_mean,eta_var,1000) #Fluttuazione dovuta all'errore di duplicazione
    a_noise = noise*eta2*a
    #try:
    a_noise = a_noise[(2*a+a_noise>min_cap).nonzero()]
    if len(a_noise)>0:
        b,c = partition_2(a,a_noise[0],p,var,min_cap,r)
    else:
        print(f"Errore, il parametro {index} è inferiore ai limiti: a = {a}")        
        if a<min_cap:
            a = min_cap
        b=a
        c=a
    return b,c

def cell_division(cell,p,var,division_mode):
    cell1=np.array([0,0,0,0,0,0])
    cell2=np.array([0,0,0,0,0,0])
    if division_mode == 'indipendente' or division_mode == 'timetest':
        cell1[0],cell2[0] = duplicate(cell[0],p,var,'SNAIL')
        cell1[1],cell2[1] = duplicate(cell[1],p,var,'mu200')
        cell1[2],cell2[2] = duplicate(cell[2],p,var,'mZEB')
        cell1[3],cell2[3] = duplicate(cell[3],p,var,'ZEB')
    elif division_mode == 'unito':
        r = np.abs(np.random.normal(p,var,1)) #Fluttuazione dovuta alla partizione 
        if r[0] > 1:
            r[0]=1
        cell1[0],cell2[0] = duplicate_2(cell[0],p,var,'SNAIL',r[0])
        cell1[1],cell2[1] = duplicate_2(cell[1],p,var,'mu200',r[0])
        cell1[2],cell2[2] = duplicate_2(cell[2],p,var,'mZEB',r[0])
        cell1[3],cell2[3] = duplicate_2(cell[3],p,var,'ZEB',r[0])
    elif division_mode == 'nsym':
        cell1[0],cell2[0] = duplicate(cell[0],0.5,var,'SNAIL')
        cell1[1],cell2[1] = duplicate(cell[1],0.5,var,'mu200')
        cell1[2],cell2[2] = duplicate(cell[2],p,var,'mZEB')
        cell1[3],cell2[3] = duplicate(cell[3],p,var,'ZEB')
    elif division_mode == 'csym':
        cell1[0],cell2[0] = duplicate(cell[0],p,var,'SNAIL')
        cell1[1],cell2[1] = duplicate(cell[1],p,var,'mu200')
        cell1[2],cell2[2] = duplicate(cell[2],0.5,var,'mZEB')
        cell1[3],cell2[3] = duplicate(cell[3],0.5,var,'ZEB')
    cell1[4] = cell[4]
    cell2[4] = cell[4]
    cell1[5] = cell[5]
    cell2[5] = cell[5]
    return cell1,cell2

#Genero il tempo della prossima divisione
#Pensavo di assumere ΔT=1 corrispondente a 10 minuti, così da avere un dt pari a 1 minuto e un rate di una divisione ogni Δt=108
def division_t():
    rate1 = 18*timescale
    #Aggiungere le parti commentate in caso di rate diversi per i diversi step
    #rate2 = 
    #rate3 =
    #rate4 =
    #rate_tot = 1/rate1 + 1/rate2 + 1/rate3 + 1/rate4
    #a=np.random.gamma(shape=4, scale=rate_tot)
    a=np.random.gamma(shape=4, scale=rate1/4)
    return a

#verifico il primo/prossimo tempo di divisione
def checktimes(dizionario):
    key, value = min(dizionario.items(), key=lambda x: x[1][4])
    t=value[4]
    j = int(key[4:])
    return t,j





#Stampo decentemente i dizionari
def printd(dizionario):
    i=1
    while i<=len(dizionario.keys()):
        print(dizionario[f"cell{i}"])
        i+=1
        
def printcolumn(vector,name):
    print(name)
    for i in range(len(vector)):
        print(vector[i])
    print("\n")



#conto il numero di cellule in uno stato rispetto ad un altro
#def count_phenotype(dizionario):
#    pheno=np.zeros(3)
#    #Idea di registrare tutti i fenotipi di ogni cellula
#    #Assegnamo i fenotipi: 0=epiteliale, 1=ibrido, 2=mesenchimale
#    for i in range(len(dizionario.keys())):
#        if dizionario[f"cell{i+1}"][2] >= 600:
#            pheno[2]+=1
#        if dizionario[f"cell{i+1}"][2] <= 125:
#            pheno[0]+=1
#        if (dizionario[f"cell{i+1}"][2] > 125) and (dizionario[f"cell{i+1}"][2] < 600):
#            pheno[1]+=1     
#    return pheno

def count_phenotype(dizionario):
    #Idea di registrare tutti i fenotipi di ogni cellula
    #Assegnamo i fenotipi: 0=epiteliale, 1=ibrido, 2=mesenchimale
    pheno=np.zeros(3)
    values = np.stack(list(dizionario.values()))

    pheno[2] = np.sum(values>= 600, axis=0)[2]
    pheno[0] = np.sum(values<= 125, axis=0)[2]
    pheno[1] = np.sum((values> 125)&(values<600), axis=0)[2]
    return pheno

def ran_remove(cells):
    r_index = rnd.randint(1,len(cells.keys()))
    cells.pop(f"cell{r_index}")
    return r_index

#evoluzione temporale del sistema
def simulazione(parameters,t):
    S=parameters[0]
    mu200=parameters[1]
    mZ=parameters[2]
    Z=parameters[3]
    #print(parameters)
    dmu200=np.empty(2)
    dmZ=np.empty(2)
    dZ=np.empty(2)
    if t<0:
        n_t_max = int(150/dt)
    else:
        n_t_max = int((t-parameters[5])/dt)
    
    for i in range(n_t_max-1):
            #calcolo y(t+1) esplicito
            dmu200[0] = ( g_mu200*H(Z,lam_Z_mu200,Z_0_mu200,n_Z_mu200)*H(S,lam_S_mu200,S_0_mu200,n_S_mu200) - mZ*Yi(gamma_mu,mu200,mu200_0,n_Z_circuit) - k_mu200*mu200 )*dt
            dmu200[1] = mu200 + dmu200[0]
            dmZ[0] = ( g_mZ*H(Z,lam_Z_mZ,Z_0_mZ,n_Z_mZ)*H(S,lam_S_mZ,S_0_mZ,n_S_mZ) - mZ*Y(gamma_m,mu200,mu200_0,n_Z_circuit) - k_mZ*mZ )*dt
            dmZ[1] = mZ + dmZ[0]
            dZ[0] = (g_Z*mZ*Y(l,mu200,mu200_0,n_Z_circuit) - k_Z*Z )*dt
            dZ[1] = Z + dZ[0]
            #calcolo f(y(t+1)) esplicito
            dmu200[0] = ( g_mu200*H(dZ[1],lam_Z_mu200,Z_0_mu200,n_Z_mu200)*H(S,lam_S_mu200,S_0_mu200,n_S_mu200) - dmZ[1]*Yi(gamma_mu,dmu200[1],mu200_0,n_Z_circuit) - k_mu200*dmu200[1] )*dt
            dmZ[0] = ( g_mZ*H(dZ[1],lam_Z_mZ,Z_0_mZ,n_Z_mZ)*H(S,lam_S_mZ,S_0_mZ,n_S_mZ) - dmZ[1]*Y(gamma_m,dmu200[1],mu200_0,n_Z_circuit) - k_mZ*dmZ[1] )*dt
            dZ[0] = (g_Z*dmZ[1]*Y(l,dmu200[1],mu200_0,n_Z_circuit) - k_Z*dZ[1] )*dt
            #calcolo y(t+1)=y(t)+dt*f(y(t+1))
            mu200 = mu200 + dmu200[0]
            mZ = mZ + dmZ[0]
            Z = Z +  dZ[0]
    #print('end_simulation')
    #plt.figure('%i' %ngrafico)
    parameters[1] = mu200 
    parameters[2] = mZ
    parameters[3] = Z
    return parameters

#produco la configurazione iniziale di cellule
def cellgen(cells,ncells,generation_mean):
    
    #Configurazione iniziale degli stati
    epi = np.array([50000.0 ,0.0    ,0.0])
    hyb = np.array([1250.0  ,420.0  ,90000.0])
    mes = np.array([0.0     ,1250.0 ,990000.0])
    
#    epi = np.array([30000.0,30.0,2000.0])
#    mes = np.array([1250.0,1000.0,750000.0])
#    hyb = np.array([9000.0,420.0,90000.0])
    
    #Intervalli di stato
    SNAILepi = [0,209000]
    SNAILhyb = [192000,225000]
    SNAILmes = [184000,400000]
    
    cell_state = np.zeros(6)
    
    #Genero cellule 
    i=1
    #Qualche riga che serve per implementare la scelta delle frazioni
    #SE SI VOGLIONO SELEZIONARE LE FRAZIONI DI PARTENZA MODIFICARE QUI
    fractions=[0.99,0.0,0.01] # 0-epi, 1-hyb, 2-mes
    if state == 4:
        fractions[1] = round(fractions[1] * ncells)
        fractions[2] = round(fractions[2] * ncells)
        fractions[0] = ncells - fractions[1] - fractions[2]
        #print("Frazioni: %i - %i - %i" %(fractions[0],fractions[1],fractions[2]))
    else:
        fractions = [ncells,ncells,ncells] 
    
    start_phenotypes = [0,0,0] # 0-epi, 1-hyb, 2-mes
    print(state,fractions)
    #GENERAZIONE CON STATO CASUALE - =3 se genero normalmente, =4 se scelgo le frazioni
    if state==3 or state==4:
        while i<=ncells:
            cell_state[0]=SNAILgen(generation_mean)
            cell_state[4]=0
            if (cell_state[0]<SNAILmes[0]):
                cell_state[1] = epi[0]
                cell_state[2] = epi[1]
                cell_state[3] = epi[2]
                cell_state=simulazione(cell_state,100)
                cell_state[4]=division_t()
                cells[f"cell{i}"] = np.copy(cell_state)
                #print(start_phenotypes)
                if start_phenotypes[0]<fractions[0]:
                    start_phenotypes[0]+=1
                    i+=1
            if (cell_state[0]>SNAILmes[0]) and (cell_state[0]<SNAILhyb[0]):
                cell_state[1] = epi[0]
                cell_state[2] = epi[1]
                cell_state[3] = epi[2]
                cell_state=simulazione(cell_state,100)
                cell_state[4]=division_t()
                cells[f"cell{i}"] = np.copy(cell_state)
                #print(start_phenotypes)
                if start_phenotypes[0]<fractions[0]:
                    start_phenotypes[0]+=1
                    i+=1
                cell_state[1] = mes[0]
                cell_state[2] = mes[1]
                cell_state[3] = mes[2]
                cell_state[4]=0
                cell_state=simulazione(cell_state,100)
                cell_state[4]=division_t()
                cells[f"cell{i}"] = np.copy(cell_state)
                #print(start_phenotypes)
                if start_phenotypes[2]<fractions[2]:
                    start_phenotypes[2]+=1
                    i+=1
            if (cell_state[0]>SNAILhyb[0]) and (cell_state[0]<SNAILepi[1]):
                cell_state[1] = epi[0]
                cell_state[2] = epi[1]
                cell_state[3] = epi[2]
                cell_state=simulazione(cell_state,100)
                cell_state[4]=division_t()
                cells[f"cell{i}"] = np.copy(cell_state) 
                #print(start_phenotypes)
                if start_phenotypes[0]<fractions[0]:
                    start_phenotypes[0]+=1
                    i+=1
                cell_state[1] = hyb[0]
                cell_state[2] = hyb[1]
                cell_state[3] = hyb[2]
                cell_state[4] = 0
                cell_state=simulazione(cell_state,100)
                cell_state[4]=division_t()
                cells[f"cell{i}"] = np.copy(cell_state) 
                #print(start_phenotypes)
                if start_phenotypes[1]<fractions[1]:
                    start_phenotypes[1]+=1
                    i+=1
                cell_state[1] = mes[0]
                cell_state[2] = mes[1]
                cell_state[3] = mes[2]
                cell_state[4] = 0
                cell_state=simulazione(cell_state,100)
                cell_state[4]=division_t()
                cells[f"cell{i}"] = np.copy(cell_state)
                #print(start_phenotypes)
                if start_phenotypes[2]<fractions[2]:
                    start_phenotypes[2]+=1
                    i+=1
            if (cell_state[0]>SNAILepi[1]) and (cell_state[0]<SNAILhyb[1]):
                cell_state[1] = hyb[0]
                cell_state[2] = hyb[1]
                cell_state[3] = hyb[2]            
                cell_state=simulazione(cell_state,100)
                cell_state[4]=division_t()
                cells[f"cell{i}"] = np.copy(cell_state) 
                #print(start_phenotypes)
                if start_phenotypes[1]<fractions[1]:
                    start_phenotypes[1]+=1
                    i+=1
                cell_state[1] = mes[0]
                cell_state[2] = mes[1]
                cell_state[3] = mes[2]
                cell_state[4] = 0
                cell_state=simulazione(cell_state,100)
                cell_state[4]=division_t()
                cells[f"cell{i}"] = np.copy(cell_state) 
                #print(start_phenotypes)
                if start_phenotypes[2]<fractions[2]:
                    start_phenotypes[2]+=1
                    i+=1
            if (cell_state[0]>SNAILhyb[1]) and (cell_state[0]<SNAILmes[1]):
                cell_state[1] = mes[0]
                cell_state[2] = mes[1]
                cell_state[3] = mes[2]
                cell_state=simulazione(cell_state,100)
                cell_state[4]=division_t()
                cells[f"cell{i}"] = np.copy(cell_state) 
                #print(start_phenotypes)
                if start_phenotypes[2]<fractions[2]:
                    start_phenotypes[2]+=1
                    i+=1
    elif state == 0: #GENERAZIONE SOLE EPITELIALI
        while i<=ncells:
            cell_state[0]=SNAILgen(generation_mean)
            cell_state[4]=0
            if (cell_state[0]<SNAILepi[1]):
                cell_state[1] = epi[0]
                cell_state[2] = epi[1]
                cell_state[3] = epi[2]
                cell_state[4] = 0
                cell_state=simulazione(cell_state,100)
                cell_state[4]=division_t()
                cells[f"cell{i}"] = np.copy(cell_state) 
                i+=1
    elif state == 1: #GENERAZIONE SOLE CELLULE IBRIDE
        while i<=ncells:
            cell_state[0]=SNAILgen(generation_mean)
            cell_state[4]=0
            if (cell_state[0]>SNAILhyb[0]) and (cell_state[0]<SNAILhyb[1]):
                cell_state[1] = hyb[0]
                cell_state[2] = hyb[1]
                cell_state[3] = hyb[2]
                cell_state[4] = 0
                cell_state=simulazione(cell_state,100)
                cell_state[4]=division_t()
                cells[f"cell{i}"] = np.copy(cell_state) 
                i+=1
    elif state == 2: #GENERAZIONE SOLE MESENCHIMALI
        while i<=ncells:
            cell_state[0]=SNAILgen(generation_mean)
            cell_state[4]=0
            if (cell_state[0]>SNAILmes[0]):
                cell_state[1] = mes[0]
                cell_state[2] = mes[1]
                cell_state[3] = mes[2]
                cell_state[4] = 0
                cell_state=simulazione(cell_state,100)
                cell_state[4]=division_t()
                cells[f"cell{i}"] = np.copy(cell_state) 
                i+=1
    else:
        print("Errore nell'inizializzazione dei parametri di generazione.")
