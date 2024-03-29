from cell_noise_LIBRERIE import *
from cell_noise_COSTANTI import *
from cell_noise_FUNZIONI import *

def main(ncells,max_population,tmax,nprints,division_mode,output_path,p_input,var_input):
    #CELLULA
    cells = defaultdict(list)
    #Stato: SNAIL=0, mu200=1, mZ=2, Z=3, tnext=4, t0=5 dove t0 è il tempo di ultima divisione e tnext il tempo della prossima

    cell1=np.array([0,0,0,0,0,0])
    cell2=np.array([0,0,0,0,0,0])

    t_phenotypes = np.empty(3)
    t_fractions  = np.empty(3)
    epi_frac = []
    hyb_frac = []
    mes_frac = []

    t = 0
    #Parametri per stampa periodica
    t_check = 0
    tprint_width = int(tmax/nprints)
    tstep = tprint_width
    tprint = []
    while t_check <= tmax:
        tprint.append(t_check)
        t_check+=tprint_width
        barslenght = 0.2
        offset = np.arange(len(tprint))
        #printd(cells)
        #print("\n")
    #TEST GRAFICO DEL SISTEMA GENERATO
    #for w in range(len(cells.keys())):
        #    a = np.copy(cells[f"cell{w+1}"][0])
        #    b = np.copy(cells[f"cell{w+1}"][2])
        #    plt.plot(a,b,'.')

    #Stampa parametri per comodità
    #print("p = %f - var = %f" %(p, var))

    #Crea un unico vettore da cui estrarre p e var cosi da rimuovere un ciclo
    pv_comb = np.array(np.meshgrid(p_vector,var_vector)).T.reshape(-1,2)
    if division_mode == 'timetest':
        pv_comb = np.array([[0.40,0.05]])
    if p_input or var_input != 0:
        pv_comb = np.array([[p_input,var_input]])
    #EVOLUZIONE E DIVISIONI
    for p, var in tqdm(pv_comb):
        print(tprint_width)
        print("\n")
        print(p, var)
        cellgen(cells,ncells,generation_mean)
        t=0
        i=ncells #aggiungere +1 se si sposta l'if e quindi l'aggiornamento dell'indice avviene dopo

        if p_input or var_input != 0:
            pbar = tqdm(total = tmax)
        
        while t < tmax:
            if p_input or var_input != 0:
                pbar.update(tstep)
            t_check,j=checktimes(cells)
            while t_check < t:
                #print(t_check, t)
                cells[f"cell{j}"] = np.copy(simulazione(cells[f"cell{j}"],t_check))
                #print(cell1)
                #cell1[0], cell2[0] = duplicate(cells[f"cell{j}"][0],p,var)
                cell1, cell2 = cell_division(cells[f"cell{j}"],p,var,division_mode)
                #print(cell1[0])
                #print(cell2[0])
                #print("\n")
                if len(cells.keys()) >= max_population:
                    i=ran_remove(cells)
                else:
                    i+=1 #Controllare bene l'indicizzazione perché è stato fornito a "simulazione" una cell vuota
                cell1[4] = t_check+division_t()
                cell1[5] = t_check
                cell2[4] = t_check+division_t()
                cell2[5] = t_check
                cells[f"cell{j}"] = np.copy(cell1)
                cells[f"cell{i}"] = np.copy(cell2)
                #Ho provato a spostare in alto l'if di controllo del numero totale
                t_check,j=checktimes(cells)
            if t%tprint_width == 0:
                for i in range(len(cells.keys())):
                    cells[f"cell{i+1}"]=np.copy(simulazione(cells[f"cell{i+1}"],t))
                t_phenotypes = count_phenotype(cells)
                print(f"\n Timestep raggiunto: t = {t} - Cellule: {len(cells.keys())}")
                t_fractions = t_phenotypes/(t_phenotypes[0]+t_phenotypes[1]+t_phenotypes[2])
                epi_frac.append(t_fractions[0])
                hyb_frac.append(t_fractions[1])
                mes_frac.append(t_fractions[2])
            t+=tstep

        if p_input or var_input != 0:
            pbar.close()
        
        for w in range(len(cells.keys())):
            cells[f"cell{w+1}"] = np.copy(simulazione(cells[f"cell{w+1}"],-1))
            #a = cells[f"cell{w+1}"][0]
            #b = cells[f"cell{w+1}"][2]
            #if a > 450000:
            #    a = 450000
            #plt.figure("SNAILplot")
            #plt.plot(a,b,'.')
            #if a < 450000:
            #    plt.figure("SNAILplot")
            #    plt.plot(a,b,'.')
        t_phenotypes = count_phenotype(cells)
        epi_frac.append(t_phenotypes[0]/(t_phenotypes[0]+t_phenotypes[1]+t_phenotypes[2]))
        hyb_frac.append(t_phenotypes[1]/(t_phenotypes[0]+t_phenotypes[1]+t_phenotypes[2]))
        mes_frac.append(t_phenotypes[2]/(t_phenotypes[0]+t_phenotypes[1]+t_phenotypes[2]))
        
        label = "_" + str(p).replace(".","") + "_" + str(var).replace(".","") + "_" + str(max_population)
        
        with open(output_path + "fractions" + label + ".txt","w") as f:
           f.write("p = %f - var = %f \n" %(p,var))
           f.write("Epiteliali \t Ibride \t Mesenchimali \n")
           for i in range(len(epi_frac)):
               f.write("%f \t %f \t %f \n" %(epi_frac[i],hyb_frac[i],mes_frac[i]))
           f.write("cellule totali finali: %i \n" %(len(cells.keys())))
           f.write("conta fenotipi:   Epi = %i, Hyb = %i, Mes = %i \n" %(t_phenotypes[0],t_phenotypes[1],t_phenotypes[2]))
           f.write("percent fenotipi: Epi = %f, Hyb = %f, Mes = %f \n" %((t_phenotypes[0]/(t_phenotypes[0]+t_phenotypes[1]+t_phenotypes[2])),(t_phenotypes[1]/(t_phenotypes[0]+t_phenotypes[1]+t_phenotypes[2])),(t_phenotypes[2]/(t_phenotypes[0]+t_phenotypes[1]+t_phenotypes[2]))))
           f.write("\n")

        save_object(offset, output_path + "barsxaxis" + label + ".pkl" )
        save_object(epi_frac, output_path + "epi_frac" + label + ".pkl")
        save_object(hyb_frac, output_path + "hyb_frac" + label + ".pkl")
        save_object(mes_frac, output_path + "mes_frac" + label + ".pkl")

        epi_frac.clear()
        hyb_frac.clear()
        mes_frac.clear()
        cells.clear()
        
        #plt.figure("histogram")
        #plt.bar(offset - barslenght, epi_frac, barslenght, label='Epiteliale', color='blue',edgecolor='black')
        #plt.bar(offset,                  hyb_frac, barslenght, label='Ibrido', color='gold',edgecolor='black')
        #plt.bar(offset + barslenght, mes_frac, barslenght, label='Mesenchimale',color='red',edgecolor='black')  
        #plt.xlabel('t')
        #plt.ylabel('Frazioni')
        #plt.title('Frazioni 0.99-0.00-0.01, p=%f, var=%f' %(p,var))
        #plt.xticks(offset, tprint)
        #plt.legend()
        #plt.savefig("./output/plot.png",dpi=1200)

if __name__ == "__main__":
    ncells, max_population, tmax, nprints = list(map(int,sys.argv[1:5]))
    p_input, var_input = list(map(float,sys.argv[5:7]))
    run_number = sys.argv[7]
    division_mode = sys.argv[8]
    output_folder = sys.argv[9]
    print(f"Parametri inseriti: ncells = {ncells} - nmax = {max_population} - max_population = {max_population} - tmax = {tmax} - nprints = {nprints} - p = {p_input} - var = {var_input} - run number = {run_number} - division_mode = {division_mode}")
    #print(f"Parametri inseriti: ncells = {ncells} - max_pop = {max_population} - tmax = {tmax} - nprints = {nprints} - run_number = {run_number} - division_mode = {division_mode}")
    modalita = ['indipendente','unito','nsym','csym','timetest']
    if division_mode in modalita:
        if output_folder == '0':
            output_path = f'./output/run{run_number}/sim_p_{str(p_input).replace(".","")}_sigma_{str(var_input).replace(".","")}/' 
        else:
            output_path = output_folder + f'/run{run_number}/sim_p_{str(p_input).replace(".","")}_sigma_{str(var_input).replace(".","")}/'
        print(output_path)
        if not os.path.exists(output_path):
            os.makedirs(output_path)
        main(ncells,max_population,tmax,nprints,division_mode,output_path,p_input,var_input)
    else:
        print(f"L'input 'division_mode' accetta solo i valori: {modalita}")
