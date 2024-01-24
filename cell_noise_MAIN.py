from cell_noise_LIBRERIE import *
from cell_noise_COSTANTI import *
from cell_noise_FUNZIONI import *

def main(ncells,max_population,tmax,nprints):
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

    #EVOLUZIONE E DIVISIONI
    for p in p_vector:
        for var in var_vector:
            cellgen(cells,ncells,generation_mean)
            t=0
            i=ncells+1
            while t < tmax:
                t_check,j=checktimes(cells)
                while t_check < t:
                    #print(t_check, t)
                    cell1 = np.copy(simulazione(cells[f"cell{j}"],t_check))
                    cell2 = np.copy(cell1)
                    #print(cell1)
                    cell1[0], cell2[0] = duplicate(cells[f"cell{j}"][0],p,var)
                    #print(cell1[0])
                    #print(cell2[0])
                    #print("\n")
                    cell1[4] = t_check+division_t()
                    cell1[5] = t_check
                    cell2[4] = t_check+division_t()
                    cell2[5] = t_check
                    cells[f"cell{j}"] = np.copy(cell1)
                    cells[f"cell{i}"] = np.copy(cell2)
                    if len(cells.keys()) >= max_population:
                        i=ran_remove(cells)
                    else:
                        i+=1
                    t_check,j=checktimes(cells)
                if t%tprint_width == 0:
                    for i in range(len(cells.keys())):
                        cells[f"cell{i+1}"]=np.copy(simulazione(cells[f"cell{i+1}"],t))
                    t_phenotypes = count_phenotype(cells)
                    t_fractions = t_phenotypes/(t_phenotypes[0]+t_phenotypes[1]+t_phenotypes[2])
                    epi_frac.append(t_fractions[0])
                    hyb_frac.append(t_fractions[1])
                    mes_frac.append(t_fractions[2])
                t+=tstep
                            
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
                                            
            with open("fractions.txt","w") as f:
               f.write("p = %f - var = %f \n" %(p,var))
               f.write("Epiteliali \t Ibride \t Mesenchimali \n")
               for i in range(len(epi_frac)):
                   f.write("%f \t %f \t %f \n" %(epi_frac[i],hyb_frac[i],mes_frac[i]))
               f.write("cellule totali finali: %i \n" %(len(cells.keys())))
               f.write("conta fenotipi:   Epi = %i, Hyb = %i, Mes = %i \n" %(t_phenotypes[0],t_phenotypes[1],t_phenotypes[2]))
               f.write("percent fenotipi: Epi = %f, Hyb = %f, Mes = %f \n" %((t_phenotypes[0]/(t_phenotypes[0]+t_phenotypes[1]+t_phenotypes[2])),(t_phenotypes[1]/(t_phenotypes[0]+t_phenotypes[1]+t_phenotypes[2])),(t_phenotypes[2]/(t_phenotypes[0]+t_phenotypes[1]+t_phenotypes[2]))))
               f.write("\n")
            epi_frac.clear()
            hyb_frac.clear()
            mes_frac.clear()
            cells.clear()


if __name__ == "__main__":
    ncells = 2
    tmax = 100
    nprints = 10

    if __name__ == "__main__":
    ncells, max_population, tmax, nprints = sys.argv[1:]
    print(f"Parametri inseriti: ncells = {ncells} - tmax = {tmax} - nprints = {nprints}")
    main(ncells,max_population,tmax,nprints)
    
