import sys
import simulation as sim
#from TCDS import simulation as sim
import os
import modelisation_v2 as v2
import numpy as np
import random as rd
import scipy.stats as stats
import pickle
import copy
import math
import os
#import matplotlib.pyplot as plt
import configparser


def fonction_modification_fichiers(INI_file,liste_utile):
    #recuperation en relatif du dossier contenant les informations en entree du genome
    #pour appliquer la casse adaptée au système d'exploitation on utilise join
    dossier_relatif=os.path.join(INI_file.split("/")[0],"tousgenesidentiques") 
    liste_fichier=os.listdir(dossier_relatif)
  
    #recuperation des informations utiles pour la construction du nouveau genome
    position_barriere=liste_utile[1]
    genome_size=liste_utile[2]
    position_debut_fin=liste_utile[0]
    
    
    fichier_barriere,fichier_gff,fichier_debut,fichier_fin="","","",""
    for fichier in liste_fichier:
        if "prot" in fichier:
            fichier_barriere=fichier           
        elif "gff" in fichier:
            fichier_gff=fichier
        elif "TSS" in fichier:
            fichier_debut=fichier
        elif "TTS" in fichier:
            fichier_fin=fichier
    
   
    #creation du fichier de barrieres proteiques    
    with open(os.path.join(dossier_relatif,fichier_barriere), "r") as fichier:
        string_barriere="prot_name\tprot_pos\n"
        for position in position_barriere: 
            position=int(position)
            string_barriere+="hns\t"+str(position)+"\n"
       
        
    with open(os.path.join(dossier_relatif,fichier_barriere), "w") as fichier:
        fichier.write(string_barriere)
        
    #modification du dossier tout gene identique
    string_tout_gene_identique=""        
    liste_position=[]
    for position in position_debut_fin.values():
        liste_position.append(position)
    with open(os.path.join(dossier_relatif,fichier_gff), "r") as fichier:       
        lignes_origines_tout_gene_identique=fichier.read().split("\n")      
        string_tout_gene_identique="\n".join(lignes_origines_tout_gene_identique[0:4])
        
        #pour la 5 eme ligne on conserve tout de la ligne sauf l'element 4  qui correspond a la taille du genome, potentiellement modifiee
        ligne4=lignes_origines_tout_gene_identique[4].split("\t")        
        string_tout_gene_identique+="\n"+"\t".join(ligne4[0:4])+"\t"+str(genome_size)+"\t"+"\t".join(ligne4[5:])+"\n"       
        temp=lignes_origines_tout_gene_identique[5:]
        for index,ligne in enumerate(temp):
            ligne=ligne.split("\t")
            #pour traiter les cas etranges de lignes incompletes
            if (len(ligne)==9):
                #fichier gff format debut toujours avant grand
                min_position=min(liste_position[index][0],liste_position[index][1])
                max_position=max(liste_position[index][0],liste_position[index][1])
                string_tout_gene_identique+="\t".join(ligne[0:3])+"\t"+str(min_position)+"\t"+str(max_position)+"\t.\t"+str(liste_position[index][2])+"\t.\t"+ligne[-1]+"\n"
        
    with open(os.path.join(dossier_relatif,fichier_gff), "w") as fichier:
        fichier.write(string_tout_gene_identique)
   
    
        
   
        
    #modification de tss 
    with open(os.path.join(dossier_relatif,fichier_debut), "r") as fichier:
        lignes_debut_position_genes=fichier.read().split("\n")
        
        string_debut_genes=lignes_debut_position_genes[0]+"\n"
        for index,ligne in enumerate(lignes_debut_position_genes[1:]):
            ligne=ligne.split("\t")
            if (len(ligne)==4):
                    string_debut_genes+=ligne[0]+"\t"+str(liste_position[index][2])+"\t"+str(liste_position[index][0])+"\t"+ligne[-1]+"\n"        
    with open(os.path.join(dossier_relatif,fichier_debut), "w") as fichier:
        fichier.write(string_debut_genes)
        
        
    with open(os.path.join(dossier_relatif,fichier_fin), "r") as fichier:
        lignes_fin_position_genes=fichier.read().split("\n")       
        string_fin_genes=lignes_fin_position_genes[0]+"\n"
        for index,ligne in enumerate(lignes_fin_position_genes[1:]):
            ligne=ligne.split("\t")
            if (len(ligne)==4):
                string_fin_genes+=ligne[0]+"\t"+str(liste_position[index][2])+"\t"+str(liste_position[index][1])+"\t"+ligne[-1]+"\n"
    with open(os.path.join(dossier_relatif,fichier_fin), "w") as fichier:
        fichier.write(string_fin_genes)
       
        
      
  
def calcul_fitness(fitness_cible,fitness_observe):
    #creation matrice cible
    fitness_cible_matrice=np.array(fitness_cible)      
    #creation matrice observe
    fitness_observe_matrice=np.array(fitness_observe)
    #transformation en frequence
    fitness_observe_matrice=fitness_observe_matrice/np.sum(fitness_observe_matrice)
    fitness_global=math.exp(-np.sum(np.abs(np.log(fitness_observe_matrice/fitness_cible_matrice))))   
    return (fitness_global)

def recuperer_fitness_cible(INI_file):
    nom_fichier=os.path.join(INI_file.split("/")[0],"environment.dat")   
    liste_fitness_optimal=[]
    with open(nom_fichier, "r") as fichier:
        lignes_fitness=fichier.read().split("\n") 
        
        for ligne in lignes_fitness:    
            ligne=ligne.split("\t")
            #on ajoute pour chaque ligne le nombre de transcrits
            if len(ligne)==2:
                liste_fitness_optimal.append(float(ligne[1]))
            
    return (liste_fitness_optimal)

def application_monte_carlo(fitness_avant_mutation,fitness_apres_mutation):
    if (fitness_apres_mutation>=fitness_avant_mutation):
        return True
    else:
        proba_garder_mutation= np.exp(-1*(fitness_avant_mutation-fitness_apres_mutation)/0.00001)
        print(proba_garder_mutation)
        garder=rd.random()<proba_garder_mutation
        return garder
            
    
def recuperation_liste_globale(nom_fichier_all_data):   
    with open(nom_fichier_all_data, 'rb') as fichier:
        mon_depickler = pickle.Unpickler(fichier)
        score_recupere = mon_depickler.load()
    return(score_recupere)
   
       

def change_parametre(config_file_path,nom_label,nom_categorie,valeur):
    config = configparser.ConfigParser()
    config.read(config_file_path)
    config.set(nom_label,nom_categorie,str(valeur))
    with open(config_file_path, 'w') as configfile:
        config.write(configfile)
  
   
def simulation(nb_simulations):
    #definition des listes de stockage de variations de fitness
    ancien_genome_optimal=[]
    variation_fitness=[]
    variation_nb_transcrits=[]
    liste_changement_position={}
    for simulation in range (nb_simulations):    
    	[position_gene,pos_barriere,taille_genome,rapport_mutation,nb_transcrits,pas_espace,taille_indel]=sim.start_transcribing(INI_file) 
    	fitness=calcul_fitness(fitness_optimal,nb_transcrits) 
    	print("on en est a la simulation ",simulation)
    	#on calcule dans tous les cas un fitness de base
    	if (simulation==0):   
        	liste_debut=[copy.deepcopy(position_gene) ,pos_barriere.copy(),taille_genome]
        	variation_nb_transcrits.append(nb_transcrits)
        	liste_changement_position["genome_originelle"]=position_gene
        	ancien_genome_optimal=[position_gene,pos_barriere,taille_genome]
    	elif application_monte_carlo(variation_fitness[-1],fitness):
        	#la mutation apporte une amelioration de la fitness ou avec une certaine probabilite est conservee meme si deletere
        	#on met alors a jour nos variables globales de variation de fitness
        	cle="generation"+str(simulation)
        	liste_changement_position[cle]=position_gene        
        	variation_nb_transcrits.append(nb_transcrits) 
        	ancien_genome_optimal=[position_gene,pos_barriere,taille_genome]
    	else:
        	#on est dans le cas ou il n'y a pas eu d'ameliorations de la fitness, dans ce cas, on revient a l'ancien genome le plus optimal
        	fonction_modification_fichiers(INI_file,ancien_genome_optimal)
    	variation_fitness.append(fitness)
    	#On cree notre nouveau genome, qui sera eventuellement conserve
    	#print("voila ce qu'on envoie a la mutation ",position_gene,pos_barriere,taille_genome,rapport_mut,pas_espace)
    	new_genome=v2.create_new_genome(position_gene,pos_barriere,int(taille_genome+1),rapport_mutation,int(pas_espace),int(taille_indel))  
    	#print( "nouveau genome apres mutation est ",new_genome, "\n\n")
    	fonction_modification_fichiers(INI_file,new_genome)
    
    return variation_fitness    
        
#############main function#################
    
           
    
    
    
        
INI_file=sys.argv[1]
#fixation par avance de certains parametres

    #on prend un seuil de 5%   
#recuperation fitness optimal
fitness_optimal=recuperer_fitness_cible(INI_file)


[position_gene,pos_barriere,taille_genome,rapport_mutation,nb_transcrits,pas_espace,taille_indel]=sim.start_transcribing(INI_file) 
liste_debut=[copy.deepcopy(position_gene) ,pos_barriere.copy(),taille_genome]

#Test à part pour le nombre de générations
    
liste_nb_simulations = [200,450]
#liste_nb_simulations = [2,5]

path = "/home/bioscience/users/Paul_Alexis_et_bastien/reseau_bio/Resultats/nb_simulations"
os.mkdir(path)
liste_fitness_end = []
#Simulations pour les valeurs d'un paramètre + création d'un time series pour chaque valeur de paramètre
for nb_simulations in liste_nb_simulations:
	variation_fitness = simulation(nb_simulations)
	liste_fitness_end.append(variation_fitness[-1])
	fonction_modification_fichiers(INI_file,liste_debut)
with open(path + '/Paramètres nombre de simulations.txt', 'w') as f:
    f.write("Paramètre       fitness\n")
    for o in range(len(liste_fitness_end)):
        f.write("%s       %s\n" %(str(liste_nb_simulations[o]), str(liste_fitness_end[o])))

"""
fig, ax = plt.subplots( nrows=1, ncols=1 )
ax.plot(liste_nb_simulations, liste_fitness_end)
ax.set_xlabel('nombre de générations')
ax.set_ylabel('fitness')
fig.suptitle('Paramètre nombre de générations', fontsize=16)	
fig.savefig(path + '/fitness.png')
plt.close(fig)
"""



#simulations pour plusieurs generations

nb_simulations=int(sys.argv[2])

test_parametres = {"['GLOBAL','DELTA_X']":[40,90], "['SIMULATION','RNAPS_NB']":[3,7],
	"['MUTATION','rapport_mutation_insert_invert']":[15,7],"['SIMULATION','GYRASE_CONC']":[0.2,0.5,0.9],
	"['MUTATION', 'taille_indel']":[1,5]}
'''

test_parametres = {"['GLOBAL','DELTA_X']":[20], "['SIMULATION','SIM_TIME']":[1000],
	"['MUTATION','rapport_mutation_insert_invert']":[0.2],"['SIMULATION','GYRASE_CONC']":[0.1],
	"['MUTATION', 'taille_indel']":[1]}

'''

valeurs_initiales = {"['GLOBAL','DELTA_X']":60, "['SIMULATION','RNAPS_NB']":4, "['MUTATION','rapport_mutation_insert_invert']":25,"['SIMULATION','GYRASE_CONC']":0.3,
	"['MUTATION', 'taille_indel']":3}

i = 0
for i in range(len(test_parametres)):
	nom_label = eval(list(test_parametres.keys())[i])[0]
	nom_categorie = eval(list(test_parametres.keys())[i])[1]
	liste_valeurs = list(test_parametres.values())[i]
	path = "/home/bioscience/users/Paul_Alexis_et_bastien/reseau_bio/Resultats/"+eval(list(test_parametres.keys())[i])[1]
	os.mkdir(path)
	liste_fitness_end = []
	print('Test pour',nom_categorie)
	#Simulations pour les valeurs d'un paramètre + création d'un time series pour chaque valeur de paramètre
	for val in liste_valeurs:
		print('valeur :',val)
		change_parametre(INI_file,nom_label,nom_categorie,val)
		variation_fitness = simulation(nb_simulations)
		liste_fitness_end.append(variation_fitness[-1])
		fonction_modification_fichiers(INI_file,liste_debut)
		with open(path + '/Paramètres' + str(val) + '.txt', 'w') as f:
			f.write("fitness\n")
			for j in range(len(variation_fitness)):
				f.write("%s\n" %str(variation_fitness[j]))
	valeur = list(valeurs_initiales.values())[i]
	print(valeur)
	change_parametre(INI_file,nom_label,nom_categorie,valeur)
	with open(path + '/Paramètres ' + str(nom_categorie) +  '.txt', 'w') as f:
		f.write("Paramètre       fitness\n")
		for k in range(len(liste_fitness_end)):
			f.write("%s       %s\n" %(str(liste_valeurs[k]), str(liste_fitness_end[k])))





    
    
#enregistrement sous forme d'objets des donnees
"""
nom_fichier_all_data=os.path.join(INI_file.split("/")[0],"donnee_globale") 
with open(nom_fichier_all_data, 'wb') as fichier:
    mon_pickler = pickle.Pickler(fichier)
    #on concatene les 3 varibales utiles
    liste_globale=[liste_changement_position,variation_fitness,variation_nb_transcrits]    
    mon_pickler.dump(liste_globale)   
"""
        
#remise a jour automatique avec le fichier de depart une fois la simulation terminee



    
    
    
    
