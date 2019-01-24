import random
import numpy as np

############################################################################
#Choisir le type de mutation (inversion ou indel), sortie : type de mutation
############################################################################
def choix_indel_inv(mut=0):
    return np.random.choice(['indel','inversion'],p=[0,1/float(mut)])


#########################################################################################
#Application de la mutation indel, sortie : nouveau dictionnaire avec positions modifiees
#########################################################################################
    
def determination_position_interdites (nb_nucleotides,liste_barrieres,dico_genes,taille_genome,inversion=False):    
    positions_interdites = []
    #On interdit d'abord les positions a +/-delta_x des barrieres si insertion, sinon +/-nb_nucleotides si deletion
    
    #espace a conserver selon le type de mutation
    if (nb_nucleotides<0):
        #on est dans le cas d'une délétion
        espace_securise=-nb_nucleotides
    else:
        #on est dans le cas d'une insertion, on peut très bien insérer près d'une barrière ou d'un gène
        espace_securise=0        
    
    
    for bar in liste_barrieres :
        #on teste aussi si on ne passe pas de l'autre cote du genome, en negatif ou en pos 
        if inversion:
            #si on est en dessous du genome (car circulaire)
            if (bar-espace_securise)<1:
                positions_interdites =positions_interdites +list(range(taille_genome-espace_securise+bar+1,taille_genome+1))+ list(range(1,bar+espace_securise+1))
            #si la position arrive apres
            elif (bar+espace_securise)>taille_genome:
                positions_interdites =positions_interdites + list(range(bar-espace_securise,taille_genome+1))+list(range(1,-taille_genome+espace_securise+bar))
            else:
                positions_interdites =positions_interdites + list(range(bar-espace_securise, bar+espace_securise+1))         
        #si pas d'inversion, on peut se content uniquement d'enlever ce qu'il y aavant
        else:
            if (bar-espace_securise)<1:
                positions_interdites =positions_interdites +list(range(taille_genome-espace_securise+bar+1,taille_genome+1))+ list(range(1,bar))
            else:
                positions_interdites =positions_interdites +  list(range(bar-espace_securise, bar+1))
            
            
            
    #Puis on interdit les positions presentes dans les genes (et a -espace_securise d'un gene pour assurer une deletion minimale)
    #on prend ici en compte min et max, puisque ces derniers peuvent avoir ete echanges par une inversion
    for info_gene in dico_genes.values():
        debut,fin=min(info_gene[0],info_gene[1]),max(info_gene[0],info_gene[1])
        if inversion:
            if (debut-espace_securise)<1:  
                positions_interdites = positions_interdites +list(taille_genome-espace_securise+debut+1,taille_genome+1)+ list(range(1,fin+espace_securise+1))
            #si la position arrive apres
            elif (fin+espace_securise)>taille_genome:
                positions_interdites =positions_interdites + list(range(debut-espace_securise, taille_genome+1))+list(range(1,-taille_genome+espace_securise+fin))
            else:
                positions_interdites =positions_interdites + list(range(debut-espace_securise, fin+espace_securise+1)) 
        else:
            if (debut-espace_securise)<1:
                positions_interdites =positions_interdites +list(range(taille_genome-espace_securise+debut+1,taille_genome+1))+ list(range(1,fin+1))
            else:
                positions_interdites =positions_interdites +  list(range(debut-espace_securise, fin+1))
                
                
    #On en deduit les positions autorisees
    positions_aut = list(range(1,taille_genome+1))
    positions_aut = [pos for pos in positions_aut if pos not in positions_interdites]   
    return(positions_aut)
    
    
    
    
def indel(dico_genes, liste_barrieres, taille_genome,delta_x,taille_indel):
    #Choisir entre insertion et deletion (random 50%) et une position.
    #Attention la position doit etre localisee en dehors des genes.

    #Choix insertion/deletion
    choix_indel = str(np.random.choice(['insertion','deletion']))
    
    #nb = random.randrange(1, 6, 1)*delta_x
    nb = taille_indel * delta_x
    #choix du nombre de nucleotides inserees
    if choix_indel == 'insertion' :
        nb=nb   
    else:
        nb=-nb
      
    positions_aut=[]
    positions_aut = determination_position_interdites (nb,liste_barrieres,dico_genes,taille_genome)
    #Et on en choisit une, mais on verifie au prealable qu'il y ait des positions possibles
    if (len(positions_aut)>=1):
        choix_pos = random.choice(positions_aut)
    #si aucune position possible on renvoit tel quel le genoem, mais en precisant l'abasence de modifications
    else:
        print("pas de modificiation possible")
        return(dico_genes, liste_barrieres, taille_genome)
        
    #Updater les positions des genes qui se trouvent apres la mutation
    #Une fois le nb adapte, on update les positions de genes et de barrieres
    for info_gene in dico_genes.values():
        if info_gene[0] > choix_pos :
            #mise a jour des positions de debut et de fin de genes
            info_gene[0] += nb
            info_gene[1] += nb
            
    #Updater les positions de barrieres qui sont apres la mutation
    
    for index,barriere_position in enumerate(liste_barrieres):
        if barriere_position> choix_pos :
            liste_barrieres[index]+=nb
            
    #mise ajour taille genome
    taille_genome += nb
   
        
    #on cree une chaine qui servira apres de cle pour le type de changement
    type_mutation=choix_indel+" "+str(choix_pos)
    print(type_mutation)

    return dico_genes, sorted(liste_barrieres), taille_genome




#############################################################################################################
#Application de l'inversion, sortie : nouveau dictionnaire avec positions et orientations des genes modifiees
#############################################################################################################
def inclus(debut_gene,pos_debut,pos_fin):
    if (pos_debut<pos_fin) and (debut_gene>pos_debut) and(debut_gene<pos_fin):
        return True
    elif (pos_debut>pos_fin) and ((debut_gene>pos_debut) or(debut_gene<pos_fin)):
        return True
    else:
        return False

      
def calcul_new_position(pos_fin, pos_debut, pos_gene_actuelle, taille_genome):
    if (pos_debut<pos_fin):
        return(pos_fin +pos_debut-pos_gene_actuelle)
    else:
        ecart=pos_fin-pos_gene_actuelle if pos_gene_actuelle<pos_fin else pos_fin+taille_genome-pos_gene_actuelle
        return((pos_debut+ecart)%taille_genome)
        
        
        
def inversion(dico_genes, liste_barrieres,taille_genome,delta_x):   
    
    #on interdit des inversions se produisant au sein des genes uniquement
    #memes zones autorisees que pour les insertions
    positions_aut =determination_position_interdites (-delta_x,liste_barrieres,dico_genes,taille_genome,True) 
    
    #on laisse une marge de delat pour les inversions, pour eviter d'avoir des genes consideres comme confondus   
    pos_debut = random.choice(positions_aut)
    #position de fin peut être apres la position de debut, car il s'agit d'un genome circulaire
    pos_fin = random.choice([pos for pos in positions_aut if pos !=pos_debut])       
    
    pos_debut,pos_fin=min( pos_debut,pos_fin),max(pos_debut,pos_fin)
    
    #Inversion : les genes presents dans la zone d'inversion ont leurs positions modifiees ainsi que leur orientation
    for info_gene in dico_genes.values():         
        #on regarde si le gene est dans la zone d'inversion, 2 cas possibles selon si debut et fin sont inverses ou non
        if inclus(pos_debut,pos_fin,info_gene[0]) : #Si le gene est inclus dans la zone d'inversion  
            
            info_gene[0] = calcul_new_position(pos_fin, pos_debut,  info_gene[0], taille_genome)
            info_gene[1] = calcul_new_position(pos_fin, pos_debut, info_gene[1], taille_genome)                   
            if info_gene[2] == '+' : #Inversion du sens
                info_gene[2] = '-'
            else:
                info_gene[2] = '+'                
                
    #Meme chose avec les barrieres
    for index,barriere_position in enumerate(liste_barrieres):        
        if inclus(pos_debut,pos_fin,barriere_position):
            liste_barrieres[index] = calcul_new_position(pos_fin, pos_debut, liste_barrieres[index], taille_genome)
            
    type_mutation="inversion entre "+str(pos_debut)+" et "+str(pos_fin)
    print(type_mutation)
    
    return (dico_genes,sorted(liste_barrieres),taille_genome)       


def create_new_genome(dico_genes,liste_barrieres,taille_genome, mut,delta_x,taille_indel):    
    choix=choix_indel_inv(mut)
    if choix == 'indel' :        
        genome_nouveau=indel(dico_genes, liste_barrieres,taille_genome,delta_x,taille_indel)
    else :        
        genome_nouveau=inversion(dico_genes,liste_barrieres,taille_genome,delta_x)
    return(genome_nouveau)
        
   
if __name__ == '__main__':
    dico_genes =  {0: [1001, 2000, '+'], 1: [4001, 5000, '+'], 2: [7001,8000, '+'], 3: [10001, 11000, '+'], 4: [13001, 14000, '+'], 5: [16001, 17000, '+'], 6: [19001, 20000, '+'], 7: [22001, 23000, '+'], 8: [25001, 26000, '+'], 9: [28001, 29000, '+']}
    liste_barrieres = [1, 3001, 6001, 9001, 12001, 15001, 18001, 21001, 24001, 27001]
    taille_genome = 30000
    mut = 1
    delta_x=60
    sim = create_new_genome(dico_genes,liste_barrieres,taille_genome, mut,delta_x)
    print(sim)
    
    for it in range(0,100):
        #print('\nSIMULATION ', it, '\n')      
        sim = create_new_genome(dico_genes,liste_barrieres,taille_genome, mut,delta_x)
        print(sim)
        
  
    
        

       



