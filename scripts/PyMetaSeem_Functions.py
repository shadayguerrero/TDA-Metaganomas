import os
import math
import gzip
import random
import numpy as np

from pathlib import Path
from random import seed
from random import randint


#######################################
# Funciones de PyMetaSeem #
#######################################

# funcion para lectura de archivos
# ingresa un archivo multifasta y arroja una lista con cada secuencia (contig)
dicc_contigs = {}
def lectura_genoma(File):
    
    lista = []
    longitudes = []
    with open(File,'r') as f:
        lines=f.read() # lectura de cada linea del archivo
        lines=lines.split('>') # identificador de '>'
        lines=['>'+ x for x in lines[1:]] # lista con cada elemento que comienza con '>'
        
        for x in lines:
            x1 = x.replace(">","@")
            x2 = x1.replace("\n",",",1) # el primer '\n' se reemplazapor una coma
            x3 = x2.replace("\n","") # los siguientes se quitan
            lista.append(x3) # lista con un contig en cada elemento y su identificador
    
        # convertir la lista en un diccionario        
        for x in lista:
            x = x.split(',')
            dicc_contigs[x[0]] = x[1] 
        return(dicc_contigs)
    
    
# funcion que cuenta los contigs de cada genoma
def cuenta_contigs(dicc_contigs):
    long_total = 0
    dicc_longitudes = {}
    num_contigs = len(dicc_contigs) # cuenta cuantos contigs contiene cada genoma (cuantas llaves)
    for key in dicc_contigs:
        longitud = len(dicc_contigs[key])# calcula la longitud de cada contig (la longitud de las claves, para cada llave)
        long_total +=  longitud # va sumando todas las longitudes de los contigs
        dicc_longitudes[key] = longitud # diccionario de longitudes 
        # key es el identificador del contig 
        # longitud -> es la longitud de cada contig
        # long_total -> lontitud total del genoma
    return (dicc_longitudes,num_contigs,long_total)   
        
      
        
        
### SE TOMA EL MINIMO DE CANTIDAD DE CONTIGS (se debe comparar por cada genoma)
### TOMAR EL NUMERO MINIMO DE CONTIGS CON MAYOR LONGITUD

# funcion para calcular la proporcion de los contigs
def proporcion_contigs(dicc_contigs,dicc_longitudes,long_total):
    dicc_proporciones = {}
    for key in dicc_longitudes:
        P = dicc_longitudes[key]/long_total
        dicc_proporciones[key] = P
    return dicc_proporciones
# numero de reads por contig se calcula en la funcion reads 
# funcion para cortar reads, dado una pocision de inicio y una longitud de corte
def cutout(contig,i,n_length): #fordward    (resive contig y entrega reads)
    read = contig[i:i+n_length]
    return(read) 

def cutout2(contig,i,n_length,inserto): #backward 
    # inserto = 400 (tamaño promedio del inserto), depende del secuenciador y lo da el usuario
    # inserto -> tamaño verdadero de la molecula de DNA que esta en el secuanciador
    j = i + inserto
    # j = i+inserto (pocision de inicio de reversa)
    # n_length = 150 (tamaño del read "cortado")
    read_cortado = contig[j-n_length:j]
    #i_read = contig[j-n_length:inserto]  #(quiero comenzar en la pocision j e ir hacia atras)
    return read_cortado


# funcion que calcula el reverso de un Read 
def reverso(read_cortado):
    temp_list = list(read_cortado)
    temp_list.reverse()
    return ''.join(temp_list)

def complementario(read_cortadoR):
    comp = []
    for i in range(len(read_cortadoR)):
        if read_cortadoR[i] == "A":
            comp.append("T")
        elif read_cortadoR[i] == "T":
            comp.append("A")
        elif read_cortadoR[i] == "G":
            comp.append("C")
        elif read_cortadoR[i] == "C":
            comp.append("G") 
            
    return "".join(comp)


def reverso_complementario(contig,i,n_length,inserto):
    read_rev_com = complementario(reverso(cutout2(contig,i,n_length,inserto))) 
    return (read_rev_com) # entrega cada read reverso complementario


# funcion para cortar los reads aleatoriamente
# Al ingresar un conjunto de genomas, este nos arroja un metagenoma con varios reads tomados de los genomas iniciales
# para esto se le entrega el numero de reads deseados por cada genoma y la longitud.

def reads(dicc_contigs,dicc_longitudes,dicc_proporciones,n_length,num_reads_total):  
    
    dicc_reads = {}  # diccionario de reads 
    dicc_reads_reverse = {} # diccionario de reverso reads
    k = 0
    kk = 0
    inserto = 400 #400 #(tamaño promedio del inserto), depende del secuenciador y lo da el usuario
     
    # LONGITUD Y NUM_READS DEPENDERA DE LAS LONGITUDES DE LOS CONTIG
    # LA LONGITUD 
    # LAS LONGITUDES DEPENDEN DE LA CALIDAD
    
    while k < num_reads_total: # número de reads total
              
        for key in dicc_contigs:
            #contig = random.choice(list(dicc_contigs.values())) # escoje un contig al azar
            num_reads_contig = round(num_reads_total* dicc_proporciones[key]) # numero de reads x proporcion #proporcion de reads por contig
            contig = dicc_contigs[key]
            # la proporcion viene del diccionario proporciones de la funcion cuenta contigs 
            for kk in range(num_reads_contig): 
                newkey = key + '_' + str(kk) # identificador de contig, con el numero de read
                i = randint(1,(len(contig)-n_length-inserto-1)) # toma una pocision de inicio al azar se debe tener en cuenta el inserto)
                dicc_reads[newkey] = cutout(contig,i,n_length) # para el contig dado anteriormente, se corta desde la pocision i de lonjitud ln
                dicc_reads_reverse[newkey] = reverso_complementario(contig,i,n_length,inserto)# reverso complementario

                kk += 1
            k += num_reads_contig # funcion para cortar los reads aleatoriamente
        # se debe calcular cuantas veces va a escoger cada contig,dependiendo del tamaño del contig inicial      

    return(dicc_reads,dicc_reads_reverse)

# funcion de creacion archivo .fastq por cada genoma
def crea_fastq(dicc_reads,dicc_reads_reverse,n_length,file_name):
    file = open(file_name+'_R1.fastq','wt') 
    for key in dicc_reads: 
        file.write(str(key))
        file.write(str('\n'))
        file.write(str(dicc_reads[key]))
        file.write(str('\n'))
        file.write(str('+'))
        file.write(str('\n'))
        file.write(str('A'*n_length)) #calidad
        file.write(str('\n'))     
    # agrega el reverso complementario  
    file = open(file_name+'_R2.fastq','wt')
    for key in dicc_reads_reverse: 
        file.write(str(key))
        file.write(str('\n'))
        file.write(str(dicc_reads_reverse[key]))
        file.write(str('\n'))
        file.write(str('+'))
        file.write(str('\n'))
        file.write(str('A'*n_length)) #calidad
        file.write(str('\n'))  
    file.close()    
    
# funcion de creacion archivo .fastq.gz por cada genoma
def crea_fastqz(dicc_reads,dicc_reads_reverse,n_length,file_name):
    file = gzip.open(file_name+'_R1.fastq.gz','wt')
    for key in dicc_reads: 
        file.write(str(key))
        file.write(str('\n'))
        file.write(str(dicc_reads[key]))
        file.write(str('\n'))
        file.write(str('+'))
        file.write(str('\n'))
        file.write(str('A'*n_length)) #calidad
        file.write(str('\n'))     
    # agrega el reverso complementario  
    file = open(file_name+'_R2.fastq.gz','wt')
    for key in dicc_reads_reverse: 
        file.write(str(key))
        file.write(str('\n'))
        file.write(str(dicc_reads_reverse[key]))
        file.write(str('\n'))
        file.write(str('+'))
        file.write(str('\n'))
        file.write(str('A'*n_length)) #calidad
        file.write(str('\n'))  
    file.close()
    
def todo(genoma, longitud_read,num_reads,name):
    dicc_contigs = lectura_genoma(genoma)
    dicc_longitudes,num_contigs,long_total = cuenta_contigs(dicc_contigs)
    dicc_proporciones = proporcion_contigs(dicc_contigs,dicc_longitudes,long_total)
    dicc_reads,dicc_reads_reverso = reads(dicc_contigs,dicc_longitudes,dicc_proporciones,longitud_read,num_reads)
    crea_fastqz(dicc_reads,dicc_reads_reverso,longitud_read,name)
    
