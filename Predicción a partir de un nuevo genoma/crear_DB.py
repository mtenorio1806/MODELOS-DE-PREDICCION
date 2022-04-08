import sys, getopt, os, errno
import numpy as np
from Bio import SeqIO
#----------------------------Blast database------------------------------------#
def blast(dire,nom,out):
  rutablast='/vault2/soft/bio/prokka-master/binaries/linux/makeblastdb'
  print('\n\t----Creando base de datos blast----')
  os.system(rutablast+' -in '+str(dire)+' -dbtype \'prot\' -title '+nom+' -out '+out+'/'+nom)
#---------------------------crear multifastas de genes-------------------------#
def crearBD(carpeta):#la carpeta contiene una carpeta por cada modelo 
  genes=os.listdir(carpeta)#contenedora deantibioticos/  
  for g in genes:
    if os.path.isdir(carpeta+'/'+g):
      print('creando archivo multifasta de '+str(g))
      listain=os.listdir(carpeta+'/'+g) 
      listagenes={}    
      for l in listain:
        if l[-4:]=='alel' or l[-4:]=='orde':
          if listagenes.get(l[:-5],None)==None:
            listagenes[l[:-5]]=''
      print('lista de genes')
      print(listagenes)
      agregar=''
      for i in listagenes:
        alelos={}
        with open(carpeta+'/'+g+'/'+i+'.alel', mode='r') as arch_alelos:
            for line in arch_alelos:
              if line[:1]=='C':#centroide
                field = line.split()
                #print('alelo'+field[1])
                #print('alelo'+field[8])
                alelos[field[8]]=field[1]#se toma de las filas C el nombre del genoma y el numero del alelo
        for contig in SeqIO.parse(carpeta+'/'+g+'/'+i+'.orde', 'fasta'):
          if alelos.get(contig.id,None)!=None:#del archivo ordenado se extraen las secuencias de aminoacidos
            if str(alelos[contig.id])=='0':
              agregar+='>'+i+'\n'+contig.seq+'\n'
            else:
              agregar+='>'+i+'_'+alelos[contig.id]+'\n'+contig.seq+'\n'
        #print(i)
      print('creando multifasta de antibiotico')
      fasta=open(carpeta+'/'+g+'/'+g+'_multi.fasta','w')   
      fasta.writelines(agregar)
      fasta.close()
      #-crear base blast
      blast(carpeta+'/'+g+'/'+g+'_multi.fasta',str(g),carpeta+'/'+g+'/BlastDB_cp')   
        

  #os.chdir(carpeta)#rgi necesita un temporal de la carpeta de entrada"""
  #for k in nocorridos:
  #  os.system('rgi main -i'+carpeta+'/'+k+'.fasta -o '+salida+'/'+k+' -t contig -a DIAMOND -n 8 --include_loose --exclude_nudge -d wgs --low_quality --clean')
    
def help_message():
  print('\Crea multifasta y base de datos por gen con todos los alelos encontrados en el entrenamiento')
  print('programado y testeado con python 3.7........... (source activate py3.7)')
  print('Usage:')
  print('\tpython3 spadesdcto.py [options]')
  print('options:')
  print('\t-f direccion de la carpeta de entrada \n ')  
def main(a):
    f=0
    try:
        opts, args = getopt.getopt(a, "f:h")#   f carpeta con archivos alel y orde para hacer base de datos
        print(opts)
    except getopt.GetoptError:
        #help_message()
        print('entra en el primer help que de error')
        sys.exit(2)
    for opt, arg in opts:
        if opt == '-h':
            help_message()
            sys.exit()            
        elif opt == "-f": #, "--folderFastqs")
            fold= os.path.realpath(arg)#elimina los caracteres simbolicos de la direccion
            print('carpeta de entrada :  '+str(fold))
            f=1
        else:
            help_message()
            sys.exit()
    if f==1:
        crearBD(fold)
    else:
        help_message()
        print('no se ha seleccionado una carpeta para ejecutar el rgi')
        sys.exit()
if __name__ == "__main__":#__name__ variable de entorno que es main si se ejecuta el mismo programa desde el compilador, pero si el script es llamado de otro script tendra el valor de "codigo.py"
  if len(sys.argv[1:]) != 0:
     main(sys.argv[1:])
  else:
      help_message()