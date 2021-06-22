import sys
import numpy as np

def Alineamineto_Local(secuencia_1, secuencia_2, penalidad):

	for i in range(1,fila):
		for j in range(1, columna):
			valor = identicalMatch
			
			if secuencia_2[i-1] != secuencia_1[j-1]:
				valor = mismatch  

			Matriz[i][j] = max(Matriz[i-1][j-1] + valor, Matriz[i-1][j] + penalidad, Matriz[i][j-1] + penalidad, 0)


	print ()
	print (Matriz)		

	mayor = -1000
	i = 0
	j = 0
	for r in range(1,fila):
		for c in range(1, columna):
			if Matriz[r][c]>mayor:
				i = r
				j = c
				mayor = Matriz[r][c]

	secuencia_alineada_1= ""
	secuencia_alineada_2= ""

	while (i > 0 or j > 0):

		valor = identicalMatch
		if secuencia_2[i-1] != secuencia_1[j-1]:
			valor = mismatch  

		if (i>0 and j>0 and Matriz[i][j] == Matriz[i-1][j-1] + valor):		
			secuencia_alineada_1= secuencia_1[j-1] + secuencia_alineada_1
			secuencia_alineada_2= secuencia_2[i-1] + secuencia_alineada_2
			i = i-1
			j = j-1

		elif (i>0 and Matriz[i][j]==Matriz[i-1][j]+penalidad):
			secuencia_alineada_1= "-" + secuencia_alineada_2	
			secuencia_alineada_2= secuencia_2[i-1] + secuencia_alineada_2
			i = i-1

		else:
			secuencia_alineada_2= "-" + secuencia_alineada_2
			secuencia_alineada_1= secuencia_1[j-1] + secuencia_alineada_1
			j = j-1	

		if Matriz[i][j]==0:
			break	


	print ()
	print("Secuencia Alineada del Archivo ", f1, " : \n")		
	print (secuencia_alineada_1)
	print ()
	print("Secuencia Alineada del Archivo ", f2, " : \n")
	print (secuencia_alineada_2)
	print("\n")
	print()			


def Obtener_Secuencia(file):
	next(file)  
	secuencia = ""

	for linea in file: 
		if len(linea.strip())!=0:
			linea=linea.rstrip('\n')
			secuencia = secuencia + linea

	file.close()

	return secuencia		


if __name__ == "__main__":

	archivo_1 = sys.argv[1]
	archivo_2 = sys.argv[2]	
	
	print("\n\tAlineamiento Local - (Data : FASTA)\n")
	penalidad = int(input("Penalidad : "))

	identicalMatch = int(input("Identical Match : "))
	mismatch = int(input("MisMatch : "))

	f1 = open(archivo_1, "r")
	f2 = open(archivo_2, "r")

	secuencia_1 = Obtener_Secuencia(f1)
	secuencia_2 = Obtener_Secuencia(f2)

	columna = len(secuencia_1)+1
	fila = len(secuencia_2)+1

	print(fila)
	print(columna)

	Matriz = np.zeros([fila, columna], dtype=int) 

	Alineamineto_Local(secuencia_1, secuencia_2, penalidad)	