import sys
import numpy as np

def Alineamiento_Global(secuencia_1, secuencia_2,penalidad):
	columna = len(secuencia_1)+1
	fila = len(secuencia_2)+1
	F = np.zeros([fila, columna], dtype=int) 
	# Agregamos 0's
	for i in range(1,columna):
		F[0][i] = i*penalidad

	for i in range(1,fila):
		F[i][0] = i*penalidad

	for i in range(1,fila):
		for j in range(1, columna):
			valor = identicalMatch
			if secuencia_2[i-1] != secuencia_1[j-1]:
				valor = mismatch  
			F[i][j]	= max(F[i-1][j-1] + valor, F[i-1][j] + penalidad, F[i][j-1] + penalidad)
	print("Matriz:\n")
	print (F)

	i = fila-1
	j = columna-1		

	secuencia_alineada_1 = ""
	secuencia_alineada_2 = ""

	while (i > 0 or j > 0): # alineamento
		valor = identicalMatch
		if secuencia_2[i-1] != secuencia_1[j-1]:
			valor = mismatch  
		if (i>0 and j>0 and F[i][j] == F[i-1][j-1] + valor):		
			secuencia_alineada_1 = secuencia_1[j-1] + secuencia_alineada_1
			secuencia_alineada_2 = secuencia_2[i-1] + secuencia_alineada_2
			i = i-1
			j = j-1
		elif (i>0 and F[i][j]==F[i-1][j]+penalidad):
			secuencia_alineada_1 = "-" + secuencia_alineada_2	
			secuencia_alineada_2 = secuencia_2[i-1] + secuencia_alineada_2
			i = i-1
		else:
			secuencia_alineada_2 = "-" + secuencia_alineada_2
			secuencia_alineada_1 = secuencia_1[j-1] + secuencia_alineada_1
			j = j-1							
	print ()
	print("Secuencia Alineada 1 :\n")		
	print (secuencia_alineada_1)
	print ()
	print("Secuencia Alineada 2:\n")
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
	file1 = sys.argv[1]
	file2 = sys.argv[2]	
	penalidad = int(sys.argv[3])
	identicalMatch = 2
	mismatch = -2
	f1 = open(file1, "r")
	f2 = open(file2, "r")
	secuencia_1 = Obtener_Secuencia(f1)
	secuencia_2 = Obtener_Secuencia(f2)
	Alineamiento_Global(secuencia_1, secuencia_2,penalidad)