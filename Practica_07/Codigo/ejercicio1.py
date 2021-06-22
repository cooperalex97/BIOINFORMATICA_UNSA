import sys
import numpy as np

def Matriz_Sustitucion(archivo):
	next(archivo)  
	r_c = {}
	s = []
	i = 0

	for line in archivo: 
		if len(line.strip())!=0:
			line=line.rstrip('\n')

			r_c[line.split('\t')[0:1][0]] = i 
			i = i+1

			s.append(list(map(int,line.split('\t')[1:])))

	archivo.close()

	return r_c, s


def Guardar_Resultado(F):
	f = open(resultadoArchivo, "w")

	f.write('            ')
	for i in range(columna-1):
		f.write(sequencia1[i]+'        ')

	f.write('\n')	

	seq2_t = ' '+ sequencia2

	for r in range(fila):
		f.write(seq2_t[r])
		for c in range(columna):			
			f.write(' ( '+ '{}'.format(F[r][c][0])+' '+'{}'.format(F[r][c][1])+') ')

		f.write('\n')
	f.close()	

	f = open(resultadoArchivo, "r")
	print(f.read())	

def Obtener_Secuencias(F, i, j, alineamiento_sequencia1 = "", alineamiento_sequencia2 = ""):

	if F[i][j][0]==0:
		print ()		
		print (alineamiento_sequencia1)
		print (alineamiento_sequencia2)	
		return
			
	
	if i > 0 or j > 0:			
		
		if (i>0 and j>0 and F[i][j][1] == 'D'):		
			alineamiento_sequencia1 = sequencia1[j-1] + alineamiento_sequencia1
			alineamiento_sequencia2 = sequencia2[i-1] + alineamiento_sequencia2
			i = i-1
			j = j-1

		elif (i>0 and F[i][j][1]=='U'):
			alineamiento_sequencia1 = "-" + alineamiento_sequencia2	
			alineamiento_sequencia2 = sequencia2[i-1] + alineamiento_sequencia2
			i = i-1

		else:
			alineamiento_sequencia2 = "-" + alineamiento_sequencia2
			alineamiento_sequencia1 = sequencia1[j-1] + alineamiento_sequencia1
			j = j-1

		Obtener_Secuencias(F, i, j, alineamiento_sequencia1, alineamiento_sequencia2)			


def Alineamiento_Local(F, i, j, maximo_total):

	diag = F[i-1][j-1][0] + s[r_c[sequencia2[i-1]]][r_c[sequencia1[j-1]]]
	up = F[i-1][j][0] + d
	left = F[i][j-1][0] + d
	diag2 = 0

	entrada_max = [diag,up,left,diag2]
	maximo = max(entrada_max)

	if maximo > maximo_total:
                maximo_total = maximo

	F[i][j][0]  = maximo
	
	if F[i][j][0]==diag:
		F[i][j][1] = 'D'	

	if F[i][j][0]==up:	
		F[i][j][1] = 'U'
	
	if F[i][j][0]==left:	
		F[i][j][1] = 'L'

	if i==fila-1 and j==columna-1:
		Guardar_Resultado(F)

		major = -1000
		s_i = 0
		s_j = 0 
		for r in range(1,fila):
			for c in range(1, columna):
				if F[r][c][0]>=major:
					s_i = i
					s_j = j
					i = r
					j = c
					major = F[r][c][0]

		Obtener_Secuencias(F, i, j)
		Obtener_Secuencias(F, s_i, s_j)
		print("Maximo Total", maximo_total)
		return

	if j<columna-1:
		Alineamiento_Local(F, i ,j+1, maximo_total)	
	else:	
		Alineamiento_Local(F, i+1 ,1, maximo_total)


if __name__ == "__main__":

	archivoNombre = sys.argv[1]	
	d = int(sys.argv[2])

	resultadoArchivo = "resultado.txt"
	
	archivo = open(archivoNombre, "r")

	
	r_c, s = Matriz_Sustitucion(archivo)	

	sequencia1 = "AAG"
	sequencia2 = "AGC"	

	
	columna = len(sequencia1)+1
	fila = len(sequencia2)+1

	F = np.zeros([fila, columna], dtype='i,O')

	maximo_total = 0

	Alineamiento_Local(F, 1, 1, maximo_total)

	
