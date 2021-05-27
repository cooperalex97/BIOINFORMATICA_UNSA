import sys
import numpy as np

def Matriz_Sustitucion(file):
	next(file)  
	r_c = {}
	s = []
	i = 0
	for line in file: 
		if len(line.strip())!=0:
			line=line.rstrip('\n')

			r_c[line.split('\t')[0:1][0]] = i 
			i = i+1
			s.append(list(map(int,line.split('\t')[1:])))
	file.close()
	return r_c, s

def Obtener_Secuencias(F, i, j, secuencia_alineada_1 = "", secuencia_alineada_2 = ""):
	if F[i][j][1]==0:
		print (secuencia_alineada_1)
		print (secuencia_alineada_2)
		Obtener_Puntaje(secuencia_alineada_1, secuencia_alineada_2)
		print ()	
		return
			
	if len(F[i][j][1])>1:
		direcciones = F[i][j][1]
		for n in range(len(direcciones)):	
			F[i][j][1] = direcciones[n]
			Obtener_Secuencias(F, i, j, secuencia_alineada_1, secuencia_alineada_2)			

	else:	
		if F[i][j][1] == 'D':		
			secuencia_alineada_1 = secuencia_1[j-1] + secuencia_alineada_1
			secuencia_alineada_2 = secuencia_2[i-1] + secuencia_alineada_2
			i = i-1
			j = j-1

		elif F[i][j][1]=='U':
			secuencia_alineada_1 = "-" + secuencia_alineada_2	
			secuencia_alineada_2 = secuencia_2[i-1] + secuencia_alineada_2
			i = i-1
		else:
			secuencia_alineada_2 = "-" + secuencia_alineada_2
			secuencia_alineada_1 = secuencia_1[j-1] + secuencia_alineada_1
			j = j-1

		Obtener_Secuencias(F, i, j, secuencia_alineada_1, secuencia_alineada_2)			


def Obtener_Puntaje(secuencia_1, secuencia_2):
	puntaje = 0
	for i in range(len(secuencia_1)):
		if secuencia_1[i] == secuencia_2[i]:
			print (s[r_c[secuencia_2[i]]][r_c[secuencia_1[i]]], end='', sep='')
			puntaje = puntaje + s[r_c[secuencia_2[i]]][r_c[secuencia_1[i]]]
		else:
			if secuencia_1[i]=='-' or secuencia_2[i]=='-':
				print ('(',d,')', end='', sep='')
				puntaje = puntaje + d
			else:
				print ('(',s[r_c[secuencia_2[i]]][r_c[secuencia_1[i]]] ,')', end='', sep='')		
				puntaje = puntaje + s[r_c[secuencia_2[i]]][r_c[secuencia_1[i]]]

		print ('+', end='', sep='')

	print("\b",end="") 		
	print ('=', puntaje)


def Alineamiento_Global(F, i, j):
	direcciones = ''

	diag = F[i-1][j-1][0] + s[r_c[secuencia_2[i-1]]][r_c[secuencia_1[j-1]]]
	up = F[i-1][j][0] + d
	left = F[i][j-1][0] + d

	F[i][j][0]	= max(diag, up, left)
	
	if F[i][j][0]==left:	
		direcciones = direcciones + 'L'

	if F[i][j][0]==diag:
		direcciones = direcciones + 'D'

	if F[i][j][0]==up:	
		direcciones = direcciones + 'U'	
	
	F[i][j][1] = direcciones

	if i==fila-1 and j==columna-1:
		Guardar_Resultado(F)
		Obtener_Secuencias(F, i, j)
		return

	elif j<columna-1:
		Alineamiento_Global(F, i ,j+1)	
	else:	
		Alineamiento_Global(F, i+1 ,1)	


def Guardar_Resultado(F):
	f = open(archivo_final, "w")

	f.write('              ')
	for i in range(columna-1):
		f.write(secuencia_1[i]+'        ')

	f.write('\n')	

	secuencia_2_t = ' '+ secuencia_2

	for r in range(fila):
		f.write(secuencia_2_t[r])
		for c in range(columna):			
			f.write(' ( '+ '{}'.format(F[r][c][0])+' '+'{}'.format(F[r][c][1])+') ')

		f.write('\n')
	f.close()	

	f = open(archivo_final, "r")
	print(f.read())

 
if __name__ == "__main__":
	filename = sys.argv[1]	
	d = int(sys.argv[2])
	archivo_final = "resultado.txt"
	
	file = open(filename, "r")

	r_c, s = Matriz_Sustitucion(file)	

	secuencia_1 = "AAG"
	secuencia_2 = "AGC"		

	# columnaa y fila  de 0's al inicio de la matriz
	columna = len(secuencia_1)+1
	fila = len(secuencia_2)+1

	F = np.zeros([fila, columna], dtype='i,O') 

	# Agregando 0's
	for i in range(1,columna):
		F[0][i][0] = i*d
		F[0][i][1] = 'L'

	for i in range(1,fila):
		F[i][0][0] = i*d
		F[i][0][1] = 'U'
		
	Alineamiento_Global(F, 1, 1)	
