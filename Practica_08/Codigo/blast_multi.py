import numpy as np
import operator
import multiprocessing 
import time
import math 

"""
Esta funcion divide el calculo de Blast en diferentes procesos
OJO : el numero de proceso esta en funcion a la cantidad de secuencias en la Base de Datos Principal
Por ejemplo en la base de datos covid_database.txt presenta  el siguiente conjunto de secuencias fasta: 

	1. >MN908947.3 Severe acute respiratory syndrome coronavirus 2 isolate Wuhan-Hu-1, complete genome
	2. >MN938384.1 Severe acute respiratory syndrome coronavirus 2 isolate 2019-nCoV_HKU-SZ-002a_2020, complete genome
	3. >JQ765563.1 Human coronavirus NL63 strain NL63/DEN/2009/9, complete genome
	4. >JQ765564.1 Human coronavirus NL63 strain NL63/DEN/2009/14, complete genome

 -> #pros = #sec

 Entonces BLAST se dividira para evaluar a dichas secuencia en un proceso diferente. 
 Al final une los resultados con el p.join() 
 Recupera los datos del diccionario final
 Luego itera para mostrar los alineamientos resultantes.

"""
def Multiprocesamiento_BLAST(sec_query, secuencias_base_datos, k, match, mismatch, gap, extension_threshold, hssp_threshold, nombres):
	procesos = []
	manager = multiprocessing.Manager()
	retornar_diccionario = manager.dict()

	for i in nombres:
		p = multiprocessing.Process(target=Blast_xSecuencia, args=(sec_query, secuencias_base_datos, k, match, mismatch, gap, extension_threshold, hssp_threshold, nombres, i, retornar_diccionario))

		procesos.append(p)

	for i in procesos:
		i.start()

	for p in procesos:
		p.join()

	total_alineamientos = 0
	salida_completa = ""

	for nombre in nombres:

		salida = retornar_diccionario.get(nombre)
		total_alineamientos += salida.count("Alineamiento")
		salida_completa += salida

	i = 0
	contador = 1
	indice = salida_completa.find("Alineamiento",i)
	while(indice != -1):
		ind2 = salida_completa.find("**",indice)
		salida_completa = salida_completa[:indice] + "Alineamiento * " + str(contador) + " " + salida_completa[ind2:]
		i = indice + 1
		indice = salida_completa.find("Alineamiento",i)
		contador+=1

	print(salida_completa)
	print("\t\t###########################################################################################")
	print("\t\t#                                                                                         #")
	print("\t\t#\t\tNumero de Procesos                       : " + str(len(procesos))+"                              #")
	print("\t\t#\t\tTotal de Alineamientos en la BD Completa : " + str(total_alineamientos)+"                              #")


# BLAST entre el query  y la secuencia especificada en la BD 
def Blast_xSecuencia(sec_query, secuencias_base_datos, k, match, mismatch, gap, extension_threshold, hssp_threshold, names, i, return_dict):
	nombre_secuencia_BD = i  # es su ID de la secuencia en la BD
	secuencia_BD = secuencias_base_datos[i]
	#secuencia_palabras = []

	tamanio_secuencia_DB = len(secuencia_BD)
	tamanio_query = len(sec_query)

	div_inicial_matriz_db = Dividir_Secuencia_DB(secuencia_BD, tamanio_secuencia_DB, k)

	fragmento_inicial_matriz_query = Dividir_Secuencia_Query(sec_query, tamanio_query, k, hssp_threshold)

	lista_pares_hssp = Encontrar_Pares_HSSP( fragmento_inicial_matriz_query, div_inicial_matriz_db)

	seed_output = Encontrar_Seed(sec_query, secuencia_BD, lista_pares_hssp, k)
	max_length = [0]
	ali_count = [1]
	output = []
	alineamientos_totales = [0]
	score_frec = {}
	for seed in seed_output:
		SMalignment(seed[3], seed[2], match, mismatch, gap, extension_threshold, nombre_secuencia_BD, seed[0], seed[1], ali_count, max_length, output, alineamientos_totales, score_frec)
	
	alineamientos_totales2 = [tamanio_query*tamanio_secuencia_DB - (tamanio_query*(tamanio_query+1))//2]

	temp_score_arr = []
	for j in score_frec:
		temp_score_arr.append([j, score_frec[j]])
	temp_score_arr = sorted(temp_score_arr, key = lambda x:x[0])
	l = len(temp_score_arr)
	s = 0
	for index in range(l-1,-1,-1):
		s+=temp_score_arr[index][1]
		score_frec[temp_score_arr[index][0]] = s
	l = len(output)


	for k in range(l):
		out = output[k]

		s = int(out[5][out[5].find("= ")+2:-1])
		# print(score_frec[s], s)
		nume = score_frec[s]
		pvalue = nume/alineamientos_totales2[0]
		evalue = alineamientos_totales[0]*pvalue
		bitscore = math.log((alineamientos_totales[0])/evalue)
		out.append("Bitscore = " + str(bitscore) + "\n")
		out.append("p-value = " + str(pvalue) + "\n")
		out.append("e-value = " + str(evalue) + "\n")

	output_string = ""


	for arr in range(l):

		for ind in range(len(output[arr])):
			output_string+=output[arr][ind]

	print(i)
	return_dict[i] = output_string

	return return_dict



#######################################################################################################################################
#
#                          FRAGMENTACION Y PREPROCESAMIENTO
#
#######################################################################################################################################
"""
Esta función se utiliza para obtener la representación numérica de la secuencia de nucleótidos.
"""
def Convertir_Palabra_Numero(palabra):
    tmp = []
    valores = {'A':1,'C':2,'G':3,'T':4}
    for w in palabra:
        tmp.append(valores[w])
    return tmp


"""
Aqui se transforma, evalúa el índice de la representación numérica de la secuencia de nucleótidos
"""
def Convertir_Palabra_Indice(palabra,len_palabra):
    tmp = 0
    valor_palabra = Convertir_Palabra_Numero(palabra)
    for i,v in enumerate(valor_palabra):
        tmp += (v-1)*(4**(len_palabra-i-1))
    return tmp   


"""
Rompe la secuencia del genoma de la BD  a sus palabras alargadas en 'W' mapeándolas sus posiciones iniciales 
en la secuencia del genoma de la BD. 
Devuelve una matriz ordenada, osea sobre la base de las representaciones numéricas de las palabras. 
A esta resultado le podemos aplicar Busqueda Binaria 

Ya Chicas si esta facil de entender xd Gaaaaaaaaaaaaaaa
OJO :  k-mer -> se realiza en funcion al valor de k 

     Ej 3-mer -> ACG  ->  k = 3

"""

def Dividir_Secuencia_DB(secuencia_BD, tamanio_secuencia_DB, k):
	diccionario_incial_dividido = {}

	i = 0
	while( i+k-1 < tamanio_secuencia_DB ):
		palabra = secuencia_BD[i:i+w]
		indice_palabra = Convertir_Palabra_Indice(palabra, k)

		if(diccionario_incial_dividido.get(indice_palabra) != None):
			diccionario_incial_dividido[indice_palabra].append(i)
		else:
			diccionario_incial_dividido[indice_palabra] = [i]
		i+=1

	array_division_inicial = []
	for j in diccionario_incial_dividido:
		array_division_inicial.append([j, diccionario_incial_dividido[j]])
	array_division_inicial = sorted(array_division_inicial, key = lambda x:x[0])

	return array_division_inicial



"""
Aqui la funcion hace lo mismo que la anterior pero para la secuencia del Query. 


"""
def Dividir_Secuencia_Query(sec_query, tamanio_query, k, hssp_threshold):
	match = 5
	mismatch = -4
	index_arr = []
	i = 0

	fragmentos_query = []
	while( i+w-1 < tamanio_query ):
		palabra = sec_query[i:i+w]
		vecinos = set([])
		Obtener_Vecinos(palabra, k, match, mismatch, 0, match*w, vecinos, hssp_threshold)
		vecinos = list(vecinos)
		for j in vecinos:
			fragmentos_query.append([ Convertir_Palabra_Indice(j,w),i])
		i+=1
	return fragmentos_query


"""

Esta función encuentra todos los vecinos de la palabra de consulta (query) con matching_score > hssp_threshold
OJO:   El HSSP Threshold se ingresa al inicio del programa. 

"""
def Obtener_Vecinos(palabra, w, match, mismatch, i, score, vecinos, hssp_threshold):
	if( score < hssp_threshold  or  i == w):
		return

	vecinos.add(palabra)

	letras = ["A","T","G","C"]

	for r in letras:
		if(r!=palabra[i]):
			palabra_vecina = palabra[:i] + r + palabra[i+1:]
			Obtener_Vecinos(palabra_vecina, w, match, mismatch, i+1, score-match+mismatch, vecinos, hssp_threshold)

	Obtener_Vecinos(palabra, w, match, mismatch, i+1, score, vecinos, hssp_threshold)



##########################################################################################################################################
#
#                        Alineación de palabras y Búsqueda de combinaciones  HSSP coincidentes 
#
##########################################################################################################################################


"""
Búsqueda binaria en la secuencia del genoma de la BD para el k-mer  de consulta.

"""
def Busqueda_Binaria_Secuencia(palabra, div_inicial_matriz_db, izquierdo, derecho):
	if derecho >= izquierdo: # Caso base

		mid = izquierdo + (derecho - izquierdo) // 2

		if div_inicial_matriz_db[mid][0] == palabra[0]:
			return div_inicial_matriz_db[mid][1]

		elif div_inicial_matriz_db[mid][0] > palabra[0]:
			return Busqueda_Binaria_Secuencia(palabra, div_inicial_matriz_db, izquierdo, mid-1)
		else:
			return Busqueda_Binaria_Secuencia(palabra, div_inicial_matriz_db, mid+1, derecho)
	else:
		return []

'''
Aqui el objetivo es encontrar los  pares de HSSP entre 
el query y las palabras de la secuencia de la BD
'''
def Encontrar_Pares_HSSP( fragmento_inicial_matriz_query, div_inicial_matriz_db):
	lista_pares_hssp = []
	for palabra in fragmento_inicial_matriz_query:
		matching_posiciones_iniciales = Busqueda_Binaria_Secuencia(palabra, div_inicial_matriz_db, 0, len(div_inicial_matriz_db)-1)

		if(matching_posiciones_iniciales != []):
			for i in matching_posiciones_iniciales:
				lista_pares_hssp.append([palabra[1], i])
	return lista_pares_hssp

'''
Aqui encontramos todos los "Seeds" de longitud exacta de coincidencias con espacios iguales a la longitud general
El resultado es el lapso de la posible unión de seeds en forma de lista 2D.
'''
def Encontrar_Seed(secuencia_query, secuencia_BD, lista_inicial, k):
    nueva_lista = []
    for i in lista_inicial:
        for j in range(1, k):
            n = [i[0]+j, i[1]+j]
            if(n not in nueva_lista and n not in lista_inicial):
                nueva_lista.append(n)

    coord_list = lista_inicial + nueva_lista

    Pos_Matriz = np.full([len(secuencia_query), len(secuencia_BD)], 0, dtype=int)
    for coord_actual in coord_list:
        Pos_Matriz[coord_actual[0]][coord_actual[1]] = 1

    visitado = np.full([len(secuencia_query), len(secuencia_BD)], 0, dtype=int)

    seed_collection = []
    rango_seed = []
    total_seeds = 0

    for coord_actual in lista_inicial:
        i, j = coord_actual
        seed_actual = []
        if(Pos_Matriz[i][j] == 1 and visitado[i][j] == 0):
            dfs(Pos_Matriz, i, j, visitado, seed_actual, secuencia_BD)
            total_seeds += 1
        # print(seed_actual)
        if(len(seed_actual) > 0):
            x_min = seed_actual[0][0]
            x_max = seed_actual[len(seed_actual)-1][0]
            y_min = seed_actual[0][1]
            y_max = seed_actual[len(seed_actual)-1][1]

            sub_secuencia_query = secuencia_query[x_min:x_max+1]
            sub_secuencia_BD = secuencia_BD[y_min:y_max+1]
            seed_collection.append(seed_actual)
            rango_seed.append([[x_min, x_max], [y_min, y_max],
                               sub_secuencia_query, sub_secuencia_BD])

    # print((seed_collection))
    # print(rango_seed)
    return rango_seed		

'''
Función para verificar si los vecinos son los correctos  para DFS
'''
def Check_Vecinos(Pos_Matriz, i, j, visitado, secuencia_BD):
    return (i >= 0 and i < len(secuencia_query) and j >= 0 and j < len(secuencia_BD) and (not visitado[i][j] and Pos_Matriz[i][j] == 1))

'''
dfs para encontrar los seeds ( semillas )
'''
def dfs(Pos_Matriz, i, j, visitado, seed_actual, secuencia_BD):
	# comprobando los bloques de la derecha , abajo y diagonalmente a la derecha
    siguiente_fila = [1, 0, 1]
    siguiente_columna = [0, 1, 1]
    visitado[i][j] = 1
    seed_actual.append([i, j])

    for k in range(3):
        if(Check_Vecinos(Pos_Matriz, i + siguiente_fila[k], j + siguiente_columna[k], visitado, secuencia_BD)):
            dfs(Pos_Matriz, i + siguiente_fila[k], j + siguiente_columna[k], visitado, seed_actual, secuencia_BD)


'''
Alineación de Smith-Waterman en las regiones donde se dan uniones de las palabras coincidentes vecinas
'''
def SMalignment(secuencia_1, secuencia_2, match, mismatch, g, extension_threshold, nombre_secuencia_BD, rango_query, rango_DB, ali_count, max_length, output, alineamientos_totales, score_frec):
	m = len(secuencia_1)
	n = len(secuencia_2)
	matriz = []

	for i in range(0, m):
		tmp = []
		for j in range(0, n):
			tmp.append(0)
		matriz.append(tmp)

	for sii in range(0, m):
		matriz[sii][0] = sii*g
	for sjj in range(0, n):
		matriz[0][sjj] = sjj*g

	maxx_score = 0
	punto_inicio = []
	punto_inicio_visitado = {}

	for siii in range(1, m):
		for sjjj in range(1, n):
			matriz[siii][sjjj] = max(matriz[siii-1][sjjj] + g, matriz[siii - 1][sjjj - 1] + Comparar_Bases(secuencia_1,secuencia_2,siii, sjjj, match, mismatch), matriz[siii][sjjj-1] + g)
			if(score_frec.get(matriz[siii][sjjj]) == None):
				score_frec[matriz[siii][sjjj]] = 1
			else:
				score_frec[matriz[siii][sjjj]] += 1
			alineamientos_totales[0] += 1
			if(matriz[siii][sjjj] >= extension_threshold):
				tup = (siii,sjjj, matriz[siii][sjjj])
				punto_inicio.append(tup)
				punto_inicio_visitado[(tup[0], tup[1])] = False


	punto_inicio = sorted(punto_inicio, key = lambda x:(x[0], x[1]))
	l = len(punto_inicio)

	for j in range(l-1,-1,-1):
		i = punto_inicio[j]
		if(not(punto_inicio_visitado[(i[0],i[1])])):
			temp = []
			a = backtracking(i[0], i[1], secuencia_1, secuencia_2, matriz, punto_inicio_visitado)
			temp.append("\n")
			temp.append("\n")
			temp.append("Alineamiento * " + str(ali_count[0]) + " **********************************************************************************************************************\n" )
			ali_count[0]+=1
			disp = Display(a[0], a[1], nombre_secuencia_BD)
			temp.append(disp + "\n")
			#temp.append("Detalles del Proceso :\n")
			temp.append("Sequence Identity : " + str(a[2]) + "\n")
			temp.append("Score = " + str(i[2]) + "\n")
			temp.append("Dtabase sequence range: " + str(rango_DB[0] + a[3]) + " - " + str(rango_DB[0] + a[5]) + "\n")
			temp.append("Query sequence range: " + str(rango_query[0] + a[4]) + " - " + str(rango_query[0] + a[6]) + "\n")
			temp.append("Length :" + str(len(a[0])) + "\n")
			max_length[0] = max(max_length[0], len(a[0]))
			output.append(temp)


def Comparar_Bases(secuencia_1,secuencia_2,i,j, match, mismatch):
    if secuencia_1[i] == secuencia_2[j]:
        return match
    else:
        return mismatch

'''
Backtracking en SW de los puntos que tienen mayor score que el threshold
'''
def backtracking(m, n, sec_1, sec_2, matriz, punto_inicio_visitado):
	secu_1 = [sec_1[m-1]]
	secu_2 = [sec_2[n-1]]

	sec_final_1 = m-1
	sec_final_2 = n-1
	sec_inicial_1 = m-1
	sec_inicial_2 = n-1
	while m > 1 and n > 1:
		if max(matriz[m-1][n-2], matriz[m-2][n-2], matriz[m-2][n-1]) == matriz[m-2][n-2]:
			m -= 1
			n -= 1
			if(punto_inicio_visitado.get((m,n)) != None):
				punto_inicio_visitado[(m,n)] = True

			sec_inicial_1 = m
			sec_inicial_2 = n

			secu_1.append(sec_1[m-1])
			secu_2.append(sec_2[n-1])
		elif max(matriz[m-1][n-2], matriz[m-2][n-2], matriz[m-2][n-1]) == matriz[m-1][n-2]:
			n -= 1
			if(punto_inicio_visitado.get((m,n)) != None):
				punto_inicio_visitado[(m,n)] = True
			sec_inicial_2 = n
			secu_1.append('-')
			secu_2.append(sec_2[n-1])
		else:
			m -= 1
			if(punto_inicio_visitado.get((m,n)) != None):
				punto_inicio_visitado[(m,n)] = True
			sec_inicial_1 = m
			secu_1.append(sec_1[m-1])
			secu_2.append('-')
	secu_1.reverse()
	secu_2.reverse()
	alinear_sec_1 = ''.join(secu_1)
	alinear_sec_2 = ''.join(secu_2)

	alinear_score = 0
	for k in range(0, len(alinear_sec_1)):
		if alinear_sec_1[k] == alinear_sec_2[k]:
			alinear_score += 1
	alinear_score = float(alinear_score)/len(alinear_sec_1)

	return [alinear_sec_1, alinear_sec_2, alinear_score, sec_inicial_1, sec_inicial_2, sec_final_1, sec_final_2]


##################################################################################################################################
#
#                 MOSTRANDO  EL  ALINEAMIENTO  BLAST  
#
#################################################################################################################################

def Display(secuencia_01, secuencia_02, nombre_secuencia_BD):
    le = 0
    l = len(secuencia_01)
    output = ""
    output+="\n"
    output += "Secuencia 1 = " + nombre_secuencia_BD + "\n"
    output += "Secuencia 2 = Secuencia del Query ( Consulta )" + "\n"
    output+="\n"
    output += '\t\t\t\tSecuencia_1 : '
    for a in list(secuencia_01):
        output += a
    output += "\n"
    output += "\n"
    output += '\t\t\t\t              '
    for k in range(l):
        if secuencia_01[k] == secuencia_02[k]:
            output += '|'
        elif(secuencia_01[k] !="-" and  secuencia_02[k]!="-"):
        	output += ':'
        else:
            output += ' '
    output += "\n"
    output += "\n"
    output += '\t\t\t\tSecuencia_2 : '
    for b in list(secuencia_02):
        output += b
    output += "\n"

    return output


def Base_Datos():
	archivo_inicial = open("covid_database.txt","r")
	secuencias_base_datos = {}
	secuencia = ""
	nombre = ""
	nombres = []
	for linea in archivo_inicial:
		if(linea[0] == ">"):
			if(nombre!=""):
				secuencias_base_datos[nombre] = secuencia
			secuencia = ""
			nombre = linea[:-1]
			nombres.append(nombre)
		else:
			secuencia += linea[:-1]
	secuencias_base_datos[nombre] = secuencia
	return secuencias_base_datos,nombres

def Presentacion():
	print("\t\t\t##############################################")
	print("\t\t\t#                                            #")
	print("\t\t\t#\t\t *** BLAST ***               #")
	print("\t\t\t#                                            #")
	print("\t\t\t##############################################\n")
	print(">>> Ejemplos de Querys (Consulta) :\n ")
	print("\t Q1 : CGTGAGTCAGCTATTGAACTGGCCGCGCAATGGAAGAGTTGTTAATCCGCAAAATCTGGCAAC\n")
	print("\t Q2 : ATTAAAGGTTTATACCTTCCCAGGTAACAAACCAACCAACTTTCGATCTCTTGTAGATCTGTTCTCTAAA\n")
	print("\t Q3 : CGAACTTTAAAATCTGTGTGGCTGTCACTCGGCTGCATGCTTAGTGCACTCACGCAGTATAATTAATAAC\n")
	print("\n\tPuede colocar q# o Q#")
	print("\n\tEjemplo q1 o Q1")
	print("\n\tSi desea otro query solo presione enter e ingrese su Query\n")

def Entrada_Datos():
	q1 = "CGTGAGTCAGCTATTGAACTGGCCGCGCAATGGAAGAGTTGTTAATCCGCAAAATCTGGCAAC"
	q2 = "ATTAAAGGTTTATACCTTCCCAGGTAACAAACCAACCAACTTTCGATCTCTTGTAGATCTGTTCTCTAAA"
	q3 = "CGAACTTTAAAATCTGTGTGGCTGTCACTCGGCTGCATGCTTAGTGCACTCACGCAGTATAATTAATAAC"
	secuencia_query = input("\tIngrese la Secuencia del Query : ")
	if secuencia_query == 'q1' or secuencia_query=='Q1':
		secuencia_query = q1
	elif secuencia_query == 'q2' or secuencia_query=='Q2':
		secuencia_query = q2
	elif secuencia_query == 'q3' or secuencia_query=='Q3':
		secuencia_query = q3
	else :
		secuencia_query = input("\tIngrese su Secuencia del Query : ")

	print("\n>>> Entrada de Datos\n")
	w = int(input("\t>> Ingrese la longitud del k-mer (Ej - k = 11) : "))
	hssp_threshold = int(input("\t>> Ingrese HSSP Threshold (Ej - 45)            : ")) ## 50 permite coincidencias exactas, 45 permite 1 mismatch
	extension_threshold = int(input("\t>> Ingrese Extension Threshold (Ej = 60)       : "))
	match = int(input("\t>> Ingrese Match                               : ")) #5
	mismatch = int(input("\t>> Ingrese MissMatch                           : ")) #-4
	gap = int(input("\t>> Ingrese GAP                                 : "))#-3

	return secuencia_query,w,hssp_threshold,extension_threshold,match,mismatch,gap

def Tiempo_Ejecucion():
	tiempo_inicio = time.time()
	Multiprocesamiento_BLAST(secuencia_query, secuencias_bd, w, match, mismatch, gap, extension_threshold, hssp_threshold, nombres)
	tiempo_final = time.time()
	return tiempo_final - tiempo_inicio


if __name__=="__main__":
	secuencias_bd, nombres = Base_Datos()
	Presentacion()
	secuencia_query,w,hssp_threshold,extension_threshold,match,mismatch,gap = Entrada_Datos()
	T_ejecucion = Tiempo_Ejecucion()
	print("\t\t#\t\tTiempo de Ejecucion                      : " + str(T_ejecucion) +" Segundos    #")
	print("\t\t#                                                                                         #")
	print("\t\t###########################################################################################")
