# Autor: Luis Peláez
# Inicio 02/05/2020
# Proyecto para el calculo del punto isoelectrico


from sys import exit, argv
#from os import system
#system("cls")
#system("clear")


# acid indica que se ioniza después de su pk
# basic indica que se ioniza antes de su pk
aminoacidos ={
	"G": ["GLY", 2.34, 9.6],
	"A": ["ALA", 2.34, 9.69],
	"P": ["PRO", 1.99, 10.96],
	"V": ["VAL", 2.32, 9.62],
	"L": ["LEU", 2.36, 9.60],
	"I": ["ILE", 2.36, 9.68],
	"M": ["MET", 2.28, 9.21],
	"F": ["PHE", 1.83, 9.13],
	"Y": ["TYR", 2.2, 9.11, 10.07, "acid"],
	"W": ["TRP", 2.38, 9.39],
	"S": ["SER", 2.21, 9.15],
	"T": ["THR", 2.11, 9.62],
	"C": ["CYS", 1.96, 10.28, 8.18, "acid"],
	"N": ["ASN", 2.02, 8.8],
	"Q": ["GLN", 2.17, 9.13],
	"D": ["ASP", 1.88, 9.60, 3.65, "acid"],
	"E": ["GLU", 2.19, 9.67, 4.25, "acid"],
	"K": ["LYS", 2.18, 8.95, 10.53, "basic"],
	"H": ["HIS", 1.82, 9.17, 6.0, "basic"],
	"R": ["ARG", 2.17, 9.04, 12.48, "basic"]
}

def ayuda(h):
	if h == "-h".lower():
		print('Para ejecutar:')
		print('user@user: python isopoint.py')
	else:
		print('Opcion inválida!')


def presentacion():
	print("************************************************************")
	print("Isopoint es una aplicación creada con el fin de facilitar el\n"
		  "cálculo del punto isoelétrico de péptidos y proteínas.")
	print("************************************************************\n")


def ingreso_cadena():
	cadena = input("Ingresa la secuencia de aminoácidos:\n")
	cadena = cadena.upper()
	verificar_peptido(cadena)
	return cadena

def verificar_peptido(cadena_aminoacidos):
	if len(cadena_aminoacidos) == 0:
		print("\n*********************************************************")
		print("* No ingresaste una secuencia de aminoácidos, verifica! *")
		print("*********************************************************\n")
		exit(1)

	for aminoacido in cadena_aminoacidos:
		if aminoacido in aminoacidos:
			pass
		else:
			print("\n****************************************")
			print(f'*** {aminoacido} no es un aminoácido, verifica! ***')
			print("****************************************\n\n")
			exit(1)

def pks_peptide(peptide_chain):
	"""
	Recibe una cadena de aminoacidos (peptido)
	devuelve un dicionario de pks ordenados de
	menor a mayor
	"""

	largo_cadena = len(peptide_chain)
	pk = []
	pks = {}
	ion = []
	ionizado = []
	#Funciona------------------------------------------------------------------
	for aa in range(largo_cadena):
		if aa == 0:
			amino_terminal = peptide_chain[aa]
			pk.append(aminoacidos[amino_terminal][2])
			ion.append("basic")
			ionizado.append(aminoacidos[amino_terminal][0])
			# Agregar cadena lateral del primer aa
			if len(aminoacidos[amino_terminal]) == 5:
				pk.append(aminoacidos[amino_terminal][3])
				ion.append(aminoacidos[amino_terminal][4])
				ionizado.append(aminoacidos[amino_terminal][0] + "_r")
			
			# Para los casos con un solo aminoácido
			if largo_cadena == 1:
				carboxi_terminal = peptide_chain[aa]
				pk.append(aminoacidos[carboxi_terminal][1])
				ion.append("acid")
				ionizado.append(aminoacidos[carboxi_terminal][0] + "_ct")
				break

		elif aa == largo_cadena - 1:
			carboxi_terminal = peptide_chain[aa]
			pk.append(aminoacidos[carboxi_terminal][1])
			ion.append("acid")
			ionizado.append(aminoacidos[carboxi_terminal][0])
			# Agregar cadena lateral del último aa
			if len(aminoacidos[carboxi_terminal]) == 5:
				pk.append(aminoacidos[carboxi_terminal][3])
				ion.append(aminoacidos[carboxi_terminal][4])
				ionizado.append(aminoacidos[carboxi_terminal][0] + "_r")

		elif aa > 0 and aa < largo_cadena:
			no_terminales = peptide_chain[aa]
			if len(aminoacidos[no_terminales]) == 5:
				pk.append(aminoacidos[no_terminales][3])
				ion.append(aminoacidos[no_terminales][4])
				ionizado.append(aminoacidos[no_terminales][0])

	tamaño_cuadro = len(pk)

	print("\n********** INFORMACION DE LA SECUENCIA **********\n")
	print(f"Secuencia de aminoácidos --> {peptide_chain}\n")
	c = 1
	print("Aminoácidos con cargas...\n")
	for residuos in ionizado:
		print(f"{c} --> {residuos}")
		c += 1
	#print(ion)
	#print(pk)
	#print("tamano ", tamaño_cuadro)
	print("\n*************************************************")
	#Funciona------------------------------------------------------------------

	pk_ordenado = sorted(pk)

	intervalos_pk, cantidad_intervalos = intervalos(pk_ordenado)
	pks["pks"] = intervalos_pk


	# creación de tabla que indica la ionización
	for intervalo in range(cantidad_intervalos):
		relleno = [0 for espacio in range(cantidad_intervalos)]
		if intervalo == tamaño_cuadro:
			break
		if ion[intervalo] == "basic":
			for i in range(tamaño_cuadro + 1):
				inferior, superior = intervalos_pk[i]
				#print("verificando si", pk[intervalo], " está en -->", inferior, superior)
				if pk[intervalo] == superior:
					relleno[i] = 1
				elif pk[intervalo] == inferior:
					relleno[i] = 0
				elif pk[intervalo] not in intervalos_pk[i]:
					if pk[intervalo] > superior:
						relleno[i] = 1
			pks[ionizado[intervalo]] = relleno

		elif ion[intervalo] == "acid":
			for i in range(tamaño_cuadro + 1):
				inferior, superior = intervalos_pk[i]
				#print("verificando si", pk[intervalo], " está en -->", inferior, superior)
				if pk[intervalo] == superior:
					relleno[i] = 0
				elif pk[intervalo] == inferior:
					relleno[i] = -1
				elif pk[intervalo] not in intervalos_pk[i]:
					if pk[intervalo] < superior:
						relleno[i] = -1
			pks[ionizado[intervalo]] = relleno

	# Detecta el punto del zwitterion
	for i in range(cantidad_intervalos):
		zwitterion = 0
		for j in pks:
			if j == "pks":
				pass
			else:
				zwitterion += pks[j][i]
		if zwitterion == 0:
			pk_intervalo = pks["pks"][i]
	# calcular punto isoelectrico
	#print(pks)
	punto_iso = sum(pk_intervalo) / 2
	return(round(punto_iso, 2))


def aminoacido(aminoacido):
	verificar_peptido(aminoacido)


def intervalos(pk_ordenado):
	n_intervalos = len(pk_ordenado) + 1
	intervalos = []
	for p in range(n_intervalos):
		if p == 0:
			intervalos.append([0, pk_ordenado[p]])
		elif p == n_intervalos - 1:
			intervalos.append([pk_ordenado[p-1], 14])
		else:
			intervalos.append([pk_ordenado[p-1], pk_ordenado[p]])
	return (intervalos, n_intervalos)


# Pruebas
#cadena = "ACKDG"	#5.92
#cadena = "MADRE" 	#3.95
#cadena = "KED" 	#3.95
#cadena = "AVDKQW" 	#6.89
#cadena = "DFRKTGH" #10.06
#cadena = "MA" 		#5.78
#cadena = "RCAEIY"	#6.21

#En consola
if len(argv) > 1:
	ayuda(argv[1])

else:
	presentacion()
	cadena = ingreso_cadena()
	pep = pks_peptide(cadena)
	print(f"\n--> pI: {pep}\n")

input("Presiona enter para salir")