from Bio.PDB import MMCIFParser
import numpy as np

parser = MMCIFParser(QUIET=True)
estrutura = parser.get_structure('1BNA', '1BNA.cif')

# Vamos pegar os resíduos das duas cadeias
cadeia_A = estrutura[0]['A']  # primeira modelo, cadeia A
cadeia_B = estrutura[0]['B']  # primeira modelo, cadeia B

# Vamos percorrer os resíduos da cadeia A e pegar o complementar da B
pares = {}
for residuo_A in cadeia_A.get_residues():
    num = residuo_A.get_id()[1]  # número do resíduo (ex: 1, 2, ...)
    if num in [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12]:  # são 12 pares no total
        # Pega o resíduo correspondente na cadeia B
        residuo_B = cadeia_B[num]  # mesmo número
        # Determina o par: se A for G ou C, é GC; se for A ou T, é AT
        base_A = residuo_A.get_resname()
        base_B = residuo_B.get_resname()
        # No 1BNA, a sequência da cadeia A é: C G C G A A T T C G C G
        # Portanto, os pares são:
        # 1: C-G, 2: G-C, 3: C-G, 4: G-C, 5: A-T, 6: A-T, 7: T-A, 8: T-A, 9: C-G, 10: G-C, 11: C-G, 12: G-C
        if (base_A == 'C' and base_B == 'G') or (base_A == 'G' and base_B == 'C'):
            tipo = 'GC'
        elif (base_A == 'A' and base_B == 'T') or (base_A == 'T' and base_B == 'A'):
            tipo = 'AT'
        else:
            continue  # não deveria ocorrer
        
        # Se for o primeiro GC ou primeiro AT, guardamos as coordenadas de todos os átomos de ambos
        if tipo == 'GC' and not 'GC' in pares:
            coords = []
            for atomo in list(residuo_A.get_atoms()) + list(residuo_B.get_atoms()):
                coords.append({
                    'elemento': atomo.get_element(),
                    'nome': atomo.get_name(),
                    'x': atomo.get_coord()[0],
                    'y': atomo.get_coord()[1],
                    'z': atomo.get_coord()[2]
                })
            pares['GC'] = coords
        if tipo == 'AT' and not 'AT' in pares:
            coords = []
            for atomo in list(residuo_A.get_atoms()) + list(residuo_B.get_atoms()):
                coords.append({
                    'elemento': atomo.get_element(),
                    'nome': atomo.get_name(),
                    'x': atomo.get_coord()[0],
                    'y': atomo.get_coord()[1],
                    'z': atomo.get_coord()[2]
                })
            pares['AT'] = coords

# Agora temos as coordenadas dos dois pares
print(f"Átomos no par GC: {len(pares['GC'])}")
print(f"Átomos no par AT: {len(pares['AT'])}")

# Você pode salvar em arquivos para usar no seu gerador
import json
with open('par_GC.json', 'w') as f:
    json.dump(pares['GC'], f)
with open('par_AT.json', 'w') as f:
    json.dump(pares['AT'], f)
