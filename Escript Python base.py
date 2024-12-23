import numpy as np
import math

class MatrizZ:
    def __init__(self):
        self.atomos = []

    def adicionar_atomo(self, elemento, x, y, z):
        self.atomos.append({'elemento': elemento, 'x': x, 'y': y, 'z': z})

    def calcular_distancia(self, i, j):
        atomo_i = self.atomos[i]
        atomo_j = self.atomos[j]
        distancia = math.sqrt((atomo_j['x'] - atomo_i['x'])**2 +
                              (atomo_j['y'] - atomo_i['y'])**2 +
                              (atomo_j['z'] - atomo_i['z'])**2)
        return distancia

    def calcular_angulo(self, i, j, k):
        atomo_i = self.atomos[i]
        atomo_j = self.atomos[j]
        atomo_k = self.atomos[k]

        vec_ij = np.array([atomo_j['x'] - atomo_i['x'], atomo_j['y'] - atomo_i['y'], atomo_j['z'] - atomo_i['z']])
        vec_jk = np.array([atomo_k['x'] - atomo_j['x'], atomo_k['y'] - atomo_j['y'], atomo_k['z'] - atomo_j['z']])

        produto_escalar = np.dot(vec_ij, vec_jk)
        magnitude_ij = np.linalg.norm(vec_ij)
        magnitude_jk = np.linalg.norm(vec_jk)

        cos_theta = produto_escalar / (magnitude_ij * magnitude_jk)
        angulo = np.degrees(np.arccos(cos_theta))

        return angulo

    def calcular_diedro(self, i, j, k, l):
        atomo_i = self.atomos[i]
        atomo_j = self.atomos[j]
        atomo_k = self.atomos[k]
        atomo_l = self.atomos[l]

        vec_ij = np.array([atomo_j['x'] - atomo_i['x'], atomo_j['y'] - atomo_i['y'], atomo_j['z'] - atomo_i['z']])
        vec_jk = np.array([atomo_k['x'] - atomo_j['x'], atomo_k['y'] - atomo_j['y'], atomo_k['z'] - atomo_j['z']])
        vec_kl = np.array([atomo_l['x'] - atomo_k['x'], atomo_l['y'] - atomo_k['y'], atomo_l['z'] - atomo_k['z']])

        n1 = np.cross(vec_ij, vec_jk)
        n2 = np.cross(vec_jk, vec_kl)

        produto_escalar = np.dot(n1, n2)
        magnitude_n1 = np.linalg.norm(n1)
        magnitude_n2 = np.linalg.norm(n2)

        cos_theta = produto_escalar / (magnitude_n1 * magnitude_n2)
        angulo_diedro = np.degrees(np.arccos(cos_theta))

        return angulo_diedro

# Exemplo de uso
mol = MatrizZ()
mol.adicionar_atomo('C', 0.0, 0.0, 0.056) 
mol.adicionar_atomo('O', 1.2, 0.094, 0.067) 
mol.adicionar_atomo('O', -1.2, 1.5, 0.098)
mol.adicionar_atomo('O', -1.2, 1.0, 0.0094)

# Calcular distâncias, ângulos e ângulos diedros
distancia_mol = mol.calcular_distancia(0, 1)
angulo_mol = mol.calcular_angulo(1, 0, 2)
diedro_mol = mol.calcular_diedro(1, 0, 2, 3)

print(f"Distância: {distancia_mol:.2f} Å")
print(f"Ângulo: {angulo_mol:.2f}°")
print(f"Diedro: {diedro_mol:.2f}°")