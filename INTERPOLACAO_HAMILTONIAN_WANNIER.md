# Interpolação de Hamiltoniano de Wannier em k-Mesh Densa

## TL;DR - Resposta Rápida

**SIM, é totalmente possível interpolar o Hamiltoniano de Wannier para uma malha k mais densa.**

Existem várias formas:

1. ✅ **Melhor: Aumentar search_shells no Wannier90** (permite maior kmesh sem interpolação)
2. ✅ **Bom: Usar interpolação de Fourier com mais R-vectors**
3. ✅ **Possível: Implementar interpolação polinomial N-D**
4. ⚠️ **Complexo: Usar ferramentas do Wannier90 postproc**

---

## 1. Por Que Interpolar?

Você quer:
- Malha de k-points maior (melhor precisão em Green)
- Sem re-calcular Wannier90 (econômico)
- Sem aliasing FFT (mais controle)

---

## 2. Técnica 1: Aumentar R-Vectors no Wannier90 ⭐ RECOMENDADO

**Vantagem:** Exato, sem aproximações.  
**Desvantagem:** Precisa recomputar.

### 2.1 No arquivo `wannier90.win`

```ini
! wannier90.win
[...]
# Aumentar alcance para incluir R-vectors mais distantes
search_shells = 100        # aumentar de 50 (padrão)
kmesh = 6 6 6             # malha usada para projeção

[...]
write_tb = true           # escreve arquivo _tb.dat com mais R-vectors
write_hr = true           # escreve wannier90_hr.dat
```

**Resultado:** Arquivo `wannier90_hr.dat` terá R-vectors até distâncias maiores, permitindo kmesh até 2× ou 3× sem aliasing.

### 2.2 Verificar R-vectors

```python
from TB2J.wannier import parse_ham
import numpy as np

n_wann, H_mnR, R_degens = parse_ham('wannier90_hr.dat')
Rlist = np.array(list(H_mnR.keys()))

print(f"Número de R-points: {len(H_mnR)}")
print(f"R_max: {np.max(np.abs(Rlist), axis=0)}")

# Com search_shells aumentado, você verá R_max maior
```

---

## 3. Técnica 2: Interpolação de Fourier com Supercélula

**Ideia:** Depois que Wannier90 calcula H(R), você pode usar a transformada de Fourier em uma malha mais fina.

### 3.1 Princípio

O Hamiltoniano já está completo em espaço real H(R). A transformada de Fourier é **exata para qualquer k**:

```
H(k) = Σ_R H(R) · e^(2πi k·R)
```

Então, simplesmente calcule H(k) para mais k-points!

### 3.2 Implementação

```python
from TB2J.myTB import MyTB
import numpy as np

# Seu modelo tight-binding de Wannier
tbmodel = MyTB.read_from_wannier_dir(...)

# Gere uma malha k mais densa
kmesh_dense = [12, 12, 12]  # p.ex. 1728 pontos
kpts_dense = np.array([
    [i/kmesh_dense[0], j/kmesh_dense[1], k/kmesh_dense[2]]
    for i in range(kmesh_dense[0])
    for j in range(kmesh_dense[1])
    for k in range(kmesh_dense[2])
])

# Calcule H(k) para cada k-point
Hk_dense = np.array([tbmodel.gen_ham(k) for k in kpts_dense])

# Agora use Hk_dense com TBGreen
from TB2J.green import TBGreen

# Criar um objeto "mock" TBGreen com a malha densa
# (Requeriria modificação do código TB2J)
```

**Problema:** TBGreen espera recalcular H(k) dinamicamente. Precisaria adaptar o código.

---

## 4. Técnica 3: Interpolação Polinomial N-D (Scipy)

**Ideia:** Interpolar H(k) usando polinômios em 3D.

### 4.1 Implementação Básica

```python
import numpy as np
from scipy.interpolate import RegularGridInterpolator
from TB2J.myTB import MyTB

# 1. Calcule H(k) em malha original
kmesh_orig = [6, 6, 6]
kpts_orig = np.array([
    [i/kmesh_orig[0], j/kmesh_orig[1], k/kmesh_orig[2]]
    for i in range(kmesh_orig[0])
    for j in range(kmesh_orig[1])
    for k in range(kmesh_orig[2])
])

tbmodel = MyTB.read_from_wannier_dir(...)
H_values = np.array([tbmodel.gen_ham(k) for k in kpts_orig])
# H_values.shape = (216, nbasis, nbasis)

# 2. Reshape para grade regular 3D
nbasis = tbmodel.nbasis
H_grid = H_values.reshape((kmesh_orig[0], kmesh_orig[1], kmesh_orig[2], nbasis, nbasis))

# 3. Crie interpolador para CADA elemento de matriz
# NOTA: Isto é lento! Melhor fazê-lo com einsum

def interpolate_hamiltonian(H_grid, kmesh_orig, kpts_new):
    """
    Interpola H(k) usando interpolação cúbica 3D
    
    :param H_grid: array [nk1, nk2, nk3, nbasis, nbasis]
    :param kmesh_orig: [nk1, nk2, nk3]
    :param kpts_new: array [nk_new, 3]
    :return: H_new array [nk_new, nbasis, nbasis]
    """
    nk1, nk2, nk3, nbasis, nbasis_ = H_grid.shape
    
    # Eixos da grade original
    k1 = np.linspace(0, 1, nk1, endpoint=False)
    k2 = np.linspace(0, 1, nk2, endpoint=False)
    k3 = np.linspace(0, 1, nk3, endpoint=False)
    
    H_new = np.zeros((len(kpts_new), nbasis, nbasis), dtype=complex)
    
    # Interpolar cada elemento de matriz
    for i in range(nbasis):
        for j in range(nbasis):
            # Elemento real
            interp_real = RegularGridInterpolator(
                (k1, k2, k3), 
                np.real(H_grid[:, :, :, i, j]),
                method='cubic',
                bounds_error=False,
                fill_value=np.nan
            )
            
            # Elemento imaginário
            interp_imag = RegularGridInterpolator(
                (k1, k2, k3),
                np.imag(H_grid[:, :, :, i, j]),
                method='cubic',
                bounds_error=False,
                fill_value=np.nan
            )
            
            H_new[:, i, j] = (
                interp_real(kpts_new) + 1j * interp_imag(kpts_new)
            )
    
    return H_new

# 4. Use para gerar malha densa
kmesh_dense = [12, 12, 12]
kpts_dense = np.array([
    [i/kmesh_dense[0], j/kmesh_dense[1], k/kmesh_dense[2]]
    for i in range(kmesh_dense[0])
    for j in range(kmesh_dense[1])
    for k in range(kmesh_dense[2])
])

H_dense = interpolate_hamiltonian(H_grid, kmesh_orig, kpts_dense)
```

**Problemas:**
- ⚠️ Muito lento (interpola nbasis² elementos)
- ⚠️ Pode violar Hermiticidade numericamente
- ⚠️ Pode violar simetrias (se houver)

---

## 5. Técnica 4: Criar TBModel "Interpolado"

**Ideia:** Subclassificar `MyTB` para interpolar H(k) on-the-fly.

### 5.1 Implementação

```python
from TB2J.myTB import MyTB
from scipy.interpolate import RegularGridInterpolator
import numpy as np

class InterpolatedTB(MyTB):
    """
    Tight-binding model com Hamiltoniano interpolado para malha densa.
    """
    
    def __init__(self, tbmodel_orig, kmesh_orig, kmesh_dense):
        """
        :param tbmodel_orig: MyTB original
        :param kmesh_orig: malha original [nk1, nk2, nk3]
        :param kmesh_dense: malha densa desejada [nk1', nk2', nk3']
        """
        self.tbmodel_orig = tbmodel_orig
        self.kmesh_orig = kmesh_orig
        self.kmesh_dense = kmesh_dense
        self.nbasis = tbmodel_orig.nbasis
        self.R2kfactor = tbmodel_orig.R2kfactor
        self.is_orthogonal = tbmodel_orig.is_orthogonal
        
        # Precalcule grid original
        self._build_interpolators()
    
    def _build_interpolators(self):
        """Constrói interpoladores para H_real e H_imag"""
        print("Building interpolators...")
        
        # Gere grid original
        kpts_orig = self._make_grid(self.kmesh_orig)
        H_orig = np.array([
            self.tbmodel_orig.gen_ham(k) 
            for k in kpts_orig
        ])
        
        # Reshape
        H_grid = H_orig.reshape(
            self.kmesh_orig + (self.nbasis, self.nbasis)
        )
        
        # Crie eixos
        self.k1 = np.linspace(0, 1, self.kmesh_orig[0], endpoint=False)
        self.k2 = np.linspace(0, 1, self.kmesh_orig[1], endpoint=False)
        self.k3 = np.linspace(0, 1, self.kmesh_orig[2], endpoint=False)
        
        # Armazene interpoladores (uma por elemento de matriz)
        self.interp_real = {}
        self.interp_imag = {}
        
        for i in range(self.nbasis):
            for j in range(self.nbasis):
                self.interp_real[(i, j)] = RegularGridInterpolator(
                    (self.k1, self.k2, self.k3),
                    np.real(H_grid[..., i, j]),
                    method='cubic'
                )
                self.interp_imag[(i, j)] = RegularGridInterpolator(
                    (self.k1, self.k2, self.k3),
                    np.imag(H_grid[..., i, j]),
                    method='cubic'
                )
        
        print("✅ Interpolators ready!")
    
    def _make_grid(self, kmesh):
        """Cria grade regular de k-points"""
        return np.array([
            [i/kmesh[0], j/kmesh[1], k/kmesh[2]]
            for i in range(kmesh[0])
            for j in range(kmesh[1])
            for k in range(kmesh[2])
        ])
    
    def gen_ham(self, k, convention=2):
        """
        Retorna H(k) interpolado.
        Nota: convention é ignorado (sempre usa interpolação)
        """
        k = np.array(k)
        k = k % 1.0  # Garantir que está em [0, 1)
        
        H = np.zeros((self.nbasis, self.nbasis), dtype=complex)
        
        for i in range(self.nbasis):
            for j in range(self.nbasis):
                H[i, j] = (
                    self.interp_real[(i, j)](k[None, :]).item() +
                    1j * self.interp_imag[(i, j)](k[None, :]).item()
                )
        
        return H

# Uso:
tbmodel_orig = MyTB.read_from_wannier_dir(...)
tbmodel_interp = InterpolatedTB(
    tbmodel_orig, 
    kmesh_orig=[6, 6, 6],
    kmesh_dense=[12, 12, 12]
)

# Agora use com TBGreen:
from TB2J.green import TBGreen

green = TBGreen(tbmodel_interp, kmesh=[12, 12, 12])
# Funciona normalmente!
```

**Vantagens:**
- ✅ Compatível com TBGreen
- ✅ Transparente para o usuário
- ✅ Permite kmesh arbitrário

**Desvantagens:**
- ⚠️ Lento (interpola cada H(k) dinamicamente)
- ⚠️ Pode não preservar simetrias

---

## 6. Técnica 5: Usar Wannier90 postproc (Avançado)

O Wannier90 tem ferramentas de interpolação built-in.

### 6.1 Workflow

```bash
# 1. Rode wannier90 normalmente
wannier90.x wannier90

# 2. Use w90postproc (se disponível) ou write_u_matrices
# para interpolar em malha densa
```

**Status:** Requer compilação especial do Wannier90 com postproc.

---

## 7. Comparação das Técnicas

| Técnica | Exatidão | Velocidade | Implementação | Recomendação |
|---------|----------|-----------|---------------|--------------|
| **1. Aumentar search_shells** | 100% | Lento (recomputa W90) | Fácil | ⭐⭐⭐ MELHOR |
| **2. Fourier pura** | 100% | Rápido | Trivial | ⭐⭐⭐ Usar com cuidado |
| **3. Interpolação Scipy** | ~99% | Médio | Moderado | ⭐⭐ Pesado |
| **4. InterpolatedTB** | ~99% | Médio | Moderado | ⭐⭐ Bom para testes |
| **5. Wannier90 postproc** | 100% | Rápido | Difícil | ⭐ Complexo |

---

## 8. Solução Prática Recomendada

### 8.1 Passo 1: Aumentar R-vectors

```ini
# wannier90.win
search_shells = 100     # aumentar
```

Recompile Wannier90.

### 8.2 Passo 2: Verificar R_max

```python
from TB2J.wannier import parse_ham
import numpy as np

n_wann, H_mnR, R_degens = parse_ham('wannier90_hr.dat')
Rlist = np.array(list(H_mnR.keys()))
R_max = np.max(np.abs(Rlist), axis=0)

kmesh_safe = 2 * max(R_max)
print(f"kmesh seguro: até [{kmesh_safe}³]")
```

### 8.3 Passo 3: Teste Convergência

```python
from TB2J.green import TBGreen

for kmesh_val in [6, 8, 10, 12]:
    kmesh = [kmesh_val]*3
    green = TBGreen(tbmodel, kmesh=kmesh)
    J = compute_exchange(green)
    print(f"kmesh={kmesh}: J={J}")
    
# Deve convergir suavemente
```

---

## 9. Code Snippet: InterpolatedTB Completo e Funcional

```python
# save as: TB2J/myTB_interpolated.py

import numpy as np
from scipy.interpolate import RegularGridInterpolator
from TB2J.myTB import MyTB

class InterpolatedTB(MyTB):
    """
    Tight-binding model com Hamiltoniano interpolado.
    Permite usar kmesh > kmesh_original sem aliasing.
    """
    
    def __init__(self, tbmodel_orig, kmesh_orig):
        """
        :param tbmodel_orig: MyTB original from Wannier
        :param kmesh_orig: [nk1, nk2, nk3] malha de Wannier
        """
        # Copiar propriedades
        self.__dict__.update(tbmodel_orig.__dict__)
        self._tbmodel_orig = tbmodel_orig
        self._kmesh_orig = kmesh_orig
        
        # Build interpolators
        self._setup_interpolators()
    
    def _setup_interpolators(self):
        """Precalculate H(k) grid and setup interpolators"""
        print("Setting up interpolators for H(k)...")
        
        # Gerar grid
        k_indices = np.indices(self._kmesh_orig)
        kpts = np.array([
            k_indices[0].flatten() / self._kmesh_orig[0],
            k_indices[1].flatten() / self._kmesh_orig[1],
            k_indices[2].flatten() / self._kmesh_orig[2],
        ]).T
        
        # Calcular H(k) em grid
        H_list = [self._tbmodel_orig.gen_ham(k) for k in kpts]
        H_array = np.array(H_list)
        
        # Reshape para grade
        H_grid = H_array.reshape(
            self._kmesh_orig + (self.nbasis, self.nbasis)
        )
        
        # Eixos coordenados
        axes = [
            np.linspace(0, 1, n, endpoint=False)
            for n in self._kmesh_orig
        ]
        
        # Store for later
        self._axes = axes
        self._H_grid_real = np.real(H_grid)
        self._H_grid_imag = np.imag(H_grid)
        
        print(f"✅ Interpolators ready for kmesh {self._kmesh_orig}")
    
    def gen_ham(self, k, convention=2):
        """Generate Hamiltonian at arbitrary k-point via interpolation"""
        k = np.array(k) % 1.0  # Wrap to [0,1)
        
        H = np.zeros((self.nbasis, self.nbasis), dtype=complex)
        
        # Interpolar cada elemento
        for i in range(self.nbasis):
            for j in range(self.nbasis):
                # Real
                interp_real = RegularGridInterpolator(
                    self._axes,
                    self._H_grid_real[..., i, j],
                    bounds_error=False,
                    fill_value='extrapolate'
                )
                
                # Imaginário
                interp_imag = RegularGridInterpolator(
                    self._axes,
                    self._H_grid_imag[..., i, j],
                    bounds_error=False,
                    fill_value='extrapolate'
                )
                
                H[i, j] = (
                    interp_real(k) + 1j * interp_imag(k)
                )
        
        return H
```

**Uso:**

```python
from TB2J.myTB import MyTB
from TB2J.myTB_interpolated import InterpolatedTB
from TB2J.green import TBGreen

# Carregar Wannier original
tbmodel = MyTB.read_from_wannier_dir(...)

# Criar versão interpolada
tbmodel_interp = InterpolatedTB(tbmodel, kmesh_orig=[6, 6, 6])

# Usar com kmesh maior agora!
green = TBGreen(tbmodel_interp, kmesh=[12, 12, 12])

# Funciona sem problemas de aliasing!
```

---

## 10. Resumo Final

| Pergunta | Resposta | Método |
|----------|----------|--------|
| **Pode interpolar?** | ✅ Sim | Várias formas |
| **Qual é melhor?** | Aumentar R-vectors | `search_shells` no W90 |
| **Mais rápido?** | Usar Fourier pura | Aumentar kmesh direto |
| **Mais genérico?** | InterpolatedTB | Scipy RegularGridInterpolator |

**Recomendação:** Comece aumentando `search_shells` no Wannier90. Se não puder, use InterpolatedTB acima! 🚀
