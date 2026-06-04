# Análise: Integração das Funções de Green em green.py

## 1. O que são Funções de Green (Green's Functions)

As funções de Green são ferramentas matemáticas para resolver equações diferenciais. No contexto de TB2J (Tight Binding to Magnons and Exchange parameters):

$$G(k, \epsilon) = \frac{1}{\epsilon + E_F - H(k)}$$

Onde:
- $\epsilon$: energia relativa ao nível de Fermi
- $E_F$: energia de Fermi
- $H(k)$: Hamiltoniano em espaço-k

No espaço real:
$$G(R, \epsilon) = \sum_k G(k, \epsilon) \cdot e^{-2\pi i \mathbf{R} \cdot \mathbf{k}} \cdot w_k$$

---

## 2. Fluxo de Cálculo das Funções de Green

```
┌─────────────────────────────────────┐
│ TBGreen.__init__()                   │
├─────────────────────────────────────┤
│ 1. Prepara k-points (prepare_kpts)   │
│    ├─ Monkhorst-Pack grid            │
│    ├─ IBZ (Irreducible Brillouin Z.) │
│    └─ k-points customizados          │
└────────┬────────────────────────────┘
         │
         ▼
┌─────────────────────────────────────┐
│ _prepare_eigen()                    │
├─────────────────────────────────────┤
│ 2. Calcula H, S, evals, evecs       │
│    em todos os k-points             │
│                                      │
│ 3. Calcula energia de Fermi          │
│    se não fornecida                 │
│                                      │
│ 4. Filtra bandas por energia         │
│    emin < E < emax                  │
└────────┬────────────────────────────┘
         │
         ▼
┌─────────────────────────────────────┐
│ get_Gk(ik, energy)                  │
├─────────────────────────────────────┤
│ 5. Calcula G(k, ε) via:             │
│    eigen_to_G(evals, evecs, ε)      │
│                                      │
│    G = Σ_i v_i v_i†/(E-ε_i)         │
└────────┬────────────────────────────┘
         │
         ▼
┌─────────────────────────────────────┐
│ get_Gk_all(energy)                  │
├─────────────────────────────────────┤
│ 6. Calcula G(k,ε) para todos os     │
│    k-points [nkpts × nbasis × nbasis]│
└────────┬────────────────────────────┘
         │
         ▼
┌─────────────────────────────────────┐
│ compute_GR(Rpts, kpts, Gks)         │
├─────────────────────────────────────┤
│ 7. Integração Fourier inversa:      │
│    G(R) = Σ_k G(k)·e^(-2πiR·k)·w_k │
│                                      │
│    Usa einsum otimizado             │
└────────┬────────────────────────────┘
         │
         ▼
┌─────────────────────────────────────┐
│ G(R, ε) [nR × nbasis × nbasis]      │
│ Funções de Green no espaço real      │
└─────────────────────────────────────┘
```

---

## 3. Cálculo Detalhado de G(k, ε)

### 3.1 Fórmula de Autovalores/Autovetores

```python
def eigen_to_G(evals, evecs, efermi, energy):
    """
    Calcula: G = Σ_i |v_i⟩ 1/(E - ε_i) ⟨v_i|
    
    Parâmetros:
    - evals: autovalores [nbasis]
    - evecs: autovetores [nbasis × nbasis]
    - efermi: energia de Fermi
    - energy: E - Ef (energia relativa)
    """
    return np.einsum(
        "ib, b, jb-> ij",  # |v_i⟩ · 1/(E-ε_i) · ⟨v_i|
        evecs,              # autovetores |v_i⟩
        1.0 / (-evals + (energy + efermi)),  # fatores 1/(E-ε_i)
        evecs.conj(),       # conjugado dos autovetores
        optimize=True,
    )
```

**Operação Einstein notation**: 
- `evecs[i,b]` × `1/(-evals[b] + (energy + efermi))` × `evecs.conj()[b,j]`
- Resultado: matriz `[i,j]` de tamanho `[nbasis, nbasis]`

### 3.2 Fórmula Alternativa (Direta)

```python
def eigen_to_G2(H, S, efermi, energy):
    """
    Calcula: G = ((E + Ef)·S - H)^(-1)
    
    Parâmetros:
    - H: Hamiltoniano na base dos eigenvetores
    - S: matriz de overlap (ortogonalidade)
    - efermi: energia de Fermi
    - energy: E - Ef
    """
    return np.linalg.inv((energy + efermi) * S - H)
```

---

## 4. Integração para Espaço Real: `compute_GR()`

```python
def compute_GR(self, Rpts, kpts, Gks):
    """
    Transforma G(k) → G(R) via transformada de Fourier inversa
    
    G(R,ε) = (1/Nk) Σ_k G(k,ε) · e^(-2πi R·k)
    """
    Rvecs = np.array(Rpts)  # [nR, 3]
    
    # Calcula matriz de fases: e^(-2πi R·k)
    # shape: [nR, nk]
    phase = np.exp(self.k2Rfactor * np.einsum("ni,mi->nm", Rvecs, kpts))
    
    # Multiplica por pesos dos k-points
    # k2Rfactor = -R2kfactor = -2π (convencional)
    phase *= self.kweights[None]  # broadcasting: [nR, nk]
    
    # Integração final usando einsum otimizado
    # "kij,rk->rij" significa:
    # GR[r,i,j] = Σ_k Gks[k,i,j] * phase[r,k]
    GR = np.einsum("kij,rk->rij", Gks, phase, optimize="optimal")
    
    return GR  # shape: [nR, nbasis, nbasis]
```

---

## 5. Parâmetros que Influenciam a Integração

### 5.1 Parâmetros de k-points (maior impacto)

| Parâmetro | Tipo | Padrão | Efeito |
|-----------|------|--------|--------|
| **kmesh** | tuple | None | Grade de Monkhorst-Pack [nk1, nk2, nk3]. Maior = mais k-points = melhor resolução |
| **ibz** | bool | False | Se True: usa zona de Brillouin irreduzível (reduz nkpts) |
| **gamma** | bool | False | Se True: inclui ponto Gamma (0,0,0) na grade |
| **kpts** | array | None | K-points customizados (sobrescreve kmesh) |
| **kweights** | array | None | Pesos para k-points customizados |

**Impacto na integração:**
```
Nkpts ↑ → integração mais precisa → maior custo computacional
Nkpts ↓ → integração menos precisa → mais rápido
```

**Exemplo:**
```python
# Caso 1: Baixa resolução
TBGreen(tbmodel, kmesh=[3,3,3])  # 27 k-points

# Caso 2: Alta resolução  
TBGreen(tbmodel, kmesh=[12,12,12])  # 1728 k-points ⚠️ caro!

# Caso 3: IBZ reduz k-points efetivamente
TBGreen(tbmodel, kmesh=[6,6,6], ibz=True)  # ~46 k-points (com simetria)
```

### 5.2 Parâmetros de Energia

| Parâmetro | Tipo | Padrão | Efeito |
|-----------|------|--------|--------|
| **efermi** | float | Auto | Energia de Fermi. Se None, calculada automaticamente |
| **initial_emin** | float | -25 | Energia mínima relativa a EF para filtrar bandas |
| **energy** | float | - | Energia em que G(ε) é calculada no método `get_Gk()` |

**Impacto na integração:**
```python
# Em _prepare_eigen():
emax = self.efermi + 5.1  # Fixo (linha 269!)
emin = self.efermi + self.adjusted_emin  # Automático baseado em gaps

# Reduz autovetores/autovalores para bandas no intervalo [emin, emax]
self.evals, self.evecs = self._reduce_eigens(
    self.evals,
    self.evecs,
    emin=emin,
    emax=emax,
)
```

**Número de bandas (nbasis) ↓ → G calcula mais rápido**

### 5.3 Parâmetros de Integração Fourier

| Parâmetro | Valor | Significado | Impacto |
|-----------|-------|-------------|---------|
| **k2Rfactor** | -2π | Fator de fase = -R2kfactor | Define escala da transformada |
| **kweights** | 1/Nk | Pesos normalizados | Importante para IBZ com pesos ≠ 1/Nk |

**A normalização:**
```python
if ibz:
    # Pesos assimétricos (maior nos pontos mais simétricos)
    self.kweights = weights_from_IBZ  # variam por simetria
else:
    # Pesos uniformes
    self.kweights = np.array([1.0 / self.nkpts] * self.nkpts)
```

### 5.4 Parâmetros de Otimização

| Parâmetro | Tipo | Padrão | Efeito |
|-----------|------|--------|--------|
| **nproc** | int | 1 | Número de processos para calcular H,S,evals,evecs |
| **use_cache** | bool | False | Se True: usa memmap para dados grandes |
| **cache_path** | str | /dev/shm | Diretório para cache (RAM virtual) |

---

## 6. Integração com Derivadas

### 6.1 `get_GR_and_dGRdx()`

```python
def get_GR_and_dGRdx(self, Rpts, energy, dHdx):
    """
    Calcula G(R) e dG(R)/dx simultaneamente
    
    dG(k)/dx = G(k) · (dH(k)/dx) · G(k)   [Regra da cadeia]
    dG(R)/dx = Σ_k dG(k)/dx · e^(-2πiR·k) · w_k
    """
    for ik, kpt in enumerate(self.kpts):
        Gk = self.get_Gk(ik, energy)      # G(k,ε)
        Gkp = Gk * self.kweights[ik]      # G(k,ε) · w_k
        dHk = dHdx.gen_ham(tuple(kpt))    # dH/dx em espaço-k
        
        # Derivada da função de Green
        dG = Gk @ dHk @ Gkp               # G · (dH/dx) · G
        
        for iR, R in enumerate(Rpts):
            phase = np.exp(self.k2Rfactor * np.dot(R, kpt))
            
            # Acumula G(R) e dG(R)/dx
            GR[R] += Gkp * (phase * self.kweights[ik])
            dGRdx[R] += dG * (phase * self.kweights[ik])
    
    return GR, dGRdx
```

---

## 7. Exemplo Completo de Uso

```python
from TB2J.green import TBGreen

# Criar objeto Green com diferentes configurações
green = TBGreen(
    tbmodel=my_tb_model,
    
    # k-points
    kmesh=[6, 6, 6],      # Grade 6×6×6 = 216 pontos
    ibz=False,            # Usar zona de Brillouin completa
    gamma=True,           # Incluir ponto Gamma
    
    # Energia
    efermi=None,          # Calcular automaticamente
    initial_emin=-25,     # Filtrar bandas de -25 até 5.1 eV
    
    # Paralelização
    nproc=4,              # Usar 4 processos
    use_cache=False,      # Não usar cache (memmap)
)

# Calcular G(k, ε) para energia ε = 0.5 eV relativa ao EF
energy = 0.5
Gk_all = green.get_Gk_all(energy)  # shape: [216, nbasis, nbasis]

# Integrar para espaço real em pontos R específicos
Rpts = [[0,0,0], [1,0,0], [1,1,0], [1,1,1]]
GR = green.get_GR(Rpts, energy)  # shape: [4, nbasis, nbasis]

# Com derivadas (para cálculo de parâmetros de exchange)
from TB2J.exchange_params import DerivativeOperator
dHdx = DerivativeOperator(...)
GR, dGRdx = green.get_GR_and_dGRdx(Rpts, energy, dHdx)
```

---

## 8. Sumário dos Parâmetros Críticos

### Afetam **Precisão da Integração:**
1. **kmesh** / **ibz**: Determina número de k-points
2. **kweights**: Pesos na integração Fourier
3. **efermi**: Centro da integração em energia
4. **initial_emin / emax**: Intervalo de bandas

### Afetam **Velocidade do Cálculo:**
1. **nproc**: Paralelização de cálculo de H,S,evals
2. **kmesh** (número de k-points)
3. **use_cache**: Trade-off memória ↔ I/O
4. **initial_emin** (filtra bandas desnecessárias)

### Afetam **Convergência:**
- **kmesh** (convergência k-point)
- **efermi** (precisa estar correto)
- **energy** em `get_Gk()` (onde G é avaliada)

---

## 9. Linha 269: O Valor Mágico 5.1

```python
emax=self.efermi + 5.1,  # Linha 269
```

**Por que 5.1?**
- Valor fixo que define o limite superior de energia
- Provavelmente baseado em observações de que bandas acima de 5.1 eV acima do EF
  não contribuem significativamente para propriedades magnéticas
- Reduz tamanho de `evecs` e `evals` mantendo precisão
- Específico para o problema TB2J (magnetismo)

**Pode ser ajustado se:**
- Estudar níveis de energia muito altos
- Trabalhar com bandas de condução profundas
