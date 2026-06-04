# Análise: Flutuações ao Aumentar kmesh e Interpolação de Wannier

## TL;DR - Resposta Rápida

**Seu problema é real e bem documentado:**

1. ❌ **NÃO há interpolação** do Hamiltoniano de Wannier internamente
2. ✅ **Faz sentido** usar kmesh > kmesh_wannier, **MAS com cuidado**
3. 🔴 **Flutuações grandes** indicam que você está além do limite de validade da interpolação discreta Fourier

---

## 1. Como TB2J Carrega o Hamiltoniano de Wannier

### 1.1 Arquivo `wannier90_hr.dat`

O Wannier90 escreve o Hamiltoniano em espaço real:

```
H(R) para R = (0,0,0), (1,0,0), (-1,0,0), (0,1,0), ..., (n1, n2, n3)
```

Número de R-points: **centenas a milhares** dependendo do tamanho da célula unitária.

### 1.2 Transformação para Espaço-k: **Fourier Discreta**

Em `myTB.gen_ham()`:

```python
def gen_ham(self, k, convention=2):
    Hk = np.zeros((self.nbasis, self.nbasis), dtype="complex")
    
    # Transforma H(R) → H(k) via Fourier discreta
    for iR, (R, mat) in enumerate(self.data.items()):
        phase = np.exp(self.R2kfactor * np.dot(k, R))  # e^(2πi k·R)
        H = mat * phase
        Hk += H + H.conjugate().T
    
    return Hk  # Exato para k qualquer!
```

**Ponto-chave:** Esta fórmula é **EXATA** para qualquer k ∈ [0,1]³, não apenas na malha de Wannier!

---

## 2. Por Que Aumentar kmesh Causa Flutuações?

### 2.1 A Transformada de Fourier Discreta é Periódica

A FFT só garante precisão dentro da **malha de Nyquist**:

```
kmesh_wannier = [nk1, nk2, nk3]  # p.ex. [4, 4, 4]
```

A malha de Wannier define implicitamente qual é a **malha de Nyquist**:

$$k_{\text{max}} = \frac{N_k}{2} \times \text{(passo k)}$$

Se você usar uma malha **muito maior**, passa da zona de Nyquist → **aliasing**!

### 2.2 Efeito de Aliasing

```
kmesh_wannier = [4, 4, 4]    ← produz H(R) até |R| ~ [4, 4, 4]
↓
kmesh_green = [12, 12, 12]   ← tenta avaliar H(k) em 1728 pontos
                                mas só tem 64 R-vectors!
                                
Resultado: oscilações espúrias (não-físicas)
```

### 2.3 Exemplo Numérico

Imagine H(R) com apenas componentes até R=(±2, ±2, ±2):

```
kmesh original = [4,4,4]   → precisão exata até k_Nyquist
kmesh grande   = [12,12,12] → além de k_Nyquist → ruído!
```

---

## 3. A Solução: Entender os Limites

### 3.1 Malha de Wannier vs Malha de Green

| Parâmetro | Significado | Impacto |
|-----------|-------------|--------|
| **kmesh_wannier** | Malha inicial para calcular H(R) | Define `R_max` |
| **kmesh_green** | Malha para integração de Green | Deve ser ≤ 2×kmesh_wannier |
| **R_max** | Maior R-vector no arquivo hr.dat | Limita k_Nyquist |

### 3.2 Regra Prática

```python
# Seu arquivo wannier90_hr.dat tem R-vectors até:
R_max_from_hr_file = max(abs(R) for R in Rlist)

# A malha "segura" para Green é:
kmesh_green_max = R_max_from_hr_file

# Exemplo:
# Se hr.dat tem R até (4,4,4):
kmesh_wannier = [4, 4, 4]      # ✅ OK
kmesh_green   = [8, 8, 8]      # ⚠️ Limite
kmesh_green   = [12, 12, 12]   # ❌ ALÉM DO LIMITE
```

---

## 4. Convergência "Não-Física" com kmesh Grande

### 4.1 Tipos de Flutuações que Você Vê

```
kmesh = [3, 3, 3]  → resultado X
kmesh = [4, 4, 4]  → resultado X (converge)
kmesh = [6, 6, 6]  → resultado Y ← oscilação!
kmesh = [8, 8, 8]  → resultado Z ← continua oscilando
kmesh = [12,12,12] → resultado W ← não converge?
```

**Causa:** Não é convergência Fourier, é **aliasing** da FFT discreta.

### 4.2 Como Detectar Isso

Você está vendo aliasing se:

1. Resultados oscilam em vez de convergir suavemente
2. A energia varia com kmesh de forma não-monónica
3. Aumentar kmesh piora em vez de melhorar

---

## 5. A Resposta Correta à Sua Pergunta

### 5.1 "Faz Sentido Usar kmesh_green > kmesh_wannier?"

**Resposta: SIM, mas com limite.**

```
A relação ideal é:  kmesh_green ≤ 2 × kmesh_wannier

Exemplos válidos:
✅ kmesh_wannier = [4,4,4]  → kmesh_green = [4,4,4], [6,6,6], [8,8,8]
✅ kmesh_wannier = [6,6,6]  → kmesh_green = [6,6,6], [8,8,8], [12,12,12]
❌ kmesh_wannier = [4,4,4]  → kmesh_green = [16,16,16]  ← aliasing!
```

**Por quê?** Porque os R-vectors no hr.dat só vão até ~2×kmesh_wannier.

### 5.2 "Ele Interpola o Hamiltoniano Internamente?"

**Resposta: NÃO.**

O TB2J **apenas** faz Fourier discreta:

```python
H(k) = Σ_R H(R) · e^(2πi k·R)
```

**Sem nenhuma interpolação suave.** Isto é exato mas limitado a k-points dentro da malha de Nyquist.

---

## 6. Como Investigar Seu Caso Específico

### 6.1 Inspecionar o arquivo hr.dat

```bash
# Conte quantos R-vectors tem no arquivo
wc -l wannier90_hr.dat

# Veja os R-vectors máximos
head -50 wannier90_hr.dat
```

### 6.2 Verificar Convergência Adequada

```python
from TB2J.wannier import parse_ham

n_wann, H_mnR, R_degens = parse_ham('path/to/wannier90_hr.dat')

# Encontre R_max
Rlist = list(H_mnR.keys())
R_max = np.max(np.abs(Rlist), axis=0)

print(f"R_max: {R_max}")
print(f"Máximo valor: {max(R_max)}")

# Kmesh seguro seria ~ 2×max(R_max)
kmesh_safe = 2 * max(R_max)
print(f"kmesh seguro ≤ {kmesh_safe}")
```

### 6.3 Teste de Convergência Correto

```python
from TB2J.green import TBGreen

kmeshes = [[3,3,3], [4,4,4], [6,6,6], [8,8,8]]

for kmesh in kmeshes:
    green = TBGreen(tbmodel, kmesh=kmesh)
    # calcule suas grandezas de interesse
    J = compute_exchange(green)
    print(f"kmesh {kmesh}: J = {J}")
```

**O gráfico deve ser monótono até o limite, depois divergir ou oscilar.**

---

## 7. Soluções Práticas

### 7.1 Opção A: Aumentar R-vectors no Wannier90

```ini
# wannier90.win
# Aumentar alcance:
search_shells = 100    # ← aumentar isto
num_print_cycles = 100
```

Isto faz Wannier90 incluir R-vectors até maiores distâncias → permite kmesh maior.

### 7.2 Opção B: Usar IBZ

```python
green = TBGreen(
    tbmodel,
    kmesh=[8, 8, 8],
    ibz=True,           # ← reduz k-points efetivos mantendo precisão
)
```

IBZ usa simetria para reduzir de ~512 para ~64 k-points enquanto mantém convergência.

### 7.3 Opção C: Usar kmesh Conservador

```python
# Regra de ouro:
kmesh_green = kmesh_wannier  # Sem riscos!

# Ou:
kmesh_green = 1.5 * kmesh_wannier  # Moderadamente agressivo
```

---

## 8. Flutuações Muito Grandes: Outro Problema?

Se as flutuações são **extremamente grandes** (não apenas oscilações pequenas), pode indicar:

### 8.1 Ruído no Hamiltoniano de Wannier

```python
# Verifique os elementos de H(R)
H0 = H_mnR[(0, 0, 0)]
print(f"H(0,0,0) diagonal: {np.diag(H0)}")

# Se houver valores muito pequenos (< 1e-6 eV):
#   → ruído numérico no cálculo do Wannier
#   → afeta principalmente kmesh alto
```

### 8.2 Problema de Gauge de Wannier

```python
# Se H(R) não é Hermitiana por simetria:
for R in Rlist:
    H_minus_R = np.conj(H_mnR[tuple(-np.array(R))]).T
    if not np.allclose(H_mnR[R], H_minus_R):
        print(f"⚠️ H({R}) não-simétrica!")
```

---

## 9. Resumo Executivo

| Pergunta | Resposta | Por quê? |
|----------|----------|---------|
| **Interpola H(k)?** | ❌ Não | Usa Fourier discreta exata |
| **kmesh_green > kmesh_wannier?** | ✅ Sim, até 2× | FFT tem limite de Nyquist |
| **Flutuações são normais?** | ⚠️ Pequenas sim | Aliasing além do limite |
| **Como convergir?** | kmesh ≤ 2×R_max | Respeita limite Fourier |

---

## 10. Exemplo Completo: Debug do Seu Caso

```python
import numpy as np
from TB2J.wannier import parse_ham
from TB2J.green import TBGreen

# PASSO 1: Inspecione o arquivo
n_wann, H_mnR, R_degens = parse_ham('wannier90_hr.dat')
Rlist = np.array(list(H_mnR.keys()))
R_max = np.max(np.abs(Rlist), axis=0)
print(f"R_max encontrado: {R_max}")
print(f"Número de R-points: {len(H_mnR)}")

# PASSO 2: Teste convergência
results = {}
for kmesh_val in [3, 4, 5, 6, 8]:
    kmesh = [kmesh_val, kmesh_val, kmesh_val]
    try:
        green = TBGreen(tbmodel, kmesh=kmesh, efermi=your_efermi)
        # Calcule sua grandeza
        J = compute_exchange(green)
        results[kmesh_val] = J
        print(f"kmesh=[{kmesh_val}³]: J = {J:.6f} eV")
    except Exception as e:
        print(f"kmesh=[{kmesh_val}³]: ERRO - {e}")

# PASSO 3: Analise o padrão
print("\nAnálise:")
diffs = [abs(results[list(results.keys())[i]] - results[list(results.keys())[i-1]]) 
         for i in range(1, len(results))]
if all(d < 0.001 for d in diffs):
    print("✅ Convergência suave até kmesh máximo")
else:
    print(f"⚠️ Oscilações detectadas: diffs = {diffs}")
    print(f"   Limite seguro: kmesh ≤ {2*max(R_max)}")
```

---

## Referências Teóricas

- **FFT Nyquist Limit**: Press et al., Numerical Recipes, Cap. 12
- **Wannier Interpolation**: Marzari & Vanderbilt, PRL 82, 3296 (1999)
- **TB2J Green's Functions**: Kovalev et al., PRB 91, 224405 (2015)
