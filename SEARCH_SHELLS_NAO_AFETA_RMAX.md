# Por Que R_max Não Muda com search_shells?

## 🎯 Entendimento Correto

**`search_shells` NÃO controla o R_max!**

O Wannier90 determina R-vectors automaticamente baseado em:

```
R_max ≈ mp_grid / 2
```

Para seu caso:
```
mp_grid = [8, 8, 8]
    ↓
R_max = [4, 4, 4]  ← Automático!
```

---

## 🔧 O Que Cada Parâmetro Faz

| Parâmetro | Função | Afeta R_max? |
|-----------|--------|-------------|
| **mp_grid** | Malha k inicial para Wannier | ✅ SIM (R_max ≈ mp_grid/2) |
| **search_shells** | Número de shells na busca por máximos localizados | ❌ NÃO |
| **hr_cutoff** | Cutoff de energia para escrever H(R) | ⚠️ Remove componentes, não adiciona |
| **num_wann** | Número de funções Wannier | ❌ NÃO |

---

## 💡 Por Que search_shells Não Funcionou?

`search_shells` é usado para:
- Determinação automática de ponto inicial para maximização
- **NÃO afeta** quantos R-vectors são gerados

É como mudar um parâmetro de otimização numérica que não tem nada a ver com a estrutura dos dados.

---

## ✅ Para Aumentar R_max: Aumente mp_grid

Você **precisa aumentar o próprio mp_grid**:

```ini
! wannier90.win

! ANTES (gera R_max ~ 4):
mp_grid = 8 8 8

! DEPOIS (gera R_max ~ 6):
mp_grid = 12 12 12

! OU (gera R_max ~ 8):
mp_grid = 16 16 16

! Escrever arquivo de Hamiltonian
write_hr = true
```

---

## ⚠️ Tradeoff

Aumentar `mp_grid` significa:

```
mp_grid = [12, 12, 12]
    ↓
12³ = 1,728 k-points (antes era 8³ = 512)
    ↓
Wannier90 leva ~3-4× mais tempo
    ↓
Gera R_max ~ 6 (antes era ~4)
```

**É computacionalmente mais caro!**

---

## 🔍 Alternativa: Verificar Estrutura do hr.dat

Talvez o Wannier90 **já esteja dando R_max maior** em certos pontos. Verifique:

```python
from TB2J.wannier import parse_ham
import numpy as np

n_wann, H_mnR, R_degens = parse_ham('wannier90_hr.dat')
Rlist = np.array(list(H_mnR.keys()))

print(f"Todos os R-points únicos:")
print(f"  Min: {np.min(Rlist, axis=0)}")
print(f"  Max: {np.max(Rlist, axis=0)}")

print(f"\nDistribuição de R:")
print(f"  R_x: {sorted(set(Rlist[:, 0]))}")
print(f"  R_y: {sorted(set(Rlist[:, 1]))}")
print(f"  R_z: {sorted(set(Rlist[:, 2]))}")

# Cound quantos R de cada "distância"
distances = np.linalg.norm(Rlist, axis=1)
for d in sorted(set(distances)):
    count = np.sum(distances == d)
    print(f"  R com norma {d}: {count} pontos")
```

---

## 🧮 O que Você Tem vs O Que Precisa

### Seu Caso Atual

```
mp_grid = [8, 8, 8]
R_max = [4, 4, 4]
Número de R-points: 9×9×9 = 729

TB2J pode usar kmesh até:
  kmesh ≤ 2×R_max ≈ [8, 8, 8]  ✅ Seguro
  kmesh = [12, 12, 12]         ❌ Aliasing!
```

### Se Aumentar para mp_grid = [12, 12, 12]

```
mp_grid = [12, 12, 12]
R_max ≈ [6, 6, 6]
Número de R-points: 13×13×13 = 2,197

TB2J pode usar kmesh até:
  kmesh ≤ 2×R_max ≈ [12, 12, 12]  ✅ Seguro
```

---

## 📋 Workflow Correto

```
PASSO 1: Seu wannier90.win ATUAL
├─ mp_grid = 8 8 8
├─ search_shells = 50 (ou qualquer valor)
└─ write_hr = true

PASSO 2: Resultado Atual
├─ R_max = 4
└─ kmesh seguro ≤ 8

PASSO 3: Para Melhorar
├─ AUMENTAR mp_grid → 12 12 12
├─ search_shells = 50 (NÃO muda isto)
└─ write_hr = true

PASSO 4: Recomputar Wannier90

PASSO 5: Novo wannier90_hr.dat
├─ R_max = 6 (aproximadamente)
└─ kmesh seguro ≤ 12
```

---

## 🔴 Aviso: Nem Sempre Funciona

Às vezes, aumentar `mp_grid` não aumenta `R_max` na proporção esperada porque:

1. **Wannier90 pode ter seu próprio cutoff interno**
   - Elementos muito pequenos de H(R) são descartados
   - Use `hr_cutoff = 1.0e-12` para ser mais agressivo

2. **A projeção pode não convergir bem em mp_grid muito grande**
   - Às vezes melhorar parâmetros de otimização ajuda

3. **Simetria do material limita R-points**
   - Materiais altamente simétricos podem ter R_max limitado

---

## ✨ Dica: Checar Primeiro Antes de Recomputar

```bash
# Se você tem arquivo hr.dat existente, simplesmente inspecione:
head -20 wannier90_hr.dat

# Procure pelos valores de R (primeiras 3 colunas de cada linha)
# São sempre números inteiros
```

---

## 📊 Resumo

| O Que Você Pensou | Realidade | Solução |
|-------------------|-----------|---------|
| `search_shells` controla R_max | Não, é só para otimização | Aumentar `mp_grid` |
| Mudei search_shells e nada muda | Correto, porque search_shells ≠ R | Aumentar `mp_grid` |
| R_max deve ir de 4 para algo maior | Sim, possível com mp_grid ↑ | mp_grid = [12,12,12] |

---

## 🎯 Sua Próxima Ação

1. **Editar wannier90.win:**
   ```ini
   mp_grid = 12 12 12    ! ← Aumentar isto
   search_shells = 50    ! ← Deixar como está
   write_hr = true
   ```

2. **Recomputar Wannier90** (vai levar mais tempo)

3. **Verificar novo R_max:**
   ```python
   from TB2J.wannier import parse_ham
   import numpy as np
   
   n_wann, H_mnR, R_degens = parse_ham('wannier90_hr.dat')
   Rlist = np.array(list(H_mnR.keys()))
   print(f"Novo R_max: {np.max(np.abs(Rlist), axis=0)}")
   ```

4. **Agora use kmesh=[12,12,12] sem aliasing! ✅**

---

## ⚡ Quick Reference

```
Para cada mp_grid, você APROXIMADAMENTE obtém:
mp_grid = [4,4,4]     → R_max ≈ 2
mp_grid = [6,6,6]     → R_max ≈ 3
mp_grid = [8,8,8]     → R_max ≈ 4      ← Seu caso atual
mp_grid = [10,10,10]  → R_max ≈ 5
mp_grid = [12,12,12]  → R_max ≈ 6
mp_grid = [14,14,14]  → R_max ≈ 7
mp_grid = [16,16,16]  → R_max ≈ 8
```

(Exato depende do algoritmo interno de Wannier90)
