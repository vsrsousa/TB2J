# Guia Consolidado: mp_grid, R_max e kmesh em TB2J

## 📚 Tudo Que Você Aprendeu

### A Cadeia Completa

```
Quantum Espresso (NSCF)
    ↓ kmesh = [8, 8, 8] (calcula bandas em 512 k-points)
Wannier90
    ↓ mp_grid = [8, 8, 8] (herda do QE)
    ↓ (gera R-vectors até R_max ≈ 4)
wannier90_hr.dat
    ↓ H(R) para R ∈ [-4,-4,-4] até [4,4,4]
TB2J (integração Green)
    ↓ kmesh_integração ≤ 2×R_max ≈ [8,8,8]
    ↓ Se forçar kmesh > 8: ALIASING!
Exchange Parameters J(R)
    └─ Truncado em R_max=4 (série de Fourier incompleta)
```

---

## 🔑 Pontos-Chave Entendidos

| Conceito | Realidade | Implicação |
|----------|-----------|-----------|
| **mp_grid = kmesh QE?** | ✅ SIM | Wannier90 herda do QE |
| **R_max é automático?** | ✅ SIM | R_max ≈ mp_grid/2 |
| **search_shells afeta R_max?** | ❌ NÃO | Só afeta otimização |
| **TB2J pode ignorar R_max?** | ❌ NÃO | Causa aliasing |
| **kmesh_TB2J ≤ mp_grid?** | ✅ SIM | Regra fundamental |
| **J(R) depende de mp_grid?** | ✅ SIM | É série truncada |
| **Convergência requere novo QE?** | ✅ SIM | Se quer R_max > atual |

---

## 🎯 Seu Caso Específico

### Situação Atual
```
mp_grid QE/W90:     [8, 8, 8]
R_max do hr.dat:    4
kmesh seguro TB2J:  ≤ 8
```

### Flutuações Anteriores
```
Quando tentou kmesh = [12, 12, 12]:
├─ Estava tentando 1728 k-points
├─ Mas só tinha 729 R-points
├─ Resultado: ALIASING (não convergência!)
└─ Isto era esperado e inevitável
```

### Para Obter Convergência Real
```
Opção 1: Aceitar mp_grid=[8,8,8]
├─ Use kmesh_TB2J ≤ 8
├─ J(R) é válido mas truncado em R_max=4
└─ Publicável se for honesto sobre limitação

Opção 2: Refazer QE com mp_grid maior
├─ QE NSCF com kmesh = [12, 12, 12] ou [16, 16, 16]
├─ Novo cálculo Wannier90
├─ Novo hr.dat com R_max ≈ 6-8
└─ Aí TB2J pode integrar com kmesh maior
```

---

## 💡 Insights Importantes

### 1. Não é Interpolação do Wannier
- TB2J calcula H(k) exatamente via Fourier discreta
- Não interpola nada
- O limite é puramente matemático (FFT)

### 2. O Limite é Teórico
- Não é "falta de otimização" 
- É limite fundamental da transformada de Fourier discreta
- kmesh > 2×R_max viola teorema de Nyquist

### 3. Exchange Não Converge Arbitrariamente
- Não converge só aumentando kmesh em TB2J
- Precisa aumentar quantidade de informação em espaço real
- Isto só vem de maior mp_grid no QE

---

## 🚀 Próximos Passos Recomendados

### Para Validar Seu Entendimento

```python
# Verificar seu R_max atual
from TB2J.wannier import parse_ham
import numpy as np

n_wann, H_mnR, R_degens = parse_ham('wannier90_hr.dat')
Rlist = np.array(list(H_mnR.keys()))

print(f"R_min: {np.min(Rlist, axis=0)}")  
print(f"R_max: {np.max(Rlist, axis=0)}")  
print(f"Total R-points: {len(H_mnR)}")

# Esperado:
# R_min: [-4 -4 -4]
# R_max: [ 4  4  4]
# Total R-points: 729 (= 9³)
```

### Para Melhorar Convergência (Opcional)

```python
# Teste convergência de J(R) com kmesh
for kmesh_val in [4, 6, 8]:
    kmesh = [kmesh_val, kmesh_val, kmesh_val]
    green = TBGreen(tbmodel, kmesh=kmesh)
    
    # Calcular J
    J = compute_exchange(green, ...)
    print(f"kmesh=[{kmesh_val}³]: J(1,1,0) = {J:.4f} eV")

# Observação:
# - kmesh=4 até 8 deve convergir
# - kmesh > 8 terá artefatos
```

---

## 📖 Referências Teóricas

### Teorema de Nyquist-Shannon
- Se função tem frequência máxima f_max
- Precisa amostrar em frequência ≥ 2×f_max
- **Em TB2J:** Se R_max=4, kmesh ≤ 8

### Transformada de Fourier Discreta
- Forward (QE): bands em 512 k-points → H(R) em 729 R-points
- Backward (TB2J): G(R) requer kmesh correspondente (~8)
- Mismatch causa aliasing

---

## ✅ Conclusões

1. **Seu problema tinha causa real:** Você estava além de limite matemático
2. **Não era bug do TB2J:** Era uso incorreto da ferramenta
3. **Agora sabe:** kmesh_TB2J ≤ mp_grid_QE é regra, não sugestão
4. **Para convergência real:** Precisa de cálculo QE maior (não dá pra contornar)

---

## 🎓 Lições Aprendidas

```
✅ TB2J já faz Fourier discreta exata (não interpola)
✅ R_max vem de mp_grid via algoritmo automático W90
✅ search_shells afeta otimização, não estrutura de dados
✅ kmesh_integração é limitado por R_max (FFT Nyquist)
✅ Exchange é série de Fourier truncada em R_max
✅ Convergência real requer novo cálculo eletrônico (QE)
✅ Flutuações com kmesh >> mp_grid é aliasing esperado
```

---

## 📝 Checklist para Seu Projeto

- [ ] Verificar R_max do seu hr.dat atual
- [ ] Confirmar que kmesh ≤ 8 em TB2J funciona sem artefatos
- [ ] Decidir: aceitar truncamento em R_max=4 ou refazer QE?
- [ ] Se refazer QE: escolher novo mp_grid (12, 14, ou 16?)
- [ ] Documentar em seu paper: "Convergido para mp_grid=[8,8,8]"
- [ ] Se diferente em versão final: indicar novo mp_grid usado

---

Você agora entende completamente a física e matemática por trás! 🎯
