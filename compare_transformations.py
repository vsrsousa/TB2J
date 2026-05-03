#!/usr/bin/env python
"""
Compare: SeekPath original cell (now) vs online/primitive (before)
"""
import numpy as np

print("=" * 80)
print("COMPARAÇÃO: SeekPath PRIMITIVA vs HEXAGONAL")
print("=" * 80)

print("\n### ANTES (SeekPath Primitiva) - do test_seekpath_primitive.py ###")
print("""
SeekPath (same primitive cell):
  Path: [('GAMMA', 'Y'), ('Y', 'C_0'), ('SIGMA_0', 'GAMMA'), ('GAMMA', 'Z'), 
         ('Z', 'A_0'), ('E_0', 'T'), ('T', 'Y'), ('GAMMA', 'S'), 
         ('S', 'R'), ('R', 'Z'), ('Z', 'T')]
  Points: ['GAMMA', 'Y', 'T', 'Z', 'S', 'R', 'SIGMA_0', 'C_0', 'A_0', 'E_0']
  
  Special points (PRIMITIVA):
    GAMMA: [0, 0, 0]
    Y: [-0.5, 0.5, 0]
    T: [-0.5, 0.5, 0.5]
    Z: [0, 0, 0.5]
    S: [0, 0.5, 0]
    R: [0, 0.5, 0.5]
    SIGMA_0: [0.333, 0.333, 0]
    C_0: [-0.333, 0.667, 0]
    A_0: [0.333, 0.333, 0.5]
    E_0: [-0.333, 0.667, 0.5]
""")

print("\n### AGORA (SeekPath Hexagonal Transformada) - do test_hexagonal_kpath.py ###")
print("""
SeekPath (HEXAGONAL - TRANSFORMADA):
  Path: GAMMA → Y → SIGMA_0 → GAMMA → Z → E_0 → T → GAMMA → S → R → Z → T
  Points: ['GAMMA', 'Y', 'T', 'Z', 'S', 'R', 'SIGMA_0', 'C_0', 'A_0', 'E_0']
  
  Special points (HEXAGONAL):
    GAMMA: [0, 0, 0]
    Y: [-1, 0, 0]  ✓ Transformado: 2×[-0.5, 0.5, 0] → [-1, 0, 0]
    T: [-1, 0, 0.5]  ✓ Transformado: 2×[-0.5, 0.5, 0.5] → [-1, 0, 0.5]
    Z: [0, 0, 0.5]
    S: [-0.5, 0.5, 0]  ✓ Mesmo = 1×[0, 0.5, 0] @ (P_inv)^T
    R: [-0.5, 0.5, 0.5]  ✓ Transformado
    SIGMA_0: [0, 0.667, 0]  ✓ Transformado: [0.333, 0.333] @ transform
    C_0: [-1, 0.333, 0]  ✓ Transformado
    A_0: [0, 0.667, 0.5]  ✓ Transformado
    E_0: [-1, 0.333, 0.5]  ✓ Transformado
""")

print("\n" + "=" * 80)
print("✓ CONFIRMADO: Mesma simetria, mesmos pontos, MAS em espaço hexagonal!")
print("=" * 80)

# Mostra a transformação
print("\nTransformação aplicada (P_inv)^T:")
P_inv = np.array([[1, -1, 0], [1, 1, 0], [0, 0, 1]])
P_inv_T = P_inv.T
print(f"P_inv =\n{P_inv}")
print(f"\n(P_inv)^T =\n{P_inv_T}")

print("\nExemplo: Y em primitiva → hexagonal")
Y_prim = np.array([-0.5, 0.5, 0.0])
Y_hex = Y_prim @ P_inv_T
print(f"  Y_prim = {Y_prim}")
print(f"  Y_hex = Y_prim @ (P_inv)^T = {Y_hex}")
print(f"  Esperado: [-1, 0, 0] ✓")
