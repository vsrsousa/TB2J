#!/usr/bin/env python
"""
Inspect TB2J.pickle to understand its structure
"""
import pickle

work_dir = "/Users/vinicius/Documents/gdni4al/TB2J_results_Gd_rcut10_kpt9"

print("=" * 90)
print("INSPECTING TB2J.pickle")
print("=" * 90)

with open(f"{work_dir}/TB2J.pickle", "rb") as f:
    data = pickle.load(f)

print(f"\nType: {type(data)}")

if isinstance(data, dict):
    print(f"Dictionary with {len(data)} keys:")
    for i, key in enumerate(list(data.keys())[:20]):
        print(f"  {i+1}. {key}: {type(data[key])}")
    
    # Check for magnetic structure
    if 'spinham' in data:
        spinham = data['spinham']
        print(f"\nFound 'spinham':")
        print(f"  Type: {type(spinham)}")
        if hasattr(spinham, 'cell'):
            print(f"  cell: {spinham.cell.shape}")
        if hasattr(spinham, 'pos'):
            print(f"  pos: {spinham.pos.shape}")
        if hasattr(spinham, 'zion'):
            print(f"  zion: {spinham.zion}")
    
    # Check for other relevant keys
    relevant_keys = ['cell', 'pos', 'positions', 'structure', 'atoms']
    for key in relevant_keys:
        if key in data:
            print(f"\nFound '{key}': {type(data[key])}")

elif hasattr(data, 'spinham'):
    print(f"\nObject with spinham attribute")
    spinham = data.spinham
    if hasattr(spinham, 'cell'):
        print(f"  spinham.cell: {spinham.cell.shape}")
