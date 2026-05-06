"""Patch ASE k-point path generation to use SeekPath.

This module monkey-patches `ase.dft.kpoints.bandpath` and
the `bandpath` method on `ase.cell.Cell` to call SeekPath
when available. Importing this module applies the patch.

Usage: import TB2J.seekpath_patch  # done automatically from TB2J.__init__
"""
from __future__ import annotations

import warnings
import numpy as np
import math
from contextlib import contextmanager
from ase.dft.kpoints import kpoint_convert

try:
    import seekpath
except Exception:
    seekpath = None

try:
    import ase
    from ase.dft import kpoints as ase_kpoints
    from ase.cell import Cell
except Exception:
    ase = None
    ase_kpoints = None
    Cell = None


# Global context for structure information
_structure_context = {
    'cell': None,
    'positions': None,
    'atomic_numbers': None,
}


@contextmanager
def set_structure_context(cell, positions=None, atomic_numbers=None):
    """Context manager to set structure information for seekpath.
    
    Parameters
    ----------
    cell : array-like
        Unit cell vectors
    positions : array-like, optional
        Atomic positions
    atomic_numbers : array-like, optional
        Atomic numbers
    """
    old_context = _structure_context.copy()
    _structure_context['cell'] = np.array(cell) if cell is not None else None
    _structure_context['positions'] = np.array(positions) if positions is not None else None
    _structure_context['atomic_numbers'] = atomic_numbers
    try:
        yield
    finally:
        _structure_context.update(old_context)


def _dist(a, b, cell=None):
    """Distance between two k-points.
    
    Uses the proper metric (kpoint_convert) like ASE does, not simple Euclidean.
    If cell is provided, converts to Cartesian k-space; otherwise uses simple distance.
    """
    if cell is not None:
        # Convert k-points to Cartesian (reciprocal space with proper metric)
        a_cart = kpoint_convert(cell, skpts_kc=np.array(a))
        b_cart = kpoint_convert(cell, skpts_kc=np.array(b))
        return np.linalg.norm(b_cart - a_cart)
    else:
        # Fallback to simple Euclidean distance
        return math.sqrt(sum((a[i] - b[i]) ** 2 for i in range(3)))


def _normalize_kpoint_names(names):
    """Convert SeekPath k-point names to standard notation.
    
    Maps GAMMA→Γ, SIGMA_N→Σ$_N$, etc.
    """
    greek_map = {
        'GAMMA': 'Γ',
        'SIGMA': 'Σ',
        'DELTA': 'Δ',
        'LAMBDA': 'Λ',
    }
    
    normalized = []
    for name in names:
        result = name
        
        # Replace greek names
        for eng, greek in greek_map.items():
            if result.startswith(eng):
                remainder = result[len(eng):]
                if remainder.startswith('_'):
                    subscript = remainder[1:]
                    result = f"{greek}$_{subscript}$"
                else:
                    result = greek
                break
        
        # Handle other subscripted letters (X_0 → X$_0$)
        if '_' in result and not result.startswith('$'):
            parts = result.split('_')
            if len(parts) == 2 and parts[1].isdigit():
                result = f"{parts[0]}$_{parts[1]}$"
        
        # Convert G to Γ
        if result == 'G':
            result = 'Γ'
        
        normalized.append(result)
    return normalized


def _seekpath_bandpath_from_cell(cell, npoints=100, path=None):
    """Return a BandPath-compatible object using seekpath for a given cell.

    Parameters
    ----------
    cell : array-like
        3x3 unit cell vectors
    npoints : int
        Number of k-points along the path
    path : str or list, optional
        Path specification (for API compatibility, ignored)

    Returns
    -------
    object with get_linear_kpoint_axis() method and special_points attribute
        Compatible with ASE BandPath interface
    """
    if seekpath is None:
        warnings.warn("seekpath not available; falling back to ASE")
        return ase_kpoints.bandpath(None, cell, npoints)

    cell = np.array(cell, dtype=float)
    
    # Use structure context if available
    try:
        if (_structure_context['positions'] is not None and 
            _structure_context['atomic_numbers'] is not None):
            structure = (
                cell.tolist(),
                _structure_context['positions'].tolist(),
                _structure_context['atomic_numbers'],
            )
        else:
            # Dummy atom fallback
            structure = (cell.tolist(), [[0.0, 0.0, 0.0]], [1])
        
        data = seekpath.get_path(structure)
        point_coords = data.get("point_coords", {})
        path_segments = data.get("path", [])
    except Exception as e:
        warnings.warn(f"seekpath failed ({e}); falling back to ASE")
        return ase_kpoints.bandpath(None, cell, npoints)
    
    if not path_segments:
        warnings.warn("seekpath found no path; falling back to ASE")
        return ase_kpoints.bandpath(None, cell, npoints)
    
    # SeekPath returns k-points already in the original cell coordinates

    # Use path segments directly from SeekPath (already optimized)
    # No need to reorder - SeekPath handles the optimization
    contiguous_path = path_segments

    # Build k-points along the path
    all_kpts = []
    all_x = []
    cumulative_x = 0.0
    special_points_dict = {}
    special_point_x_list = []  # List of (label, x) - includes duplicates at breaks
    
    # First pass: calculate segment lengths to distribute npoints proportionally
    segment_lengths = []
    for start_label, end_label in contiguous_path:
        start_k = np.array(point_coords[start_label], dtype=float)
        end_k = np.array(point_coords[end_label], dtype=float)
        # Use proper metric for k-space distance
        d = _dist(end_k, start_k, cell=cell)
        segment_lengths.append(d)
    
    total_length = sum(segment_lengths)
    
    # Distribute npoints proportionally to segment lengths
    segment_npoints = []
    points_assigned = 0
    for i, length in enumerate(segment_lengths):
        if i == len(segment_lengths) - 1:
            # Last segment: use remaining points
            seg_npts = npoints - points_assigned
        else:
            # Proportional distribution
            seg_npts = max(2, int(round(npoints * length / total_length)))
            points_assigned += seg_npts
        segment_npoints.append(seg_npts)
    
    # Second pass: generate k-points with distributed counts
    for seg_idx, (start_label, end_label) in enumerate(contiguous_path):
        start_k = np.array(point_coords[start_label], dtype=float)
        end_k = np.array(point_coords[end_label], dtype=float)
        seg_npts = segment_npoints[seg_idx]
        
        # Record starting point position
        special_points_dict[start_label] = start_k
        
        # Add starting point (first segment or after a path break)
        if seg_idx == 0:
            special_point_x_list.append((start_label, cumulative_x))
        else:
            prev_end = contiguous_path[seg_idx - 1][1]
            if start_label != prev_end:
                # Path break: add this starting point as duplicate (prev_end already added as segment end)
                special_point_x_list.append((start_label, cumulative_x))
        
        # Generate k-points along this segment with distributed count
        # For contiguous segments, skip the first point to avoid duplication
        is_contiguous = (seg_idx > 0 and contiguous_path[seg_idx][0] == contiguous_path[seg_idx-1][1])
        
        if is_contiguous:
            # Generate with one extra point and skip the first one
            segment_kpts = np.array([
                start_k + (end_k - start_k) * t 
                for t in np.linspace(0, 1, seg_npts + 1)
            ])[1:]
        else:
            # Include the first point normally
            segment_kpts = np.array([
                start_k + (end_k - start_k) * t 
                for t in np.linspace(0, 1, seg_npts)
            ])
        
        # Calculate distances along segment (using proper metric with cell)
        if is_contiguous and seg_idx > 0:
            # For contiguous segments, first point connects to last point of previous segment
            prev_last_kpt = all_kpts[-1]
            first_kpt = segment_kpts[0]
            d_to_first = _dist(first_kpt, prev_last_kpt, cell=cell)
            segment_x = [cumulative_x + d_to_first]
            
            # Then continue with rest of points
            for i in range(1, len(segment_kpts)):
                d = _dist(segment_kpts[i], segment_kpts[i-1], cell=cell)
                segment_x.append(segment_x[-1] + d)
            segment_x = np.array(segment_x)
        else:
            # Normal case: starts at 0 and shifts by cumulative distance
            segment_x = [0.0]
            for i in range(1, len(segment_kpts)):
                d = _dist(segment_kpts[i], segment_kpts[i-1], cell=cell)
                segment_x.append(segment_x[-1] + d)
            segment_x = np.array(segment_x) + cumulative_x
        
        # Add all k-points from this segment
        all_kpts.extend(segment_kpts)
        all_x.extend(segment_x.tolist())
        
        # Record ending point
        special_points_dict[end_label] = end_k
        special_point_x_list.append((end_label, float(segment_x[-1])))
        
        cumulative_x = float(segment_x[-1])
    
    # Convert to numpy arrays
    all_kpts = np.array(all_kpts)
    all_x = np.array(all_x)
    
    # Build special_points dict with normalized names
    special_points = {}
    for k, v in special_points_dict.items():
        normalized_k = _normalize_kpoint_names([k])[0]
        special_points[normalized_k] = list(v)
    
    # Build path string like ASE does (e.g., "GMKGALH,LM,KH")
    path_string = ""
    for seg_idx, (start_label, end_label) in enumerate(contiguous_path):
        if seg_idx == 0:
            # First segment: add start + end
            path_string += start_label + end_label
        else:
            # Subsequent segments
            if contiguous_path[seg_idx-1][1] != start_label:
                # Path break - add comma and start + end
                path_string += ","
                path_string += start_label + end_label
            else:
                # Contiguous: only add end (start already added as previous end)
                path_string += end_label
    
    # Create a simple object that mimics BandPath
    class SeekPathBandPath:
        def __init__(self):
            self.kpts = all_kpts
            self.x_vals = all_x
            self.special_points = special_points
            self.path = path_string  # Add path attribute for compatibility
        
        def __iter__(self):
            """Compatibility with tuple unpacking like (kpts, x, X) = bandpath()"""
            x, xspecial, _ = self.get_linear_kpoint_axis()
            return iter([self.kpts, x, xspecial])
        
        def __getitem__(self, index):
            """Support indexing like bandpath()[0] and bandpath()[-1]"""
            x, xspecial, _ = self.get_linear_kpoint_axis()
            values = [self.kpts, x, xspecial]
            return values[index]
        
        def get_linear_kpoint_axis(self):
            # Return x-coords and special points (like ASE does)
            # Use special_point_x_list which has points in path order with correct cumulative X
            # NO duplicates at path breaks - the gap in cumulative_x is visualization enough
            knames = _normalize_kpoint_names([label for label, _ in special_point_x_list])
            Xs_array = np.array([x for _, x in special_point_x_list])
            return all_x, Xs_array, knames
    
    return SeekPathBandPath()



def build_custom_bandpath_from_seekpath(cell, point_names, npoints=100):
    """Build a custom band path from user-specified k-point names using SeekPath.
    
    Parameters
    ----------
    cell : array-like
        3x3 unit cell vectors
    point_names : list of str
        List of k-point labels, e.g., ['Gamma', 'X', 'P', 'N']
    npoints : int
        Total number of k-points to generate along the path
    
    Returns
    -------
    object with get_linear_kpoint_axis() and special_points attribute
        Compatible with ASE BandPath interface
    
    Raises
    ------
    ValueError
        If any point name is not found in SeekPath results or other validation fails
    """
    if seekpath is None:
        raise RuntimeError("SeekPath not available")
    
    cell = np.array(cell, dtype=float)
    
    # Get SeekPath data
    try:
        if (_structure_context['positions'] is not None and 
            _structure_context['atomic_numbers'] is not None):
            structure = (
                cell.tolist(),
                _structure_context['positions'].tolist(),
                _structure_context['atomic_numbers'],
            )
        else:
            structure = (cell.tolist(), [[0.0, 0.0, 0.0]], [1])
        
        data = seekpath.get_path(structure)
        point_coords = data.get("point_coords", {})
    except Exception as e:
        raise RuntimeError(f"SeekPath failed: {e}")
    
    # Validate all requested points exist (case-insensitive matching)
    # Create lowercase mapping to actual point names in seekpath
    point_coords_lower = {k.upper(): k for k in point_coords.keys()}
    
    # Map user input to actual point names
    mapped_point_names = []
    for name in point_names:
        name_upper = name.upper()
        if name_upper not in point_coords_lower:
            available = ", ".join(sorted(point_coords.keys()))
            raise ValueError(f"Point '{name}' not found. Available points: {available}")
        mapped_point_names.append(point_coords_lower[name_upper])
    
    point_names = mapped_point_names
    
    # Build segments from consecutive point pairs
    segments = []
    for i in range(len(point_names) - 1):
        segments.append((point_names[i], point_names[i+1]))
    
    if not segments:
        raise ValueError("Need at least 2 points to define a path")
    
    # Calculate segment lengths for proportional distribution of npoints
    segment_lengths = []
    for start_label, end_label in segments:
        start_k = np.array(point_coords[start_label], dtype=float)
        end_k = np.array(point_coords[end_label], dtype=float)
        d = _dist(end_k, start_k, cell=cell)
        segment_lengths.append(d)
    
    total_length = sum(segment_lengths)
    if total_length == 0:
        raise ValueError("All segments have zero length")
    
    # Distribute npoints proportionally
    segment_npoints = []
    points_assigned = 0
    for i, length in enumerate(segment_lengths):
        if i == len(segment_lengths) - 1:
            seg_npts = npoints - points_assigned
        else:
            seg_npts = max(2, int(round(npoints * length / total_length)))
        segment_npoints.append(seg_npts)
        points_assigned += seg_npts
    
    # Build k-points
    all_kpts = []
    all_x = []
    cumulative_x = 0.0
    special_point_x_list = []
    
    for seg_idx, (start_label, end_label) in enumerate(segments):
        start_k = np.array(point_coords[start_label], dtype=float)
        end_k = np.array(point_coords[end_label], dtype=float)
        seg_npts = segment_npoints[seg_idx]
        
        # First point of segment
        if seg_idx == 0:
            all_kpts.append(start_k)
            all_x.append(cumulative_x)
            special_point_x_list.append((start_label, cumulative_x))
        
        # Interpolate between start and end
        for i in range(1, seg_npts + 1):
            frac = i / seg_npts
            kpt = (1.0 - frac) * start_k + frac * end_k
            all_kpts.append(kpt)
            
            # x-coordinate (distance along path)
            dist = _dist(kpt, start_k, cell=cell)
            x = cumulative_x + dist
            all_x.append(x)
            
            # Record special point at end of segment
            if i == seg_npts:
                special_point_x_list.append((end_label, x))
                cumulative_x = x
    
    all_kpts = np.array(all_kpts)
    all_x = np.array(all_x)
    
    # Normalize k-point names
    knames = _normalize_kpoint_names([label for label, _ in special_point_x_list])
    Xs_array = np.array([x for _, x in special_point_x_list])
    
    # Create BandPath-compatible object with ORIGINAL point coordinates for special_points
    class CustomBandPath:
        def __init__(self, kpts, x, special_points_x, kpoint_names, original_point_coords, special_point_labels):
            self.kpts = kpts
            self.x = x
            self.special_points_x = special_points_x
            self.kpoint_names = kpoint_names
            # Use original coordinates from SeekPath, not interpolated
            self.special_points = {name: np.array(original_point_coords[label], dtype=float) 
                                  for name, label in zip(kpoint_names, special_point_labels)}
        
        def __iter__(self):
            return iter([self.kpts, self.x, self.special_points_x])
        
        def __getitem__(self, index):
            return [self.kpts, self.x, self.special_points_x][index]
        
        def get_linear_kpoint_axis(self):
            return self.x, self.special_points_x, self.kpoint_names
    
    # Get original labels (before normalization) for mapping
    original_labels = [label for label, _ in special_point_x_list]
    
    return CustomBandPath(all_kpts, all_x, Xs_array, knames, point_coords, original_labels)


def _apply_patch():
    """Apply monkey-patches to ASE to use SeekPath."""
    if ase_kpoints is None or seekpath is None:
        return

    _original_bandpath = ase_kpoints.bandpath

    def bandpath_wrapper(kpoints_or_path, cell, npoints=100):
        """Wrapper using SeekPath when possible."""
        try:
            # If explicit k-vectors, use original ASE
            if hasattr(kpoints_or_path, "__iter__") and not isinstance(kpoints_or_path, str):
                arr = np.array(kpoints_or_path)
                if len(arr.shape) > 1:
                    return _original_bandpath(kpoints_or_path, cell, npoints)
            
            # Use SeekPath for paths or automatic generation
            return _seekpath_bandpath_from_cell(cell, npoints=npoints)
        except Exception as e:
            warnings.warn(f"SeekPath failed: {e}; falling back to ASE")
            return _original_bandpath(kpoints_or_path, cell, npoints)

    ase_kpoints.bandpath = bandpath_wrapper

    # Patch Cell.bandpath method
    if Cell is not None:
        _original_cell_bandpath = Cell.bandpath

        def cell_bandpath_wrapper(self, path=None, npoints=100, **kwargs):
            """Wrapper using SeekPath."""
            try:
                cell_array = (np.array(self.get_cell()) if hasattr(self, "get_cell") 
                             else np.array(self))
                return _seekpath_bandpath_from_cell(cell_array, npoints=npoints)
            except Exception as e:
                warnings.warn(f"SeekPath failed: {e}; falling back to ASE")
                return _original_cell_bandpath(self, path=path, npoints=npoints, **kwargs)

        Cell.bandpath = cell_bandpath_wrapper


_apply_patch()
