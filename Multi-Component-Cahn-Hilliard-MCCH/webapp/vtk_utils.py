"""
High-performance VTK Structured Points Parser and Image Renderer for Cahn-Hilliard simulations.
Optimized for low-latency web serving with in-memory caching and fast PIL/NumPy generation.
"""

import os
import re
import io
import time
from collections import OrderedDict
import numpy as np
from PIL import Image, ImageDraw, ImageFont
import matplotlib
matplotlib.use('Agg')
import matplotlib.cm as cm

# In-memory LRU Cache for parsed VTK files to provide instantaneous scrubbing
_VTK_CACHE = OrderedDict()
_MAX_CACHE_ENTRIES = 32

def parse_vtk_header(filepath):
    """Quickly read dimensions and scalar names without parsing all float data."""
    if not os.path.exists(filepath):
        return None
    nx, ny, nz = 1, 1, 1
    fields = []
    with open(filepath, 'r') as f:
        for _ in range(30):
            line = f.readline()
            if not line:
                break
            line_str = line.strip()
            if line_str.startswith("DIMENSIONS"):
                parts = line_str.split()
                nx, ny, nz = int(parts[1]), int(parts[2]), int(parts[3])
            elif line_str.startswith("SCALARS"):
                parts = line_str.split()
                fields.append(parts[1])
    return {"nx": nx, "ny": ny, "nz": nz, "fields": fields}


def parse_vtk_file(filepath):
    """
    Fast VTK Structured Points ASCII parser using block string searching and NumPy fromstring.
    Cached in memory for rapid animation playback.
    """
    global _VTK_CACHE
    if filepath in _VTK_CACHE:
        # Move to end (LRU)
        _VTK_CACHE.move_to_end(filepath)
        return _VTK_CACHE[filepath]

    if not os.path.exists(filepath):
        return None

    with open(filepath, 'r') as f:
        content = f.read()

    # Find DIMENSIONS
    dim_idx = content.find('DIMENSIONS')
    if dim_idx == -1:
        return None
    dim_line = content[dim_idx:content.find('\n', dim_idx)]
    parts = dim_line.split()
    nx, ny, nz = int(parts[1]), int(parts[2]), int(parts[3])

    # Find step and time if present in header
    step = 0
    sim_time = 0.0
    time_match = re.search(r'MCCH Simulation Step (\d+) Time ([\d\.\+eE\-]+)', content[:500])
    if time_match:
        step = int(time_match.group(1))
        sim_time = float(time_match.group(2))

    fields = {}
    pos = 0
    while True:
        s_idx = content.find('SCALARS', pos)
        if s_idx == -1:
            break
        s_end = content.find('\n', s_idx)
        scalar_line = content[s_idx:s_end].split()
        name = scalar_line[1]

        # Find LOOKUP_TABLE
        lut_idx = content.find('LOOKUP_TABLE', s_end)
        if lut_idx == -1:
            break
        data_start = content.find('\n', lut_idx) + 1

        # Next SCALARS or end of file
        next_s = content.find('SCALARS', data_start)
        if next_s == -1:
            data_str = content[data_start:]
            pos = len(content)
        else:
            data_str = content[data_start:next_s]
            pos = next_s

        arr = np.fromstring(data_str, dtype=np.float64, sep=' ')
        # Shape: (nz, ny, nx)
        if arr.size == nx * ny * nz:
            fields[name] = arr.reshape((nz, ny, nx))

    result = {
        "filepath": filepath,
        "nx": nx,
        "ny": ny,
        "nz": nz,
        "dim": 3 if nz > 1 else 2,
        "step": step,
        "time": sim_time,
        "fields": fields,
        "field_names": list(fields.keys())
    }

    # Cache management
    if len(_VTK_CACHE) >= _MAX_CACHE_ENTRIES:
        _VTK_CACHE.popitem(last=False)
    _VTK_CACHE[filepath] = result

    return result


def extract_slice_2d(parsed_vtk, field_name, plane='xy', slice_idx=None):
    """
    Extracts a 2D slice from the parsed 3D/2D field.
    plane: 'xy', 'xz', or 'yz'
    Returns (2d_array, extent_info)
    """
    if not parsed_vtk or field_name not in parsed_vtk["fields"]:
        return None, None

    data_3d = parsed_vtk["fields"][field_name]
    nz, ny, nx = parsed_vtk["nz"], parsed_vtk["ny"], parsed_vtk["nx"]

    if nz == 1 or plane == 'xy':
        if slice_idx is None:
            slice_idx = nz // 2
        slice_idx = max(0, min(nz - 1, int(slice_idx)))
        # data_3d is (nz, ny, nx) -> slice is (ny, nx)
        slice_2d = data_3d[slice_idx, :, :]
        extent = {"width": nx, "height": ny, "xlabel": "X", "ylabel": "Y", "slice_pos": slice_idx, "max_slice": nz}
    elif plane == 'xz':
        if slice_idx is None:
            slice_idx = ny // 2
        slice_idx = max(0, min(ny - 1, int(slice_idx)))
        # data_3d[:, slice_idx, :] -> (nz, nx)
        slice_2d = data_3d[:, slice_idx, :]
        extent = {"width": nx, "height": nz, "xlabel": "X", "ylabel": "Z", "slice_pos": slice_idx, "max_slice": ny}
    elif plane == 'yz':
        if slice_idx is None:
            slice_idx = nx // 2
        slice_idx = max(0, min(nx - 1, int(slice_idx)))
        # data_3d[:, :, slice_idx] -> (nz, ny)
        slice_2d = data_3d[:, :, slice_idx]
        extent = {"width": ny, "height": nz, "xlabel": "Y", "ylabel": "Z", "slice_pos": slice_idx, "max_slice": nx}
    else:
        slice_2d = data_3d[0, :, :]
        extent = {"width": nx, "height": ny, "xlabel": "X", "ylabel": "Y", "slice_pos": 0, "max_slice": 1}

    return slice_2d, extent


def draw_slab_overlay(img, nx, num_ranks):
    """Draw vertical dashed lines representing MPI 1D slab decomposition boundaries."""
    if not num_ranks or num_ranks <= 1:
        return img

    w, h = img.size
    scale_x = w / nx
    base = nx // num_ranks
    rem = nx % num_ranks

    draw = ImageDraw.Draw(img)
    x_pos = 0

    # Overlay line style
    for r in range(num_ranks):
        loc_nx = base + (1 if r < rem else 0)
        x_pos += loc_nx
        if r < num_ranks - 1:
            px = int(x_pos * scale_x)
            # Draw dashed line
            dash_len = 6
            for y in range(0, h, dash_len * 2):
                draw.line([(px, y), (px, min(h, y + dash_len))], fill=(255, 60, 60, 220), width=2)
            # Rank label
            draw.text((px - 32, 8), f"R{r}|R{r+1}", fill=(255, 220, 0, 240))

    return img


def render_slice_png(data_2d, colormap_name='viridis', vmin=None, vmax=None, overlay_slabs=None, orig_nx=None, target_size=512):
    """
    Renders 2D numpy array directly to PNG bytes using matplotlib colormap and PIL.
    Ultra fast (~5-10ms).
    """
    if data_2d is None:
        return b""

    # Determine normalization bounds
    if vmin is None:
        vmin = float(np.min(data_2d))
    if vmax is None:
        vmax = float(np.max(data_2d))

    if abs(vmax - vmin) < 1e-12:
        normed = np.zeros_like(data_2d, dtype=np.float32)
    else:
        normed = np.clip((data_2d - vmin) / (vmax - vmin), 0.0, 1.0)

    # Get colormap
    try:
        cmap = matplotlib.colormaps.get_cmap(colormap_name)
    except Exception:
        try:
            cmap = cm.get_cmap(colormap_name)
        except Exception:
            cmap = cm.get_cmap('viridis')

    # Apply colormap (H, W, 4)
    rgba = (cmap(normed) * 255).astype(np.uint8)

    # Invert Y axis for standard scientific lower origin
    rgba_flipped = rgba[::-1, :, :]
    img = Image.fromarray(rgba_flipped, mode='RGBA')

    # Upscale smoothly with nearest/bilinear interpolation
    h, w = data_2d.shape
    if target_size and (w < target_size or h < target_size):
        scale = max(1, target_size // max(w, h))
        img = img.resize((w * scale, h * scale), resample=Image.NEAREST)

    # Draw slab overlay if requested
    if overlay_slabs and overlay_slabs > 1:
        domain_nx = orig_nx if orig_nx else w
        img = draw_slab_overlay(img, domain_nx, overlay_slabs)

    buf = io.BytesIO()
    img.save(buf, format='PNG', optimize=False)
    return buf.getvalue(), float(vmin), float(vmax)


def render_rgb_composite_png(c0_2d, c1_2d, c2_2d, overlay_slabs=None, orig_nx=None, target_size=512):
    """
    Renders RGB composite image for 3-component systems (R=c0, G=c1, B=c2).
    """
    if c0_2d is None or c1_2d is None or c2_2d is None:
        return b""

    r = np.clip(c0_2d * 255.0, 0, 255).astype(np.uint8)
    g = np.clip(c1_2d * 255.0, 0, 255).astype(np.uint8)
    b = np.clip(c2_2d * 255.0, 0, 255).astype(np.uint8)
    rgb = np.stack([r, g, b], axis=-1)

    # Invert Y axis for lower origin
    rgb_flipped = rgb[::-1, :, :]
    img = Image.fromarray(rgb_flipped, mode='RGB')

    h, w = c0_2d.shape
    if target_size and (w < target_size or h < target_size):
        scale = max(1, target_size // max(w, h))
        img = img.resize((w * scale, h * scale), resample=Image.NEAREST)

    if overlay_slabs and overlay_slabs > 1:
        domain_nx = orig_nx if orig_nx else w
        img = draw_slab_overlay(img, domain_nx, overlay_slabs)

    buf = io.BytesIO()
    img.save(buf, format='PNG', optimize=False)
    return buf.getvalue()


def extract_line_cut(parsed_vtk, axis='x', coord=None, z_slice=0):
    """
    Extracts 1D concentration cut across the computational domain.
    axis='x': horizontal line cut at fixed y=coord -> returns profiles along X [0..nx-1]
    axis='y': vertical line cut at fixed x=coord -> returns profiles along Y [0..ny-1]
    """
    if not parsed_vtk:
        return None

    nx, ny, nz = parsed_vtk["nx"], parsed_vtk["ny"], parsed_vtk["nz"]
    z = max(0, min(nz - 1, int(z_slice)))

    comp_keys = [k for k in sorted(parsed_vtk["fields"].keys()) if k.startswith('c')]

    if axis == 'x':
        # Cut at fixed Y
        if coord is None:
            coord = ny // 2
        y = max(0, min(ny - 1, int(coord)))
        coords = list(range(nx))
        profiles = {}
        for key in comp_keys:
            profiles[key] = parsed_vtk["fields"][key][z, y, :].tolist()
        return {"axis": "x", "cut_at_y": y, "z": z, "positions": coords, "profiles": profiles}
    else:
        # Cut at fixed X
        if coord is None:
            coord = nx // 2
        x = max(0, min(nx - 1, int(coord)))
        coords = list(range(ny))
        profiles = {}
        for key in comp_keys:
            profiles[key] = parsed_vtk["fields"][key][z, :, x].tolist()
        return {"axis": "y", "cut_at_x": x, "z": z, "positions": coords, "profiles": profiles}
