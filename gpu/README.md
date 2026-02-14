# WebGPU Barnes-Hut N-Body

This demo implements a Barnes-Hut–style N-body simulation on the GPU using a fixed quadtree hierarchy. It runs entirely in WebGPU compute shaders and renders 10k particles with a GPU instanced quad renderer.

## Features

- 10,000+ particles in a 2D periodic domain (providing your graphic card can handle it)
- Barnes-Hut approximation using a fixed multilevel quadtree
- GPU-only mass aggregation, force evaluation, and integration
- Mouse/touch pan and zoom

## Requirements

- A browser with WebGPU enabled (Chrome/Edge).
- Serve the demo over HTTP (local file URLs won’t reliably access WebGPU).

## Run

From the repo root, start a static server:

```bash
python -m http.server 8000
```

Then open:

```
http://localhost:8000/javascript/gpu/index.html
```

## Notes

This uses a fixed quadtree hierarchy (32→1) to avoid CPU tree construction. The force evaluation includes a Barnes-Hut style opening criterion and skips child cells when a parent is accepted.
