# Notebooks

The maintained notebooks are English, reader-facing examples of the current
public Python API. They are reproducible companions to the API guides, not a
compatibility layer or an alternative package loader.

The full Python contract is in [`docs/python-api.md`](../docs/python-api.md),
and the canonical scientific definitions for saliency are in
[`docs/saliency.md`](../docs/saliency.md).

## Prepare the environment

Install the package and optional notebook dependencies outside Jupyter, from
the repository root:

```bash
python -m pip install -e ".[notebooks]"
jupyter lab notebooks
```

Every notebook then imports the active environment normally:

```python
import mmcfilters
```

Notebook cells must not install packages, search a CMake build directory,
modify `sys.path`/import hooks, or load a checkout-specific extension. This
keeps the notebooks identical for editable installs, wheels, and CI.

## Authoring policy

- Maintain one English version of each notebook. Do not add parallel language
  copies or Portuguese titles, prose, labels, comments, and output text.
- Put `import mmcfilters` in the first code cell. Optional version output is
  allowed, but loader and installation logic is not.
- Use package-provided, synthetic, or repository-relative data. Never store
  `/Users/...`, `/Volumes/...`, or another workstation path in source or output.
- Tutorials may store bounded, sanitized outputs only after a clean top-to-bottom
  run. Experiment and diagnostic notebooks keep outputs cleared.
- Stored execution counts must cover every non-empty code cell and be exactly
  `1..N`. Error outputs are never committed.
- Keep algorithmic assertions near the result they verify. Saliency notebooks
  must distinguish formal edge maps, Cousty persistence, monotone LCA
  projection, raster `contourMap`, and Xu shape-space projection.

Run the static policy check before committing:

```bash
python scripts/validate_notebooks.py
```

Run the two small canonical smoke notebooks in isolated output copies:

```bash
python scripts/execute_notebooks.py --profile smoke
```

Use `--profile full` for the complete maintained suite. Pull requests execute
the smoke profile; the scheduled and opt-in manual CI profiles execute the full
suite.

## Recommended reading order and roles

| Order | Notebook | Mode and scope | Data | CI/output policy |
| ---: | --- | --- | --- | --- |
| 1 | `SimpleExamples.ipynb` | broad morphology introduction | repository images and synthetic fixtures | full / cleared |
| 2 | `Attribute_Filters.ipynb` | attribute computation and filter rules | repository images | full / cleared |
| 3 | `MaxDistExample.ipynb` | maximum-distance attribute tutorial | repository images | full / cleared |
| 4 | `Higra_MaxDist_Filtering.ipynb` | Higra interoperability diagnostic | repository images | full / cleared |
| 5 | `ToS_Contour_Example.ipynb` | Tree-of-Shapes contour tutorial | synthetic image | full / cleared |
| 6 | `UAO_Examples.ipynb` | Ultimate Attribute Opening tutorial | repository images | full / cleared |
| 7 | `Comparative_Morphological_Tree_Filtering.ipynb` | component, shape, and residual-tree comparison | synthetic image | full / cleared |
| 8 | `Filter.ipynb` | filtering, Cousty persistence, and Xu shaping comparison | repository image | full / cleared |
| 9 | `Saliency_Maps_Tutorial.ipynb` | canonical hierarchy-saliency tutorial | hand-computable synthetic image | smoke / clean stored output |
| 10 | `Hierarchy_Saliency_Map_Python_API.ipynb` | Python binding reference for hierarchy edge maps | `skimage.data.camera()` | smoke / cleared |
| 11 | `Hierarchy_Saliency_Map_Experiments.ipynb` | saliency-policy and projection diagnostics | `skimage.data.camera()` | full / cleared |
| 12 | `Extinction_Values_Experiments.ipynb` | extinction records, selection, persistence, and cuts | `skimage.data.camera()` | full / cleared |
| 13 | `Shape_Space_Saliency_Xu.ipynb` | Xu shaping tutorial | `skimage.data.coins()` | full / cleared |
| 14 | `Coins_Circularity_Contour_Saliency.ipynb` | contour-weighted QFZ experiment without shape space | `skimage.data.coins()` | full / cleared |

## Preferred tree model

Build public trees through `mmcfilters.MorphologicalTreeFactory`. Python image
factories currently require C-contiguous `np.uint8` arrays. Use the returned
`WeightedMorphologicalTree` for topology queries, altitude access,
reconstruction, and Higra export.

Common builders are:

- `createMaxTree(image, radius=1.5)`;
- `createMinTree(image, radius=1.5)`;
- `createTreeOfShapes(image)`.

Use `computeTopologyAttributes` for topology/support-only descriptors and
`computeSingleAttribute`/`computeAttributes` when weighted context is required.
Descriptor semantics are documented in
[`docs/attributes.md`](../docs/attributes.md).
