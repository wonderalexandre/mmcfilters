# MorphologicalAttributeFilters

MorphologicalAttributeFilters é uma biblioteca C++/Python para filtragem de imagens conectadas baseada em árvores morfológicas (component trees e trees of shapes). O código fornece implementações de alto desempenho e bindings Python para exploração interativa.

## Principais funcionalidades

* Construção de árvores morfológicas (component tree, tree of shapes) com diferentes conectividades.
* Cálculo incremental de atributos geométricos, topológicos e radiométricos.
* Filtros baseados em atributos (regras direta, subtrativa, poda, aberturas, ultimate attribute opening, etc.).
* Utilidades para valores de extinção, famílias de primitivas, MSER e Bit-Quads.
* Bindings Pybind11 cobrindo operações de alto nível no Python.

## Instalação

```bash
pip install mmcfilters
```

## Exemplo rápido

```python
import numpy as np
import mmcfilters

img = np.random.randint(0, 255, size=(128, 128), dtype=np.uint8)
tree = mmcfilters.MorphologicalTree(img.ravel(), img.shape[0], img.shape[1], True, 1.5)

attrs = mmcfilters.AttributeComputedIncrementally.computeSingleAttribute(tree, mmcfilters.Attribute.AREA)
area = np.asarray(attrs[1]).reshape(-1)

flt = mmcfilters.AttributeFilters(tree)
filtered = flt.filteringDirectRule(area > 50).reshape(img.shape)
```

## Documentação e comentários

Todos os `class` e `struct` do projeto possuem agora docstrings e comentários atualizados para facilitar a navegação do código C++. Consulte os headers em `mmcfilters/include` e `mmcfilters/pybind` para detalhes sobre cada componente.
