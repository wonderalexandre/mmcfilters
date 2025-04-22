

# Agora importa a versão e outros módulos
from .version import __version__

 
# Importar o módulo Pybind11 nativo
import mmcfilters

# Expor as classes do módulo Pybind11
MorphologicalTree = mmcfilters.MorphologicalTree
AdjacencyRelation = mmcfilters.AdjacencyRelation
AttributeFilters = mmcfilters.AttributeFilters
AttributeOpeningPrimitivesFamily = mmcfilters.AttributeOpeningPrimitivesFamily
IteratorNodesDescendants = mmcfilters.IteratorNodesDescendants
IteratorNodesOfPathToRoot = mmcfilters.IteratorNodesOfPathToRoot
IteratorPixelsOfCC = mmcfilters.IteratorPixelsOfCC
NodeMT = mmcfilters.NodeMT
ResidualTree = mmcfilters.ResidualTree
UltimateAttributeOpening = mmcfilters.UltimateAttributeOpening
Attribute = mmcfilters.Attribute
# Expor tudo no pacote
__all__ = [
    '__version__',
    # Classes Pybind11
    'MorphologicalTree',
    'AdjacencyRelation',
    'AttributeFilters',
    'AttributeOpeningPrimitivesFamily',
    'IteratorNodesDescendants',
    'IteratorNodesOfPathToRoot',
    'IteratorPixelsOfCC',
    'NodeMT',
    'Attribute',
    'ResidualTree',
    'UltimateAttributeOpening',
]
