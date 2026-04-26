"""Compatibility helpers for legacy notebooks.

The notebooks in this directory predate the NodeId-first Python API. This module
keeps those notebooks executable while routing their calls through the current
public API.
"""

from __future__ import annotations

from dataclasses import dataclass
import sys
import types
from typing import Any

import numpy as np


def install(mmcfilters: Any) -> None:
    if getattr(mmcfilters, "_notebook_compat_installed", False):
        return

    native_attribute_filters = mmcfilters.AttributeFilters
    native_extinction_values = mmcfilters.ExtinctionValues
    native_aopf = mmcfilters.AttributeOpeningPrimitivesFamily
    native_uao = mmcfilters.UltimateAttributeOpening
    native_contours = mmcfilters.ContourComputation

    if "colorama" not in sys.modules:
        colorama = types.ModuleType("colorama")

        class _Color:
            BLACK = ""
            RED = ""
            GREEN = ""
            YELLOW = ""
            BLUE = ""
            MAGENTA = ""
            CYAN = ""
            WHITE = ""
            RESET = ""
            LIGHTBLACK_EX = ""
            LIGHTRED_EX = ""
            LIGHTGREEN_EX = ""
            LIGHTYELLOW_EX = ""
            LIGHTBLUE_EX = ""
            LIGHTMAGENTA_EX = ""
            LIGHTCYAN_EX = ""
            LIGHTWHITE_EX = ""

        class _Style:
            RESET_ALL = ""

        colorama.Fore = _Color
        colorama.Back = _Color
        colorama.Style = _Style
        colorama.__path__ = []
        sys.modules["colorama"] = colorama
        initialise = types.ModuleType("colorama.initialise")
        initialise.init = lambda *args, **kwargs: None
        initialise.deinit = lambda *args, **kwargs: None
        initialise.reinit = lambda *args, **kwargs: None
        sys.modules["colorama.initialise"] = initialise

    def unwrap(value: Any) -> Any:
        return value._native if isinstance(value, LegacyTree) else value

    def node_id(value: Any) -> int:
        return value.id if isinstance(value, LegacyNode) else int(value)

    @dataclass(frozen=True)
    class LegacyNode:
        tree: "LegacyTree"
        id: int

        @property
        def children(self) -> list["LegacyNode"]:
            return [LegacyNode(self.tree, child_id) for child_id in self.tree._native.getChildren(self.id)]

        @property
        def parent(self) -> "LegacyNode":
            return LegacyNode(self.tree, self.tree._native.getNodeParent(self.id))

        @property
        def level(self) -> int:
            return int(self.tree._native.getAltitude(self.id))

        @property
        def residue(self) -> int:
            return int(self.tree._native.getNodeResidue(self.id))

        @property
        def cnps(self) -> list[int]:
            return list(self.tree._native.getProperParts(self.id))

        @property
        def repNode(self) -> int:
            pixels = self.cnps or self.pixelsOfCC()
            return int(pixels[0]) if pixels else -1

        @property
        def repCNPs(self) -> list[int]:
            pixels = self.cnps or self.pixelsOfCC()
            return [int(pixels[0])] if pixels else []

        @property
        def area(self) -> float:
            return float(self.tree._area()[self.id])

        def recNode(self):
            return self.tree._native.reconstructNode(self.id)

        def pixelsOfCC(self) -> list[int]:
            pixels: list[int] = []
            for subtree_node_id in self.tree._native.getNodeSubtree(self.id):
                pixels.extend(self.tree._native.getProperParts(subtree_node_id))
            return pixels

        def postOrderTraversal(self) -> list["LegacyNode"]:
            return [LegacyNode(self.tree, node_id) for node_id in self.tree._native.getPostOrderNodes(self.id)]

        def __repr__(self) -> str:
            return f"Node({self.id})"

    class LegacyTree:
        def __init__(self, native_tree: Any):
            self._native = native_tree
            self._area_cache = None

        def _area(self):
            if self._area_cache is None:
                self._area_cache = mmcfilters.Attribute.computeSingleAttribute(self._native, mmcfilters.Attribute.AREA)
            return self._area_cache

        def __getattr__(self, name: str) -> Any:
            return getattr(self._native, name)

        @property
        def root(self) -> LegacyNode:
            return LegacyNode(self, self._native.getRoot())

        @property
        def listNodes(self) -> list[LegacyNode]:
            return [LegacyNode(self, node_id) for node_id in self._native.getAliveNodeIds()]

        @property
        def leaves(self) -> list[LegacyNode]:
            return [LegacyNode(self, node_id) for node_id in self._native.getLeafNodeIds()]

        def getSC(self, pixel_id: int) -> LegacyNode:
            return LegacyNode(self, self._native.getSmallestComponent(pixel_id))

        def getNode(self, node_id_value: Any) -> LegacyNode:
            return LegacyNode(self, node_id(node_id_value))

        def getChildren(self, node_id_value: Any):
            return self._native.getChildren(node_id(node_id_value))

        def getProperParts(self, node_id_value: Any):
            return self._native.getProperParts(node_id(node_id_value))

        def getAltitude(self, node_id_value: Any):
            return self._native.getAltitude(node_id(node_id_value))

        def getNodeResidue(self, node_id_value: Any):
            return self._native.getNodeResidue(node_id(node_id_value))

        def reconstructionImage(self):
            return self._native.reconstructionImage()

    def tree_of_shapes(image, interpolation=None):
        if interpolation is None:
            interpolation = mmcfilters.ToSInterpolation.SelfDual
        elif isinstance(interpolation, str):
            interpolation = mmcfilters.ToSInterpolation.Min4cMax8c if interpolation.lower() == "4c8c" else mmcfilters.ToSInterpolation.SelfDual
        return LegacyTree(mmcfilters.WeightedMorphologicalTree.createTreeOfShapes(image, interpolation))

    def component_tree(image, isMaxtree: bool = True, radius: float = 1.5):
        return LegacyTree(mmcfilters.WeightedMorphologicalTree.createComponentTree(image, bool(isMaxtree), float(radius)))

    def morphological_tree(image, *args, **kwargs):
        interpolation = kwargs.pop("interpolation", None)
        is_maxtree = kwargs.pop("isMaxtree", None)
        radius = kwargs.pop("radius", 1.5)

        if args:
            first = args[0]
            if isinstance(first, bool):
                is_maxtree = first
                if len(args) > 1:
                    radius = args[1]
            elif isinstance(first, str):
                interpolation = first
            else:
                interpolation = interpolation or mmcfilters.ToSInterpolation.SelfDual

        if is_maxtree is not None:
            return component_tree(image, is_maxtree, radius)
        return tree_of_shapes(image, interpolation)

    morphological_tree.createComponentTree = component_tree
    morphological_tree.createTreeOfShapes = tree_of_shapes
    morphological_tree.createFromAttributeMapping = lambda attr_map, image, isMaxtree=True, radius=1.5: component_tree(
        np.asarray(attr_map, dtype=np.uint8),
        bool(isMaxtree),
        float(radius),
    )
    morphological_tree.MAX_TREE = 0
    morphological_tree.MIN_TREE = 1
    morphological_tree.TREE_OF_SHAPES = 2
    mmcfilters.MorphologicalTree = morphological_tree

    original_compute_single = mmcfilters.Attribute.computeSingleAttribute
    original_compute_attrs = mmcfilters.Attribute.computeAttributes
    original_compute_delta = mmcfilters.Attribute.computeSingleAttributeWithDelta
    original_compute_mapping = mmcfilters.Attribute.computeAttributeMapping

    mmcfilters.Attribute.computeSingleAttribute = staticmethod(lambda tree, *args, **kwargs: original_compute_single(unwrap(tree), *args, **kwargs))
    mmcfilters.Attribute.computeAttributes = staticmethod(lambda tree, *args, **kwargs: original_compute_attrs(unwrap(tree), *args, **kwargs))
    mmcfilters.Attribute.computeSingleAttributeWithDelta = staticmethod(lambda tree, *args, **kwargs: original_compute_delta(unwrap(tree), *args, **kwargs))
    mmcfilters.Attribute.computeAttributeMapping = staticmethod(lambda tree, *args, **kwargs: original_compute_mapping(unwrap(tree), *args, **kwargs))
    mmcfilters.Attribute.computerAttributeMapping = staticmethod(lambda tree, *args, **kwargs: original_compute_mapping(unwrap(tree), *args, **kwargs))
    mmcfilters.Attribute.computeAttribute = staticmethod(lambda tree, attribute, *args, **kwargs: original_compute_single(unwrap(tree), attribute, *args, **kwargs))

    class AttributeFilters:
        def __init__(self, tree):
            self._native = native_attribute_filters(unwrap(tree))

        def __getattr__(self, name: str) -> Any:
            return getattr(self._native, name)

    class AttributeOpeningPrimitivesFamily:
        def __init__(self, tree, attr, *args):
            self._native = native_aopf(unwrap(tree), attr, *args)

        def __getattr__(self, name: str) -> Any:
            return getattr(self._native, name)

    class UltimateAttributeOpening:
        def __init__(self, tree, attr):
            self._native = native_uao(unwrap(tree), attr)

        def __getattr__(self, name: str) -> Any:
            return getattr(self._native, name)

    class ExtinctionValues:
        @dataclass(frozen=True)
        class Value:
            leaf: Any
            cutoffNode: Any
            extinction: Any

            def __iter__(self):
                yield self.leaf
                yield self.cutoffNode
                yield self.extinction

            def __getitem__(self, index: int):
                return (self.leaf, self.cutoffNode, self.extinction)[index]

        def __init__(self, tree, attr):
            self._tree = tree if isinstance(tree, LegacyTree) else None
            self._native = native_extinction_values(unwrap(tree), attr)

        def getExtinctionValues(self):
            values = self._native.getExtinctionValues()
            if self._tree is None:
                return values
            return [
                self.Value(LegacyNode(self._tree, leaf), LegacyNode(self._tree, cutoff), extinction)
                for leaf, cutoff, extinction in values
            ]

        def __getattr__(self, name: str) -> Any:
            return getattr(self._native, name)

    class ContourComputation:
        @staticmethod
        def extraction(tree):
            contours = native_contours.extraction(unwrap(tree))

            class Contours:
                def contoursByNode(self):
                    return contours.contoursByNode()

                def getContour(self, node_id_value):
                    return contours.getContour(node_id(node_id_value))

                def __getattr__(self, name: str) -> Any:
                    return getattr(contours, name)

            return Contours()

    mmcfilters.AttributeFilters = AttributeFilters
    mmcfilters.AttributeOpeningPrimitivesFamily = AttributeOpeningPrimitivesFamily
    mmcfilters.UltimateAttributeOpening = UltimateAttributeOpening
    mmcfilters.ExtinctionValues = ExtinctionValues
    mmcfilters.ContourComputation = ContourComputation
    mmcfilters._notebook_compat_installed = True
