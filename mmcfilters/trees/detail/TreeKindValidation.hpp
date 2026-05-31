#pragma once

#include "../MorphologicalTree.hpp"

#include <stdexcept>
#include <string>

namespace mmcfilters::detail {

inline const char* treeKindName(MorphologicalTreeKind kind) noexcept {
    switch (kind) {
        case MorphologicalTreeKind::MAX_TREE:
            return "MAX_TREE";
        case MorphologicalTreeKind::MIN_TREE:
            return "MIN_TREE";
        case MorphologicalTreeKind::TREE_OF_SHAPES:
            return "TREE_OF_SHAPES";
        case MorphologicalTreeKind::SELF_DUAL_RESIDUAL_TREE:
            return "SELF_DUAL_RESIDUAL_TREE";
    }
    return "UNKNOWN";
}

inline bool isComponentTreeKind(MorphologicalTreeKind kind) noexcept {
    return kind == MorphologicalTreeKind::MAX_TREE || kind == MorphologicalTreeKind::MIN_TREE;
}

inline void validateComponentTreeKind(const MorphologicalTree& tree, const char* context) {
    const MorphologicalTreeKind kind = tree.getTreeType();
    if (!isComponentTreeKind(kind)) {
        throw std::invalid_argument(
            std::string(context) +
            " is defined only for MAX_TREE and MIN_TREE; got " +
            treeKindName(kind) + ".");
    }
}

} // namespace mmcfilters::detail
