#pragma once

#include "../mmcfilters/trees/TreeOfShapesProducer.hpp"
#include "../mmcfilters/trees/ValuedMorphologicalTree.hpp"

#include <cstdint>
#include <memory>
#include <stdexcept>
#include <string>
#include <utility>
#include <variant>

namespace mmcfilters::pybindings {

/**
 * @brief Owning Python boundary for the supported valued-tree altitude types.
 *
 * Python exposes one semantic `ValuedMorphologicalTree` type. This adapter
 * retains the concrete C++ altitude type so ordinary uint8 hierarchies remain
 * uint8 while trees of shapes preserve exact doubled-unit `ToSGrayLevel`
 * altitudes.
 */
class PythonValuedMorphologicalTree {
  public:
    /** @brief Native uint8 valued-tree specialization. */
    using UInt8Tree = ValuedMorphologicalTree<std::uint8_t>;
    /** @brief Native exact tree-of-shapes valued-tree specialization. */
    using TreeOfShapesTree = ValuedMorphologicalTree<ToSGrayLevel>;
    /** @brief Owning runtime storage for supported altitude types. */
    using Storage = std::variant<std::shared_ptr<UInt8Tree>, std::shared_ptr<TreeOfShapesTree>>;

  private:
    Storage storage_; ///< Active native tree owner.

    /** @brief Rejects a null native tree. @tparam Tree Concrete tree type. @param tree Candidate owner. @return Valid owner. */
    template <class Tree> static std::shared_ptr<Tree> requireTree(std::shared_ptr<Tree> tree) {
        if (!tree) {
            throw std::invalid_argument("ValuedMorphologicalTree requires a non-null native tree.");
        }
        return tree;
    }

  public:
    /** @brief Owns a uint8 valued tree. @param tree Native tree owner. */
    explicit PythonValuedMorphologicalTree(std::shared_ptr<UInt8Tree> tree) : storage_(requireTree(std::move(tree))) {}

    /** @brief Owns an exact tree-of-shapes valued tree. @param tree Native tree owner. */
    explicit PythonValuedMorphologicalTree(std::shared_ptr<TreeOfShapesTree> tree) : storage_(requireTree(std::move(tree))) {}

    /** @brief Borrows a uint8 valued tree through a non-owning shared pointer. @param tree Native tree to borrow. */
    explicit PythonValuedMorphologicalTree(UInt8Tree& tree)
        : storage_(std::shared_ptr<UInt8Tree>(std::addressof(tree), [](UInt8Tree*) {})) {}

    /** @brief Borrows a constant uint8 valued tree. @param tree Native tree to borrow. */
    explicit PythonValuedMorphologicalTree(const UInt8Tree& tree) : PythonValuedMorphologicalTree(const_cast<UInt8Tree&>(tree)) {}

    /** @brief Visits the mutable active owner. @tparam Visitor Callable type. @param visitor Callable to invoke. @return Callable result. */
    template <class Visitor> decltype(auto) visit(Visitor&& visitor) {
        return std::visit([&](auto& tree) -> decltype(auto) { return std::forward<Visitor>(visitor)(tree); }, storage_);
    }

    /** @brief Visits the constant active owner. @tparam Visitor Callable type. @param visitor Callable to invoke. @return Callable result. */
    template <class Visitor> decltype(auto) visit(Visitor&& visitor) const {
        return std::visit([&](const auto& tree) -> decltype(auto) { return std::forward<Visitor>(visitor)(tree); }, storage_);
    }

    /** @brief Returns the shared topology. @return Active tree topology. */
    [[nodiscard]] const MorphologicalTree& topology() const noexcept {
        return visit([](const auto& tree) -> const MorphologicalTree& { return tree->topology(); });
    }

    /** @brief Requires the uint8 specialization. @param context Diagnostic operation name. @return uint8 tree owner. */
    [[nodiscard]] std::shared_ptr<UInt8Tree> requireUInt8Tree(const char* context) const {
        if (const auto* tree = std::get_if<std::shared_ptr<UInt8Tree>>(&storage_)) {
            return *tree;
        }
        throw std::invalid_argument(std::string(context) + " requires a uint8 component tree.");
    }

    /** @brief Prunes one node from the active tree. @param nodeId Live non-root node. */
    void pruneNode(NodeId nodeId) { visit([&](auto& tree) { tree->pruneNode(nodeId); }); }

    /** @brief Merges one node into its parent. @param nodeId Live non-root node. */
    void mergeNodeIntoParent(NodeId nodeId) { visit([&](auto& tree) { tree->mergeNodeIntoParent(nodeId); }); }
};

} // namespace mmcfilters::pybindings
