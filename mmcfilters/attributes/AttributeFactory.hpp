#pragma once

#include "AttributeComputer.hpp"
#include "AreaComputer.hpp"
#include "BoundingBoxComputer.hpp"
#include "GrayLevelStatsComputer.hpp"
#include "MomentBasedAttributeComputer.hpp"
#include "TreeTopologyComputer.hpp"
#include "VolumeComputer.hpp"
#include "BitquadsComputer.hpp"
#include "MaxDistComputer.hpp"

namespace mmcfilters {

/**
 * @brief Factory that maps public attribute requests to concrete computers.
 *
 * @details
 * The incremental attribute pipeline never hard-codes concrete computer types
 * outside this factory. `AttributeFactory` is therefore the bridge between:
 * - the public request space (`Attribute` and `AttributeGroup`);
 * - the internal execution space (`AttributeComputer` subclasses).
 *
 * Several scalar attributes map to the same computer because they share a
 * traversal or intermediate state. The factory captures that grouping policy
 * in one place, which keeps the orchestration layer independent from the
 * details of each implementation.
 */
class AttributeFactory {
private:
    template <class T>
    static const AttributeComputer& singleton() {
        static const T computer{};
        return computer;
    }

    /**
     * @brief Returns the concrete computer responsible for one scalar
     * attribute.
     */
    static const AttributeComputer& createImpl(Attribute attr) {
        switch (attr) {
            case AREA:
                return singleton<AreaComputer>();

            case RELATIVE_VOLUME:
            case VOLUME:
                return singleton<VolumeComputer>();

            case GRAY_HEIGHT:
            case LEVEL:
            case MEAN_LEVEL:
            case VARIANCE_LEVEL:
                return singleton<GrayLevelStatsComputer>();

            case BOX_COL_MIN:
            case BOX_COL_MAX:
            case BOX_ROW_MIN:
            case BOX_ROW_MAX:
            case RATIO_WH:
            case RECTANGULARITY:
            case DIAGONAL_LENGTH:
            case BOX_HEIGHT:
            case BOX_WIDTH:
                return singleton<BoundingBoxComputer>();

            case AXIS_ORIENTATION:
            case LENGTH_MAJOR_AXIS:
            case LENGTH_MINOR_AXIS:
            case ECCENTRICITY:
            case INERTIA:
            case COMPACTNESS:
            case CIRCULARITY:
                return singleton<MomentBasedAttributeComputer>();

            case CENTRAL_MOMENT_20:
            case CENTRAL_MOMENT_02:
            case CENTRAL_MOMENT_11:
            case CENTRAL_MOMENT_30:
            case CENTRAL_MOMENT_03:
            case CENTRAL_MOMENT_21:
            case CENTRAL_MOMENT_12:
                return singleton<CentralMomentsComputer>();

            case HU_MOMENT_1:
            case HU_MOMENT_2:
            case HU_MOMENT_3:
            case HU_MOMENT_4:
            case HU_MOMENT_5:
            case HU_MOMENT_6:
            case HU_MOMENT_7:
                return singleton<HuMomentsComputer>();

            case HEIGHT_NODE:
            case DEPTH_NODE:
            case IS_LEAF_NODE:
            case IS_ROOT_NODE:
            case NUM_CHILDREN_NODE:
            case NUM_SIBLINGS_NODE:
            case NUM_DESCENDANTS_NODE:
            case NUM_LEAF_DESCENDANTS_NODE:
            case LEAF_RATIO_NODE:
            case BALANCE_NODE:
            case AVG_CHILD_HEIGHT_NODE:
                return singleton<TreeTopologyComputer>();

            case BITQUADS_AREA:
            case BITQUADS_NUMBER_EULER:
            case BITQUADS_NUMBER_HOLES:
            case BITQUADS_PERIMETER:
            case BITQUADS_PERIMETER_CONTINUOUS:
            case BITQUADS_CIRCULARITY:
            case BITQUADS_PERIMETER_AVERAGE:
            case BITQUADS_LENGTH_AVERAGE:
            case BITQUADS_WIDTH_AVERAGE:
                return singleton<BitquadsComputer>();

            case MAX_DIST:
                return singleton<MaxDistComputer>();

            default:
                throw std::runtime_error("Attribute not supported.");
        }
    }

    /**
     * @brief Returns the concrete computer responsible for one public
     * attribute group.
     */
    static const AttributeComputer& createImpl(AttributeGroup group) {
        switch (group) {
            case AttributeGroup::BOUNDING_BOX:
                return singleton<BoundingBoxComputer>();
            case AttributeGroup::CENTRAL_MOMENTS:
                return singleton<CentralMomentsComputer>();
            case AttributeGroup::HU_MOMENTS:
                return singleton<HuMomentsComputer>();
            case AttributeGroup::MOMENT_BASED:
                return singleton<MomentBasedAttributeComputer>();
            case AttributeGroup::TREE_TOPOLOGY:
                return singleton<TreeTopologyComputer>();
            case AttributeGroup::BITQUADS:
                return singleton<BitquadsComputer>();
            default:
                throw std::runtime_error("Attribute group not supported.");
        }
    }

public:
    /**
     * @brief Returns the stateless singleton associated with a scalar
     * attribute or attribute group.
     */
    static const AttributeComputer& create(const AttributeOrGroup& attr) {
        return std::visit(
            [](auto&& actualAttr) -> const AttributeComputer& {
                return AttributeFactory::createImpl(actualAttr);
            },
            attr);
    }
};

} // namespace mmcfilters
