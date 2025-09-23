#pragma once

#include "AttributeComputer.hpp"
#include "AreaComputer.hpp"
#include "BoundingBoxComputer.hpp"
#include "GrayLevelStatsComputer.hpp"
#include "MomentBasedAttributeComputer.hpp"
#include "TreeTopologyComputer.hpp"
#include "VolumeComputer.hpp"
#include "BitquadsComputer.hpp"

namespace mmcfilters {

class AttributeFactory {
private:
    static std::shared_ptr<AttributeComputer> createImpl(Attribute attr) {
        switch (attr) {
            case AREA:
                return std::make_shared<AreaComputer>();

            case RELATIVE_VOLUME:
            case VOLUME:
                return std::make_shared<VolumeComputer>();

            case GRAY_HEIGHT:
            case LEVEL:
            case MEAN_LEVEL:
            case VARIANCE_LEVEL:
                return std::make_shared<GrayLevelStatsComputer>();

            case BOX_COL_MIN:
            case BOX_COL_MAX:
            case BOX_ROW_MIN:
            case BOX_ROW_MAX:
            case RATIO_WH:
            case RECTANGULARITY:
            case DIAGONAL_LENGTH:
            case BOX_HEIGHT:
            case BOX_WIDTH:
                return std::make_shared<BoundingBoxComputer>();

            case AXIS_ORIENTATION:
            case LENGTH_MAJOR_AXIS:
            case LENGTH_MINOR_AXIS:
            case ECCENTRICITY:
            case INERTIA:
            case COMPACTNESS:
            case CIRCULARITY:
                return std::make_shared<MomentBasedAttributeComputer>();

            case CENTRAL_MOMENT_20:
            case CENTRAL_MOMENT_02:
            case CENTRAL_MOMENT_11:
            case CENTRAL_MOMENT_30:
            case CENTRAL_MOMENT_03:
            case CENTRAL_MOMENT_21:
            case CENTRAL_MOMENT_12:
                return std::make_shared<CentralMomentsComputer>();

            case HU_MOMENT_1:
            case HU_MOMENT_2:
            case HU_MOMENT_3:
            case HU_MOMENT_4:
            case HU_MOMENT_5:
            case HU_MOMENT_6:
            case HU_MOMENT_7:
                return std::make_shared<HuMomentsComputer>();

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
                return std::make_shared<TreeTopologyComputer>();

            case BITQUADS_AREA:
            case BITQUADS_NUMBER_EULER:
            case BITQUADS_NUMBER_HOLES:
            case BITQUADS_PERIMETER:
            case BITQUADS_PERIMETER_CONTINUOUS:
            case BITQUADS_CIRCULARITY:
            case BITQUADS_PERIMETER_AVERAGE:
            case BITQUADS_LENGTH_AVERAGE:
            case BITQUADS_WIDTH_AVERAGE:
                return std::make_shared<BitquadsComputer>();

            default:
                throw std::runtime_error("Attribute not supported.");
        }
    }

    static std::shared_ptr<AttributeComputer> createImpl(AttributeGroup group) {
        switch (group) {
            case AttributeGroup::BOUNDING_BOX:
                return std::make_shared<BoundingBoxComputer>();
            case AttributeGroup::CENTRAL_MOMENTS:
                return std::make_shared<CentralMomentsComputer>();
            case AttributeGroup::HU_MOMENTS:
                return std::make_shared<HuMomentsComputer>();
            case AttributeGroup::MOMENT_BASED:
                return std::make_shared<MomentBasedAttributeComputer>();
            case AttributeGroup::TREE_TOPOLOGY:
                return std::make_shared<TreeTopologyComputer>();
            case AttributeGroup::BITQUADS:
                return std::make_shared<BitquadsComputer>();
            default:
                throw std::runtime_error("Attribute group not supported.");
        }
    }

public:
    static std::shared_ptr<AttributeComputer> create(const AttributeOrGroup& attr) {
        return std::visit(
            [](auto&& actualAttr) -> std::shared_ptr<AttributeComputer> {
                return AttributeFactory::createImpl(actualAttr);
            },
            attr);
    }
};

} // namespace mmcfilters

