"""Small public example of attribute filtering from Python."""

import mmcfilters
import numpy as np


def main() -> None:
    image = np.array(
        [
            [3, 3, 2, 2],
            [3, 4, 4, 2],
            [1, 4, 5, 2],
            [1, 1, 5, 0],
        ],
        dtype=np.uint8,
    )

    tree = mmcfilters.MorphologicalTreeFactory.create_max_tree(image, radius=1.5)
    area = mmcfilters.Attribute.compute(
        tree,
        mmcfilters.Attribute.AREA,
    )
    filters = mmcfilters.AttributeFilters(tree)

    keep_large = area >= 4.0
    direct = filters.filtering_by_direct_rule(keep_large)
    pruning_min = filters.filtering_by_pruning_min_rule(area, 4.0)
    pruning_max = filters.filtering_by_pruning_max_rule(area, 4.0)

    print("Direct filter:\n", direct)
    print("Pruning-min rule:\n", pruning_min)
    print("Pruning-max rule:\n", pruning_max)


if __name__ == "__main__":
    main()
