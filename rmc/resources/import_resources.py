import argparse

from gnomad.resources.import_resources import (
    get_module_importable_resources,
    get_resources_descriptions,
)
import rmc.resources.grch37 as grch37
import rmc.resources.grch38 as grch38


grch37_resources = get_module_importable_resources(grch37, "grch37")
grch38_resources = get_module_importable_resources(grch38, "grch38")
all_resources = {**grch37_resources, **grch38_resources}


def main(args):
    for resource_arg in args.resources:
        resource_name, resource = all_resources[resource_arg]
        print(f"Importing {resource_name}...")
        resource.import_resource(args.overwrite)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(formatter_class=argparse.RawTextHelpFormatter)
    parser.add_argument(
        "resources",
        choices=list(all_resources.keys()),
        metavar="resource",
        nargs="+",
        help="Resource to import. Choices are:\n\n"
        + get_resources_descriptions(all_resources),
    )
    parser.add_argument(
        "--overwrite", help="Overwrites existing files", action="store_true"
    )
    main(parser.parse_args())
