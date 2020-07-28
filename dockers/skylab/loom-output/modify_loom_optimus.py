import argparse
import loompy
def create_loom_files(args):
    """This function creates the loom file or folder structure in output_loom_path in format file_format,
       with biomaterial_id from the input
    
    Args:
        args (argparse.Namespace): input arguments for the run
    """
    ds = loompy.connect(args.input_loom_path)
    ds.attrs['biomaterial_id'] = args.biomaterial_id
    ds.close()

def main():
    description = """This script modifies the Optimus loom output and adds new attributes to the file.
                   This script can be used as a module or run as a command line script."""

    parser = argparse.ArgumentParser(description=description)

    parser.add_argument(
        "--input_path_for_loom",
        dest="input_loom_path",
        required=True,
        help="path to .loom file is to be created",
    )

    parser.add_argument(
        "--biomaterial_id",
        dest="biomaterial_id",
        default="Unknown sample",
        help="the biomaterial_id is cell_suspension.biomaterial_id defined by the user",
    )

    parser.add_argument(
        "--verbose",
        dest="verbose",
        action="store_true",
        help="whether to output verbose debugging messages",
    )
    
    args = parser.parse_args()

    create_loom_files(args)

if __name__ == "__main__":
    main()

