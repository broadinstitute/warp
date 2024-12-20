import argparse
import json
import os
import ast

def update_test_inputs(inputs_json, truth_path, results_path, update_truth, branch_name):
    with open(inputs_json, 'r') as file:
        test_inputs = json.load(file)

    # Get the sample name from the test inputs JSON
    sample_name = os.path.splitext(os.path.basename(inputs_json))[0]

    # Get the pipeline name from the test inputs JSON
    pipeline_name = next(iter(test_inputs)).split('.')[0]

    # Append "Test" in front of the pipeline name
    test_name = f"Test{pipeline_name}"

    # Update all keys and ensure arrays are preserved
    updated_inputs = {}
    for key, value in test_inputs.items():
        new_key = key.replace(pipeline_name, test_name)

        # Handle the case where value might be a string representation of a list
        if isinstance(value, list):
            # Check if any element in the list is a string representation of another list
            processed_value = []
            for item in value:
                if isinstance(item, str) and item.startswith('[') and item.endswith(']'):
                    try:
                        # Use ast.literal_eval to safely evaluate string representation of list
                        inner_list = ast.literal_eval(item)
                        processed_value.extend(inner_list)
                    except (ValueError, SyntaxError):
                        processed_value.append(item)
                else:
                    processed_value.append(item)
            updated_inputs[new_key] = processed_value
        else:
            updated_inputs[new_key] = value

    # Add the truth_path and results_path to the updated inputs
    updated_inputs[f"{test_name}.results_path"] = f"{results_path}/{sample_name}/"
    updated_inputs[f"{test_name}.truth_path"] = f"{truth_path}/{sample_name}/"
    updated_inputs[f"{test_name}.update_truth"] = update_truth

    # Save the updated test inputs JSON
    output_name = f"updated_{sample_name}_{branch_name}.json"
    with open(output_name, 'w') as file:
        json.dump(updated_inputs, file, indent=4)

    print(f"{output_name}")
    return output_name

def main():
    description = """This script updates the test inputs JSON to work with the test wrapper WDL,
    which runs the pipeline and verification"""

    parser = argparse.ArgumentParser(description=description)

    parser.add_argument(
        "--truth_path",
        dest="truth_path",
        required=True,
        help="The base path where the truth data is stored",
    )

    parser.add_argument(
        "--results_path",
        dest="results_path",
        required=True,
        help="The base path where the test data will be stored",
    )

    parser.add_argument(
        "--inputs_json",
        dest="inputs_json",
        required=True,
        help="The JSON file containing the test inputs, formatted to run the pipeline WDL. "
             "This will be updated to run the wrapper Test WDL",
    )

    parser.add_argument(
        "--update_truth",
        dest="update_truth",
        default="false",
        required=False,
        choices=["true", "false"],
        help="Boolean flag to update the truth data. If true, the truth data will be updated with the test data. ",
    )

    parser.add_argument(
        "--branch_name",
        required=True,
        help="Branch name of the current pipeline run")

    args = parser.parse_args()
    # convert the update_truth flag to a boolean
    update_truth_bool = args.update_truth.lower() == "true"

    # Update the test inputs to work with the test wrapper WDL
    update_test_inputs(args.inputs_json, args.truth_path, args.results_path, update_truth_bool, args.branch_name)


if __name__ == "__main__":
    main()