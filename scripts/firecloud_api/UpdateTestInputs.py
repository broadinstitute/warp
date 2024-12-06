import argparse
import json
import os


def update_test_inputs(inputs_json, truth_path, results_path, update_truth, commit_hash):
    # Update the test inputs JSON to work with the test wrapper WDL
    # The test wrapper WDL runs the pipeline WDL and verifies the results
    # The test wrapper WDL requires the following inputs:
    # - truth_path: The path to the truth data
    # - results_path: The path to the test data
    # - update_truth: Boolean indicating whether truth should be updated, default is False
    # The test inputs JSON will be updated to include these inputs

    with open(inputs_json, 'r') as file:
        test_inputs = json.load(file)

    # get the sample name from the test inputs JSON, this is needed for tests with multiple inputs
    sample_name = os.path.splitext(os.path.basename(inputs_json))[0]

    # Get the pipeline name from the test inputs JSON
    pipeline_name = next(iter(test_inputs)).split('.')[0]

    # Append "Test" in front of the pipeline name
    test_name = f"Test{pipeline_name}"

    # Update all keys in the json file to replace the pipeline name with the test name
    for key in list(test_inputs.keys()):
        new_key = key.replace(pipeline_name, test_name)
        test_inputs[new_key] = test_inputs.pop(key)

    # Add the truth_path and results_path to the test inputs JSON
    test_inputs[f"{test_name}.results_path"] = f"{results_path}/{sample_name}/"
    test_inputs[f"{test_name}.truth_path"] = f"{truth_path}/{sample_name}/"
    test_inputs[f"{test_name}.update_truth"] = update_truth
    test_inputs[f"{test_name}.commit_hash"] = commit_hash

    # Save the updated test inputs JSON
    output_name = f"updated_{sample_name}.json"
    with open(output_name, 'w') as file:
        json.dump(test_inputs, file, indent=4)

    #print out the contents of the updated json file
    with open(output_name, 'r') as file:
        print(file.read())

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
        "--commit_hash",
        required=True,
        help="Commit hash of the current pipeline run")


    args = parser.parse_args()
    # convert the update_truth flag to a boolean
    update_truth_bool = args.update_truth.lower() == "true"

    # Update the test inputs to work with the test wrapper WDL
    update_test_inputs(args.inputs_json, args.truth_path, args.results_path, update_truth_bool, args.commit_hash)


if __name__ == "__main__":
    main()