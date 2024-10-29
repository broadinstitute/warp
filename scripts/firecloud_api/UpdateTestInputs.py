import argparse
import json


def update_test_inputs(test_inputs, truth_path, test_path):
    # Update the test inputs JSON to work with the test wrapper WDL
    # The test wrapper WDL runs the pipeline WDL and verifies the results
    # The test wrapper WDL requires the following inputs:
    # - truth_path: The path to the truth data
    # - test_path: The path to the test data
    # The test inputs JSON will be updated to include these inputs

    # Get the pipeline name from the test inputs JSON
    pipeline_name = next(iter(test_inputs)).split('.')[0]

    #Append "Test" in front of the pipeline name
    test_name = f"Test{pipeline_name}"

    # Update all keys in the json file to replace the pipeline name with the test name
    for key in list(test_inputs.keys()):
        new_key = key.replace(pipeline_name, test_name)
        test_inputs[new_key] = test_inputs.pop(key)

    # Add the truth_path and test_path to the test inputs JSON
    test_inputs[f"{test_name}.truth_path"] = truth_path
    test_inputs[f"{test_name}.test_path"] = test_path

    # Save the updated test inputs JSON
    with open("updated_test_inputs.json", 'w') as file:
        json.dump(test_inputs, file, indent=4)

    print("Test inputs JSON updated with truth_path and test_path")


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
        "--test_path",
        dest="test_path",
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
        default=False,
        required=False,
        choices=[True, False],
        help="Boolean flag to update the truth data. If True, the truth data will be updated with the test data. ",
    )

    args = parser.parse_args()

    # Load the JSON file containing the test inputs
    with open(args.inputs_json, 'r') as file:
        test_inputs = json.load(file)

    # Update the test inputs to work with the test wrapper WDL
    update_test_inputs(test_inputs, args.truth_path, args.test_path)

if __name__ == "__main__":
    main()