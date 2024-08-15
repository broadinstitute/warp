import requests
import argparse

def get_workflow_outputs(token, namespace, workspace_name, submission_id, workflow_id, pipeline_name):
    # API endpoint to get the workflow outputs
    url = f"https://api.firecloud.org/api/workspaces/{namespace}/{workspace_name}/submissions/{submission_id}/workflows/{workflow_id}/outputs"
    print(f"Requesting URL: {url}")

    # Headers including the authorization token
    headers = {
        'accept': '*/*',
        'Authorization': f'Bearer {token}',
    }

    # Make the GET request
    response = requests.get(url, headers=headers)

    # Check if the request was successful
    if response.status_code == 200:
        json_response = response.json()  # parse the JSON response
        # extract the outputs section using the task name
        outputs = json_response.get('tasks', {}).get(pipeline_name, {}).get('outputs', {})

        # Turn the outputs dictionary into a list of values
        output_values = list(outputs.values())

        return outputs, output_values
    else:
        print(f"Failed to retrieve workflow outputs. Status code: {response.status_code}")
        return None, None

if __name__ == "__main__":
    # Define the command-line arguments
    parser = argparse.ArgumentParser(description='Fetch workflow outputs from the API.')
    parser.add_argument('--token', required=True, help='Authentication token')
    parser.add_argument('--namespace', required=True, help='Workspace namespace')
    parser.add_argument('--workspace', required=True, help='Workspace name')
    parser.add_argument('--submission_id', required=True, help='Submission ID')
    parser.add_argument('--workflow_id', required=True, help='Workflow ID')
    parser.add_argument('--pipeline_name', required=True, help='Name of the pipeline')

    # Parse the arguments
    args = parser.parse_args()

    # Call the function with the parsed arguments
    outputs, output_values = get_workflow_outputs(args.token, args.namespace, args.workspace, args.submission_id, args.workflow_id, args.pipeline_name)

    if outputs:
        print("Outputs:")
        print(outputs)

        print("\nOutput Values:")
        print(output_values)