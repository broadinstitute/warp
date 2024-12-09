"""
firecloud_api.py
Author: Kevin Palis <palis@broadinstitute.org>

This module provides an object-oriented interface for interacting with the Firecloud REST API.
It includes functionalities to submit workflows, retrieve workflow outputs, and monitor
workflow statuses.

Classes:
    - FirecloudAPI: A class to handle Firecloud API interactions.

Usage:
    Initialize the FirecloudAPI class with API token, namespace, and workspace details, and
    call its methods to interact with the Firecloud service.
"""

import requests
import time
import json
import sys
import subprocess
from urllib.parse import quote
import base64
import os
from google.auth.transport.requests import Request
from google.auth import credentials
from google.oauth2 import service_account

sa_json_b64 = os.environ.get("SA_JSON_B64")
print(f"Service account JSON: {sa_json_b64}")

def is_base64_encoded(data):
    """Checks if a string is a valid Base64-encoded string."""
    try:
        # Try to decode and re-encode the string to check its validity
        base64.b64decode(data, validate=True)
        return True
    except (base64.binascii.Error, TypeError):
        return False
# Retrieve SA_JSON_B64 from the environment variable
sa_json_b64 = os.environ.get("SA_JSON_B64")
print(f"Service account JSON (Base64): {sa_json_b64}")

# Check if the string is Base64-encoded
if sa_json_b64 and is_base64_encoded(sa_json_b64):
    try:
        # Decode the Base64 string and parse as JSON
        decoded_sa = base64.b64decode(sa_json_b64).decode('utf-8')
        sa_credentials = json.loads(decoded_sa)  # Assuming the decoded content is a valid JSON
        print("Service account credentials loaded successfully.")
    except (json.JSONDecodeError, UnicodeDecodeError) as e:
        print(f"Failed to parse decoded JSON: {e}")
else:
    print("The SA_JSON_B64 is not a valid Base64-encoded string.")


try:
    scopes = ['profile', 'email', 'openid']
    decoded_sa = base64.b64decode(sa_json_b64).decode('utf-8')
    sa_credentials = service_account.Credentials.from_service_account_info(json.loads(decoded_sa), scopes=scopes)
    print("Service account credentials loaded successfully.")
except Exception as e:
    print(f"Failed to load service account credentials: {e}")
    sa_credentials = None  # Set a fallback or exit as appropriate

class FirecloudAPI:
    def __init__(self, token, namespace, workspace_name):
        """
        Initializes the FirecloudAPI object with authentication and workspace details.

        :param token: API access token
        :param namespace: Workspace namespace
        :param workspace_name: Workspace name
        """
        self.token = token
        self.namespace = namespace
        self.workspace_name = workspace_name
        self.base_url = "https://api.firecloud.org/api"
        self.headers = {
            'accept': '*/*',
            'Authorization': f'Bearer {self.token}',
        }

    def get_workflow_outputs(self, submission_id, workflow_id, pipeline_name):
        """
        Fetches workflow outputs from the Firecloud API.

        :param submission_id: The ID of the submission
        :param workflow_id: The ID of the workflow
        :param pipeline_name: The name of the pipeline whose outputs are required
        :return: Outputs dictionary and a list of output values
        """
        # Construct the API endpoint URL for fetching workflow outputs
        url = f"{self.base_url}/workspaces/{self.namespace}/{self.workspace_name}/submissions/{submission_id}/workflows/{workflow_id}/outputs"
        response = requests.get(url, headers=self.headers)

        # Check if the API request was successful
        if response.status_code == 200:
            json_response = response.json()
            # Extract outputs for the specified pipeline name
            outputs = json_response.get('tasks', {}).get(pipeline_name, {}).get('outputs', {})
            output_values = list(outputs.values())
            return outputs, output_values
        else:
            print(f"Failed to retrieve workflow outputs. Status code: {response.status_code}")
            return None, None

    #def refresh_token(self):
    #    """
    #    Refreshes the API token using gcloud's application default credentials.
    #    :return: The new token as a string
    #    """
    #    try:
    #        # Execute the gcloud command to get the new access token
    #        result = subprocess.run(
    #            ["gcloud", "auth", "application-default", "print-access-token"],
    #            capture_output=True, text=True, check=True
    #        )
    #        return result.stdout.strip()  # Return the new token
    #    except subprocess.CalledProcessError as e:
    #        print(f"Error refreshing token: {e.stderr}", file=sys.stderr)
    #        return None

    def create_submission(self, submission_data):
        """
        Submits a workflow to the Firecloud API.

        :param submission_data: JSON data containing submission details
        :return: Submission ID if successful, None otherwise
        """
        # Construct the API endpoint URL for creating a new submission
        url = f"{self.base_url}/workspaces/{self.namespace}/{self.workspace_name}/submissions"
        response = requests.post(url, headers=self.headers, json=submission_data)

        # Check if the submission was created successfully
        if response.status_code == 201:
            submission_id = response.json().get('submissionId')
            #print(f"Submission created with ID: {submission_id}")
            return submission_id
        else:
            print(f"Failed to create submission. Status code: {response.status_code}")
            print(f"Response content: {response.text}")
            return None


    def poll_submission_status(self, submission_id):
        """
        Polls the status of a submission until it is complete and returns a dictionary of workflow IDs and their statuses.

        :param submission_id: The ID of the submission to poll
        :return: Dictionary with workflow IDs as keys and their statuses as values
        """
        # Construct the API endpoint URL for polling submission status
        status_url = f"{self.base_url}/workspaces/{self.namespace}/{self.workspace_name}/submissions/{submission_id}"
        workflow_status_map = {}

        # Continuously poll the status of the submission until completion
        while True:
            status_response = requests.get(status_url, headers=self.headers)

            # Check if the response status code is successful (200)
            if status_response.status_code != 200:
                print(f"Error: Received status code {status_response.status_code}", file=sys.stderr)
                print(f"Response content: {status_response.text}", file=sys.stderr)
                return {}

            try:
                # Parse the response as JSON
                status_data = status_response.json()
            except json.JSONDecodeError:
                print("Error decoding JSON response.", file=sys.stderr)
                print(f"Response content: {status_response.text}", file=sys.stderr)
                return {}

            # Retrieve workflows and their statuses
            workflows = status_data.get("workflows", [])
            for workflow in workflows:
                workflow_id = workflow.get("workflowId")
                workflow_status = workflow.get("status")
                if workflow_id and workflow_status:
                    workflow_status_map[workflow_id] = workflow_status

            # Check if the submission is complete
            submission_status = status_data.get("status", "")
            if submission_status == "Done":
                break

            # Wait for 60 seconds before polling again
            time.sleep(20)

        return workflow_status_map

    def quote_values(self, data):
        """
        Recursively quotes values in a dictionary or list to match Firecloud API format.
        """
        if isinstance(data, dict):
            return {key: self.quote_values(value) for key, value in data.items()}
        elif isinstance(data, list):
            return [self.quote_values(item) for item in data]
        elif isinstance(data, (str, int, float, bool)):
            return f"\"{data}\""
        else:
            return data  # Return as-is if it's not a string, int, float, or boolean


    def upload_test_inputs(self, pipeline_name, test_inputs, branch_name):
        """
        Uploads test inputs to the workspace via Firecloud API.

        :param test_inputs: JSON data containing test inputs
        :return: True if successful, False otherwise
        """
        # Construct the API endpoint URL for the method configuration
        # properly encode the space in WARP Tests as %20 using from urllib.parse import quote
        url = f"{self.base_url}/workspaces/{self.namespace}/{quote(self.workspace_name)}/method_configs/{self.namespace}/{pipeline_name}"

        print(url)
        # get the current method configuration
        response = requests.get(url, headers=self.headers)
        config = response.json()
        print(f"Current method configuration: {json.dumps(config, indent=2)}")
        # update the config with the new inputs
        print(f"Opening test inputs file: {test_inputs}")
        with open(test_inputs, 'r') as file:
            inputs_json = json.load(file)
            print("Test inputs loaded successfully.")
            inputs_json = self.quote_values(inputs_json)
            config["inputs"] = inputs_json

        # Construct the methodUri with the branch name
        base_url = "github.com/broadinstitute/warp/{pipeline_name}"
        method_uri = f"dockstore://{quote(base_url)}/{branch_name}"
        print(f"Updating methodUri with branch name: {method_uri}")
        config["methodRepoMethod"]["methodUri"] = method_uri

        print(f"Updating methodVersion with branch name: {branch_name}")
        config["methodRepoMethod"]["methodVersion"] = branch_name

        # We need to increment the methodConfigVersion by 1 every time we update the method configuration
        config["methodConfigVersion"] += 1  # Increment version number by  1
        print(f"Updated method configuration: {json.dumps(config, indent=2)}")


       # post the updated method config to the workspace
        response = requests.post(url, headers=self.headers, json=config)
        print(f"Response status code: {response.status_code}")
        print(f"Response text: {response.text}")

        # Check if the test inputs were uploaded successfully
        if response.status_code == 200:
            print("Test inputs uploaded successfully.")
            return True
        else:
            print(f"Failed to upload test inputs. Status code: {response.status_code}")
            return False


# Bash Script Interaction
if __name__ == "__main__":
    import argparse

    # Set up command-line argument parsing
    parser = argparse.ArgumentParser(description='Interact with Firecloud API.')
    parser.add_argument('--token', required=True, help='API access token')
    parser.add_argument('--namespace', required=True, help='Workspace namespace')
    parser.add_argument('--workspace', required=True, help='Workspace name')
    parser.add_argument('--action', required=True, choices=['get_outputs', 'submit', 'poll_status', 'upload_test_inputs'], help='Action to perform')
    parser.add_argument('--submission_id', help='Submission ID (required for get_outputs and poll_status)')
    parser.add_argument('--workflow_id', help='Workflow ID (required for get_outputs)')
    parser.add_argument('--pipeline_name', help='Pipeline name (required for get_outputs)')
    parser.add_argument('--submission_data_file', help='Path to submission data JSON file (required for submit)')
    parser.add_argument('--test_input_file', help='Path to test inputs JSON file (required for upload_test_inputs)')
    parser.add_argument('--branch_name', help='Branch name for the method configuration (required for upload)')

    args = parser.parse_args()

    # Initialize the FirecloudAPI instance with provided arguments
    firecloud_api = FirecloudAPI(args.token, args.namespace, args.workspace)

    # Perform actions based on the specified action argument
    if args.action == 'get_outputs':
        if not all([args.submission_id, args.workflow_id, args.pipeline_name]):
            print("For 'get_outputs', --submission_id, --workflow_id, and --pipeline_name are required.")
        else:
            outputs, output_values = firecloud_api.get_workflow_outputs(args.submission_id, args.workflow_id, args.pipeline_name)
            #print(outputs)
            # Convert the dictionary, outputs, to a JSON string and print it
            if outputs:
                print(json.dumps(outputs)) # Output the dictionary as a JSON string for bash parsing
            else:
                print("No outputs found or an error occurred.", file=sys.stderr)

    elif args.action == 'submit':
        if not args.submission_data_file:
            print("For 'submit', --submission_data_file is required.")
        else:
            # Load submission data from the specified JSON file
            with open(args.submission_data_file, 'r') as file:
                submission_data = json.load(file)
            submission_id = firecloud_api.create_submission(submission_data)
            print(submission_id)

    elif args.action == 'poll_status':
      if not args.submission_id:
          print("For 'poll_status', --submission_id is required.", file=sys.stderr)
      else:
          workflow_status_map = firecloud_api.poll_submission_status(args.submission_id)

          # Convert the dictionary to a JSON string and print it
          if workflow_status_map:
              print(json.dumps(workflow_status_map))  # Output the dictionary as a JSON string for bash parsing
          else:
              print("No workflows found or an error occurred.", file=sys.stderr)

    elif args.action == 'upload_test_inputs':
        if not all([args.pipeline_name, args.test_input_file, args.branch_name]):
            print("For 'upload_test_inputs', --pipeline_name, --test_input_file and --branch_name are required.", file=sys.stderr)
        else:
            success = firecloud_api.upload_test_inputs(args.pipeline_name, args.test_input_file, args.branch_name)
            print(success)