import base64
import json
import logging
import requests
import traceback
from time import sleep
from datetime import datetime, timezone
from urllib.parse import quote
from google.auth.transport.requests import Request
from google.oauth2 import service_account
import argparse
import sys


class FirecloudAPI:
    def __init__(self, workspace_namespace, workspace_name, sa_json_b64, user, action, method_namespace, method_name):
        self.sa_json_b64 = sa_json_b64
        self.namespace = workspace_namespace
        self.workspace_name = workspace_name
        self.user = user  # Store the user email
        self.base_url = "https://api.firecloud.org/api"
        self.headers = self._build_auth_headers()
        self.action = action
        self.method_namespace = method_namespace
        self.method_name = method_name


    def _build_auth_headers(self):
        scopes = ["profile", "email", "openid"]
        #sa_credentials = service_account.Credentials.from_service_account_info(
        #    json.loads(base64.b64decode(self.sa_json_b64).decode("utf-8")), scopes=scopes
        #)
        #TODO - Fix this, probably needs to be encoded in base64
        sa_credentials = service_account.Credentials.from_service_account_info(
            json.loads(self.sa_json_b64), scopes=scopes
        )
        delegated_credentials = sa_credentials.with_subject(self.user)
        token = self._get_user_token(delegated_credentials)
        return {
            "content-type": "application/json",
            "Authorization": f"Bearer {token}",
        }

    def _get_user_token(self, credentials):
        if not credentials.valid or (credentials.expiry and (credentials.expiry - datetime.now(timezone.utc)).total_seconds() < 60):
            logging.info("Refreshing user access token.")
            credentials.refresh(Request())
        return credentials.token

    def submit_job(self, submission_data_file):
        #logging.info(f"Submitting job for method {self.method_namespace}/{self.method_name} in workspace {self.namespace}/{self.workspace_name}.")
        uri = f"{self.base_url}/workspaces/{self.namespace}/{self.workspace_name}/submissions"
        response = requests.post(uri, json=submission_data_file, headers=self.headers)

        # Check if the submission was created successfully
        if response.status_code != 201:
            logging.error(f"Failed to submit job. Status code: {response.status_code}. Response: {response.text}")
            raise Exception("Submission failed.")
        submission_id = response.json().get("submissionId")
        logging.info(f"Job submitted successfully. Submission ID: {submission_id}")
        return submission_id

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

    @staticmethod
    def quote_values(inputs_json):
        return {key: f'"{value}"' for key, value in inputs_json.items()}

    def main(self):
        logging.info("Starting process based on action.")

        if self.action == "submit_job":
            submission_id = self.submit_job()
            logging.info(f"Job submission complete with ID: {submission_id}")
        elif self.action == "upload_test_inputs":
            success = self.upload_test_inputs(self.pipeline_name, self.test_input_file, self.branch_name)
            if success:
                logging.info("Test inputs uploaded successfully.")
            else:
                logging.error("Failed to upload test inputs.")
        else:
            logging.error(f"Unknown action: {self.action}")


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--sa-json-b64", required=True, help="Base64 encoded service account JSON")
    parser.add_argument("--user", required=True, help="User email for impersonation")
    parser.add_argument("--workspace-namespace", required=True, help="Namespace of the workspace.")
    parser.add_argument("--workspace-name", required=True, help="Name of the workspace.")
    parser.add_argument("--pipeline_name", help="Pipeline name (required for 'upload_test_inputs')")
    parser.add_argument("--test_input_file", help="Path to test input file (required for 'upload_test_inputs')")
    parser.add_argument("--branch_name", help="Branch name for the method repository (required for 'upload_test_inputs')")
    parser.add_argument(
        "action",
        choices=["submit_job", "upload_test_inputs"],
        help="Action to perform: 'submit_job' or 'upload_test_inputs'"
    )
    parser.add_argument("--method_namespace", help="Method namespace")
    parser.add_argument("--method_name", help="Method name")
    parser.add_argument('--submission_data_file', help='Path to submission data JSON file (required for submit)')
    args = parser.parse_args()

    # Pass action to the FirecloudAPI constructor
    api = FirecloudAPI(
        sa_json_b64=args.sa_json_b64,
        user=args.user,
        workspace_namespace=args.workspace_namespace,
        workspace_name=args.workspace_name,
        action=args.action,
        method_namespace=args.method_namespace,
        method_name=args.method_name
    )

    # Perform the selected action
    if args.action == "upload_test_inputs":
        if not args.pipeline_name or not args.test_input_file or not args.branch_name:
            parser.error("Arguments --pipeline_name, --test_input_file, and --branch_name are required for 'upload_test_inputs'")
        api.upload_test_inputs(args.pipeline_name, args.test_input_file, args.branch_name)
    elif args.action == "submit_job":
        api.submit_job()

    #api.main()
