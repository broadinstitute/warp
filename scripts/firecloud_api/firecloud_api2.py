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
import os


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
        # Construct the API endpoint URL for creating a new submission
        url = f"{self.base_url}/workspaces/{self.namespace}/{quote(self.workspace_name)}/submissions"
        response = requests.post(url, json=submission_data_file, headers=self.headers)
        print(f"Response received. Status code: {response.status_code}")


        # Check if the submission was created successfully
        if response.status_code == 201:
            submission_id = response.json().get('submissionId')
            print(f"Job submitted successfully. Submission ID is: {submission_id}")
            os.environ['SUBMISSION_ID'] = submission_id
            print(f"Job submitted successfully. Submission ID again is: {submission_id}")
            return submission_id
        else:
            print(f"Failed to submit job. Status code: {response.status_code}")
            print(f"Response content: {response.text}")
            return None


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

    def poll_job_status(self, submission_id, polling_interval=30, timeout=3600):
        """
        Polls the status of a job submission until it reaches a terminal state.
        :param submission_id: The ID of the job submission to monitor.
        :param polling_interval: Time (in seconds) between status checks.
        :param timeout: Maximum time (in seconds) to wait before giving up.
        :return: The final status of the job.
        """
        uri = f"{self.base_url}/workspaces/{self.namespace}/{quote(self.workspace_name)}/submissions/{submission_id}"
        start_time = datetime.now()

        while (datetime.now() - start_time).total_seconds() < timeout:
            response = requests.get(uri, headers=self.headers)
            if response.status_code != 200:
                print(f"Failed to fetch submission status. Status code: {response.status_code}. Response: {response.text}")
                raise Exception("Failed to fetch job status.")

            status = response.json().get("status")
            print(f"Current status for submission {submission_id}: {status}")

            if status in ["Done", "Failed", "Aborted"]:
                print(f"Job {submission_id} reached terminal status: {status}")
                return status

            sleep(polling_interval)

        raise TimeoutError(f"Polling timed out after {timeout} seconds for submission {submission_id}.")

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
        elif self.action == "poll_job_status":
            status = self.poll_job_status()
            logging.info(f"Final job status: {status}")
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
        choices=["submit_job", "upload_test_inputs", "poll_job_status"],
        help="Action to perform: 'submit_job' or 'upload_test_inputs' or 'poll_job_status'"
    )
    parser.add_argument("--method_namespace", help="Method namespace")
    parser.add_argument("--method_name", help="Method name")
    parser.add_argument('--submission_data_file', help='Path to submission data JSON file (required for submit)')
    parser.add_argument('--submission_id', help='Submission ID (required for poll_job_status)')
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
        # Check for required arguments for upload_test_inputs action
        if not args.pipeline_name or not args.test_input_file or not args.branch_name:
            parser.error("Arguments --pipeline_name, --test_input_file, and --branch_name are required for 'upload_test_inputs'")
        # Call the function to upload test inputs
        api.upload_test_inputs(args.pipeline_name, args.test_input_file, args.branch_name)

    elif args.action == "submit_job":
        # Check for required argument for submit_job action
        if not args.submission_data_file:
            parser.error("Argument --submission_data_file is required for 'submit_job'")
        # Load the submission data from the provided file
        with open(args.submission_data_file, 'r') as file:
            submission_data = json.load(file)
        # Submit the job with the loaded submission data
        api.submit_job(submission_data)

    elif args.action == "pollpoll_job_status":
        # Check for required argument for poll_job_status action
        if not args.submission_id:
            parser.error("Argument --submission_id is required for 'poll_job_status'")
        # Poll the job status with the provided submission ID
        api.poll_job_status(args.submission_id)




    #api.main()
