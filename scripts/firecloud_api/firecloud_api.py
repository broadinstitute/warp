import base64
import json
import requests
from datetime import datetime, timezone
from urllib.parse import quote
from google.auth.transport.requests import Request
from google.oauth2 import service_account
from google.auth import credentials
import argparse
import logging
import time
import sys

# Configure logging to display INFO level and above messages
logging.basicConfig(
    level=logging.INFO,  # This will show INFO and higher levels (INFO, WARNING, ERROR, CRITICAL)
    format='%(asctime)s - %(levelname)s - %(message)s'
)

class FirecloudAPI:
    def __init__(self, workspace_namespace, workspace_name, sa_json_b64, user, action, method_namespace, method_name):
        self.sa_json_b64 = sa_json_b64
        self.namespace = workspace_namespace
        self.workspace_name = workspace_name
        self.user = user  # Store the user email
        self.base_url = "https://api.firecloud.org/api"
        self.action = action
        self.method_namespace = method_namespace
        self.method_name = method_name

        # Setup credentials once during initialization
        scopes = ['profile', 'email', 'openid']
        decoded_sa = base64.b64decode(sa_json_b64).decode('utf-8')
        sa_credentials = service_account.Credentials.from_service_account_info(
            json.loads(decoded_sa),
            scopes=scopes
        )
        self.delegated_creds = sa_credentials.with_subject(user)

    def get_method_config_name(self, pipeline_name, branch_name, test_type):
        """
        Helper method to consistently generate method configuration names

        :param pipeline_name: Name of the pipeline
        :param branch_name: Name of the branch
        :param test_type: Type of test (Scientific or Plumbing)
        :return: Formatted method configuration name
        """
        return f"{pipeline_name}_{test_type}_{branch_name}"

    def build_auth_headers(self, token: str):
        if not self.delegated_creds.valid:
            logging.info("Refreshing credentials.")
            self.delegated_creds.refresh(Request())
        token = self.delegated_creds.token
        return {
            "content-type": "application/json",
            "Authorization": f"Bearer {token}",
        }

    def get_user_token(self, credentials: credentials):
        """
        Get test user's access token
        """
        # if token is expired or about to expire in 10 seconds, refresh and then use it
        if not credentials.valid:
            logging.info("Fetching user's new access token")
            credentials.refresh(Request())
            logging.info("Token refreshed.")
        else:
            expiry_timestamp = credentials.expiry.replace(tzinfo=timezone.utc).timestamp()
            now_timestamp = datetime.now(timezone.utc).timestamp()
            # if token is about to expire in 1 minute, refresh and then use it
            if expiry_timestamp - now_timestamp < 60:
                logging.info("Fetching user's new access token")
                credentials.refresh(Request())
                logging.info("Token refreshed.")

        return credentials.token

    def submit_job(self, submission_data_file):
        """
        Submits a job to Terra/Firecloud with retry logic for intermittent 500 errors.

        :param submission_data_file: The JSON data for the submission
        :return: The submission ID if successful, None otherwise
        """
        # Set up retry parameters
        max_retry_duration = 15 * 60  # 15 minutes in seconds
        start_time = time.time()
        retry_delay = 5  # Start with a 5-second delay between retries
        max_retry_delay = 30  # Maximum retry delay in seconds
        max_attempts = 10  # Maximum number of retry attempts

        attempts = 0
        while attempts < max_attempts:
            attempts += 1

            # Check if we've exceeded the maximum retry duration
            current_time = time.time()
            if current_time - start_time > max_retry_duration:
                logging.error(f"Exceeded maximum retry duration of {max_retry_duration/60} minutes.")
                return None

            try:
                token = self.get_user_token(self.delegated_creds)
                headers = self.build_auth_headers(token)
                url = f"{self.base_url}/workspaces/{self.namespace}/{quote(self.workspace_name)}/submissions"

                logging.info(f"Submitting job, attempt {attempts}/{max_attempts}")
                response = requests.post(url, json=submission_data_file, headers=headers)

                # Print status code and response body for debugging
                logging.info(f"Response status code for submitting job: {response.status_code}")

                # Handle different response codes
                if response.status_code == 201:  # Success
                    try:
                        # Parse the response as JSON
                        response_json = response.json()
                        logging.info(f"Response body: {response.text}")

                        # Extract the submissionId
                        submission_id = response_json.get("submissionId", None)
                        if submission_id:
                            logging.info(f"Submission ID extracted: {submission_id}")
                            return submission_id
                        else:
                            logging.error("Error: submissionId not found in the response.")
                            return None
                    except json.JSONDecodeError:
                        logging.error("Error: Failed to parse JSON response.")
                        logging.error(f"Response body: {response.text}")
                        # If we can't parse the JSON but got a 201, we might still want to retry
                        if attempts < max_attempts:
                            time.sleep(retry_delay)
                            retry_delay = min(retry_delay * 1.5, max_retry_delay)
                            continue
                        return None

                elif response.status_code == 500:  # Server error, retry
                    logging.warning(f"Received 500 error. Retrying in {retry_delay} seconds...")
                    logging.warning(f"Response body: {response.text}")
                    time.sleep(retry_delay)
                    # Implement exponential backoff with a cap
                    retry_delay = min(retry_delay * 1.5, max_retry_delay)
                    continue

                elif response.status_code >= 400 and response.status_code < 500:  # Client error
                    # For 4xx errors, only retry a few times as they might be temporary auth issues
                    logging.error(f"Client error (4xx): {response.status_code}")
                    logging.error(f"Response body: {response.text}")
                    if response.status_code == 401 or response.status_code == 403:
                        # Auth errors might be temporary, retry with token refresh
                        self.delegated_creds.refresh(Request())
                        if attempts < 3:  # Only retry auth errors a few times
                            time.sleep(retry_delay)
                            continue
                    return None

                else:  # Other error codes
                    logging.error(f"Failed to submit job. Status code: {response.status_code}")
                    logging.error(f"Response body: {response.text}")
                    if attempts < max_attempts:
                        time.sleep(retry_delay)
                        retry_delay = min(retry_delay * 1.5, max_retry_delay)
                        continue
                    return None

            except requests.exceptions.RequestException as e:
                # Handle network errors
                logging.warning(f"Network error occurred: {e}. Retrying in {retry_delay} seconds...")
                time.sleep(retry_delay)
                # Implement exponential backoff with a cap
                retry_delay = min(retry_delay * 1.5, max_retry_delay)
                continue

        logging.error(f"Failed to submit job after {max_attempts} attempts.")
        return None


    def create_new_method_config(self, branch_name, pipeline_name):
        """
        Creates a new method configuration in the workspace via Firecloud API.
        Includes a retry mechanism for 404 errors from Dockstore.

        :param branch_name: The branch name
        :param pipeline_name: The name of the pipeline
        :return: The name of the created method configuration or None if failed
        """
        # Create method config name with test type
        method_config_name = self.get_method_config_name(pipeline_name, branch_name, args.test_type)

        # Flag to track if we've already retried for a 404 error
        dockstore_404_retried = False

        # Function to create the payload
        def create_payload():
            return {
                "deleted": False,
                "inputs": {},
                "methodConfigVersion": 0,
                "methodRepoMethod": {
                    "methodUri": f"dockstore://github.com/broadinstitute/warp/{pipeline_name}/{branch_name}",
                    "sourceRepo": "dockstore",
                    "methodPath": f"github.com/broadinstitute/warp/{pipeline_name}",
                    "methodVersion": f"{branch_name}"
                },
                "name": method_config_name,
                "namespace": "warp-pipelines",
                "outputs": {},
                "prerequisites": {}
            }

        # Attempt to create the method configuration
        def attempt_creation():
            payload = create_payload()
            logging.info(f"Creating new method configuration: {json.dumps(payload, indent=2)}")

            # Construct the API endpoint URL for creating a new method configuration
            url = f"{self.base_url}/workspaces/{self.namespace}/{quote(self.workspace_name)}/method_configs/{self.namespace}/{method_config_name}"

            token = self.get_user_token(self.delegated_creds)
            headers = self.build_auth_headers(token)

            # Create the new method configuration in the workspace
            response = requests.put(url, headers=headers, json=payload)

            return response

        # First attempt
        response = attempt_creation()

        # Check if we got a 404 error (likely from Dockstore)
        if response.status_code == 404 and not dockstore_404_retried:
            error_message = response.text
            logging.warning(f"Received 404 error, possibly from Dockstore: {error_message}")
            logging.info(f"Waiting 5 minutes before retrying...")

            # Wait for 5 minutes (300 seconds)
            time.sleep(300)

            # Mark that we've retried for this error
            dockstore_404_retried = True

            # Retry the creation
            logging.info("Retrying method configuration creation after 5-minute wait")
            response = attempt_creation()

        # Final check if the method configuration was created successfully
        if response.status_code == 200:
            logging.info(f"Method configuration {method_config_name} created successfully.")
            return method_config_name
        else:
            logging.error(f"Failed to create method configuration. Status code: {response.status_code}")
            logging.error(f"Response body: {response.text}")
            raise Exception(f"Failed to create method configuration for {pipeline_name} on the branch {branch_name}")


    def upload_test_inputs(self, pipeline_name, test_inputs, branch_name, test_type):
        """
        Uploads test inputs to the workspace via Firecloud API.

        :param test_inputs: JSON data containing test inputs
        :param pipeline_name: The name of the pipeline
        :param branch_name: The name of the branch
        :param test_type: The type of test (Scientific or Plumbing)
        :return: True if successful, False otherwise
        """

        method_config_name = self.get_method_config_name(pipeline_name, branch_name, test_type)
        url = f"{self.base_url}/workspaces/{self.namespace}/{quote(self.workspace_name)}/method_configs/{self.namespace}/{method_config_name}"

        token = self.get_user_token(self.delegated_creds)
        headers = self.build_auth_headers(token)

        # get the current method configuration
        response = requests.get(url, headers=headers)

        if response.status_code == 404:
            logging.info(f"Method config {method_config_name} not found. Creating new config...")
            if not self.create_new_method_config(branch_name, pipeline_name):
                logging.error("Failed to create new method configuration.")
                return False
            response = requests.get(url, headers=headers)
            if response.status_code != 200:
                logging.error(f"Failed to get method configuration. Status code: {response.status_code}")
                return False

        config = response.json()
        print(f"Current method configuration: {json.dumps(config, indent=2)}")
        # update the config with the new inputs
        print(f"Opening test inputs file: {test_inputs}")
        with open(test_inputs, 'r') as file:
            inputs_json = json.load(file)
            print("Test inputs loaded successfully.")
            inputs_json = self.quote_values(inputs_json)
            print(f"here is test json after quote_values: {json.dumps(inputs_json, indent=2)}")
            config["inputs"] = inputs_json

        # Construct the methodUri with the branch name
        base_url = f"github.com/broadinstitute/warp/{pipeline_name}"
        method_uri = f"dockstore://{quote(base_url)}/{branch_name}"
        print(f"Updating methodUri with branch name: {method_uri}")
        config["methodRepoMethod"]["methodUri"] = method_uri

        print(f"Updating methodVersion with branch name: {branch_name}")
        config["methodRepoMethod"]["methodVersion"] = branch_name

        # We need to increment the methodConfigVersion by 1 every time we update the method configuration
        config["methodConfigVersion"] += 1  # Increment version number by  1
        print(f"Updated method configuration: {json.dumps(config, indent=2)}")


        # post the updated method config to the workspace
        response = requests.post(url, headers=headers, json=config)
        print(f"Response status code for uploading inputs: {response.status_code}")
        print(f"Response text: {response.text}")

        # Check if the test inputs were uploaded successfully
        if response.status_code == 200:
            print("Test inputs uploaded successfully.")
            return True
        else:
            print(f"Failed to upload test inputs. Status code: {response.status_code}")
            return False

    def poll_job_status(self, submission_id):
        """
        Polls the status of a submission until it is complete and returns a dictionary of workflow IDs and their statuses.
        Includes retry mechanism for handling intermittent 500 errors.

        :param submission_id: The ID of the submission to poll
        :return: Dictionary with workflow IDs as keys and their statuses as values
        """
        # Construct the API endpoint URL for polling submission status
        status_url = f"{self.base_url}/workspaces/{self.namespace}/{self.workspace_name}/submissions/{submission_id}"
        workflow_status_map = {}

        # Set up retry parameters
        start_time = time.time()
        retry_delay = 5  # Start with a 5-second delay between retries
        max_retry_delay = 30  # Maximum retry delay in seconds
        max_retry_duration = 15 * 60  # Maximum time to spend retrying server errors (15 minutes)

        # Continuously poll the status of the submission until completion
        while True:
            try:
                # Get the token and headers
                token = self.get_user_token(self.delegated_creds)
                headers = self.build_auth_headers(token)
                status_response = requests.get(status_url, headers=headers)

                # Check for 500 errors and retry if necessary
                if status_response.status_code in [500, 502, 503]:
                    elapsed_time = time.time() - start_time
                    logging.warning(f"Received {status_response.status_code} error. Retrying in {retry_delay} seconds...")
                    logging.warning(f"Response content: {status_response.text[:500]}")

                    # Check if we've exceeded the maximum retry duration
                    if elapsed_time > max_retry_duration:
                        logging.error(f"Exceeded maximum retry duration of {max_retry_duration/60} minutes for handling server errors.")
                        return {}

                    time.sleep(retry_delay)
                    # Implement exponential backoff with a cap
                    retry_delay = min(retry_delay * 1.5, max_retry_delay)
                    continue

                # Check if the response status code is successful (200)
                if status_response.status_code != 200:
                    logging.error(f"Error: Received status code {status_response.status_code}")
                    logging.info(f"Response content: {status_response.text}")
                    # For non-500 errors, wait and retry a few times
                    if time.time() - start_time <= 60:  # Only retry for the first minute for non-500 errors
                        logging.warning(f"Retrying in {retry_delay} seconds...")
                        time.sleep(retry_delay)
                        continue
                    return {}

                try:
                    # Parse the response as JSON
                    status_data = status_response.json()
                    # Reset retry delay after successful request
                    retry_delay = 5
                except json.JSONDecodeError:
                    logging.error("Error decoding JSON response.")
                    logging.info(f"Response content: {status_response.text}")
                    time.sleep(retry_delay)
                    continue

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
                    logging.info("Submission is done.")
                    break

                # Wait for 20 seconds before polling again
                time.sleep(20)

            except requests.exceptions.RequestException as e:
                # Handle network errors
                logging.warning(f"Network error occurred: {e}. Retrying in {retry_delay} seconds...")
                time.sleep(retry_delay)
                # Implement exponential backoff with a cap
                retry_delay = min(retry_delay * 1.5, max_retry_delay)

        return workflow_status_map

    def quote_values(self, inputs_json):
        """
        Format JSON values with proper handling of nested structures
        """
        def format_value(val):
            if isinstance(val, bool):
                return str(val).lower()
            elif isinstance(val, dict):
                return json.dumps(val, indent=2)
            elif isinstance(val, list):
                if all(isinstance(x, str) for x in val):
                    return json.dumps(val)
                return json.dumps([format_value(x) for x in val])
            elif isinstance(val, (int, float)):
                return str(val)
            elif val is None:
                return ""
            elif isinstance(val, str):
                if val.startswith("{") and val.endswith("}"):
                    try:
                        parsed = json.loads(val)
                        return json.dumps(parsed, indent=2)
                    except json.JSONDecodeError:
                        return f'"{val}"'
                return f'"{val}"'
            return f'"{str(val)}"'

        return {key: format_value(value) for key, value in inputs_json.items()}

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
            logging.error(f"Failed to retrieve workflow outputs. Status code: {response.status_code}")
            return None, None

    def delete_method_config(self, method_config_name):
        """
        Deletes a method configuration from the workspace.

        :param method_config_name: The name of the method configuration to delete
        :return: True if deletion is successful, False otherwise
        """
        url = f"{self.base_url}/workspaces/{self.namespace}/{quote(self.workspace_name)}/method_configs/{self.namespace}/{method_config_name}"

        token = self.get_user_token(self.delegated_creds)
        headers = self.build_auth_headers(token)

        # Send a DELETE request to delete the method configuration
        response = requests.delete(url, headers=headers)

        if response.status_code == 204:
            logging.info(f"Method configuration {method_config_name} deleted successfully.")
            print("True")
            return True
        else:
            logging.error(f"Failed to delete method configuration {method_config_name}. Status code: {response.status_code}")
            logging.error(f"Response body: {response.text}")
            return False

    def get_active_submissions(self, method_config_name=None):
        """
        Get all active workflow submissions for the workspace.
        Optionally filter by method configuration name.
        """
        url = f"{self.base_url}/workspaces/{self.namespace}/{quote(self.workspace_name)}/submissions"
        token = self.get_user_token(self.delegated_creds)
        headers = self.build_auth_headers(token)

        response = requests.get(url, headers=headers)

        if response.status_code != 200:
            logging.error(f"Failed to get submissions. Status code: {response.status_code}")
            logging.error(f"Response body: {response.text}")
            return []

        submissions = response.json()
        active_submissions = []

        for submission in submissions:
            # Check if submission is active (not Done, Aborted, or Failed)
            if submission['status'] in ['Submitted', 'Running', 'Queued']:
                config_name = submission.get('methodConfigurationName', '')
                if config_name.startswith(method_config_name):
                    active_submissions.append(submission)

        return active_submissions

    def cancel_submission(self, submission_id):
        """
        Cancel a specific workflow submission.
        """
        url = f"{self.base_url}/workspaces/{self.namespace}/{quote(self.workspace_name)}/submissions/{submission_id}"
        token = self.get_user_token(self.delegated_creds)
        headers = self.build_auth_headers(token)

        response = requests.delete(url, headers=headers)

        if response.status_code not in [204]:
            logging.error(f"Failed to cancel submission {submission_id}. Status code: {response.status_code}")
            logging.error(f"Response body: {response.text}")
            return False

        logging.info(f"Successfully cancelled submission {submission_id}")
        return True

    def cancel_old_submissions(self, pipeline_name, branch_name):
        """
        Cancel all active submissions for a pipeline's method configuration.
        Returns the number of cancelled submissions.
        """
        method_config_name = self.get_method_config_name(pipeline_name, branch_name, args.test_type)
        active_submissions = self.get_active_submissions(method_config_name)
        cancelled_count = 0

        for submission in active_submissions:
            if self.cancel_submission(submission['submissionId']):
                cancelled_count += 1
                logging.info(f"Cancelled submission {submission['submissionId']}")

        return cancelled_count


    def main(self):
        logging.info("Starting process based on action.")

        if self.action == "submit_job":
            submission_id = self.submit_job()
            logging.info(f"Job submission complete with ID: {submission_id}")
        elif self.action == "create_new_method_config":
            if not args.pipeline_name or not args.branch_name:
                parser.error("Arguments --pipeline_name and --branch_name are required for 'create_new_method_config'")
            method_config_name = self.create_new_method_config(args.branch_name, args.pipeline_name)
            print(method_config_name)
            if method_config_name:
                logging.info(f"Method configuration created with name: {method_config_name}")
            else:
                logging.error("Failed to create method configuration.")
        elif self.action == "delete_method_config":
            if not args.method_config_name:
                if not all([args.pipeline_name, args.branch_name]):
                    parser.error("Either --method_config_name or both --pipeline_name and --branch_name are required")
                method_config_name = self.get_method_config_name(args.pipeline_name, args.branch_name, args.test_type)
            else:
                method_config_name = args.method_config_name
            result = self.delete_method_config(method_config_name)
            print(str(result).lower())
        elif self.action == "upload_test_inputs":
            success = self.upload_test_inputs(self.pipeline_name, self.test_input_file, self.branch_name, self.test_type)
            if success:
                logging.info("Test inputs uploaded successfully.")
            else:
                logging.error("Failed to upload test inputs.")
        elif self.action == "poll_job_status":
            status = self.poll_job_status()
            logging.info(f"Final job status: {status}")
        elif self.action == "create_new_method_config":
            method_config_name = self.create_new_method_config(self.branch_name, self.pipeline_name)
            if method_config_name:
                logging.info("Method configuration created successfully.")
            else:
                logging.error("Failed to create method configuration.")
        elif self.action == "delete_method_config":
            if not args.method_config_name:
                parser.error("Argument --method_config_name is required for 'delete_method_config'")
            else:
                # Delete the method configuration
                result = self.delete_method_config(args.method_config_name)
                if result:
                    logging.info("Method configuration deleted successfully.")
                else:
                    logging.error("Failed to delete method configuration.")
        elif self.action == "get_workflow_outputs":
            if not args.submission_id or not args.workflow_id or not args.pipeline_name:
                parser.error("Arguments --submission_id, --workflow_id, and --pipeline_name are required for 'get_workflow_outputs'")
            # Fetch workflow outputs
            outputs, output_values = self.get_workflow_outputs(args.submission_id, args.workflow_id, args.pipeline_name)
            if outputs:
                logging.info(f"Workflow outputs: {json.dumps(outputs, indent=2)}")
                logging.info(f"Output values: {output_values}")
            else:
                logging.error("Failed to retrieve workflow outputs.")
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
    parser.add_argument("--method_namespace", help="Method namespace")
    parser.add_argument("--method_name", help="Method name")
    parser.add_argument('--submission_data_file', help='Path to submission data JSON file (required for submit)')
    parser.add_argument('--submission_id', help='Submission ID (required for poll_job_status)')
    parser.add_argument('--workflow_id', help='Workflow ID (required for get_workflow_outputs)')
    parser.add_argument("--source", help="Source GCS path for gsutil copy")
    parser.add_argument("--destination", help="Destination GCS path for gsutil copy")
    parser.add_argument("--method_config_name", help="Name of the method configuration to delete")
    parser.add_argument("--test_type", help="Test type (Scientific or Plumbing)")
    parser.add_argument("action", choices=["submit_job", "upload_test_inputs", "poll_job_status", "get_workflow_outputs", "create_new_method_config", "delete_method_config", "cancel_old_submissions"],
                        help="Action to perform: 'submit_job', 'upload_test_inputs', 'poll_job_status', 'get_workflow_outputs',  'create_new_method_config', or 'delete_method_config'")

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


    if args.action == "upload_test_inputs":
        # Check for required arguments for upload_test_inputs action
        if not args.pipeline_name or not args.test_input_file or not args.branch_name:
            parser.error("Arguments --pipeline_name, --test_input_file, and --branch_name are required for 'upload_test_inputs'")
        # Call the function to upload test inputs
        api.upload_test_inputs(args.pipeline_name, args.test_input_file, args.branch_name, args.test_type)

    elif args.action == "submit_job":
        # Check for required argument for submit_job action
        if not args.submission_data_file:
            parser.error("Argument --submission_data_file is required for 'submit_job'")
        # Load the submission data from the provided file
        else:
            with open(args.submission_data_file, 'r') as file:
                submission_data = json.load(file)
            # Submit the job with the loaded submission data
            submission_id = api.submit_job(submission_data)
            if submission_id:
              print(submission_id)
              logging.info("Submission successful.")
            else:
              logging.error("Submission failed.")
              sys.exit(1)

    elif args.action == "poll_job_status":
        if not args.submission_id:
            parser.error("Argument --submission_id is required for 'poll_job_status'")
        else:
            # Poll the job status with the provided submission ID
            workflow_status_map = api.poll_job_status(args.submission_id)

            # Convert the dictionary to a JSON string and print it
            if workflow_status_map:
                print(json.dumps(workflow_status_map))  # Output the dictionary as a JSON string for bash parsing
            else:
                print("No workflows found or an error occurred.")
 
    elif args.action == "create_new_method_config":
        # Check for required arguments for create_new_method_config action
        if not args.pipeline_name or not args.branch_name:
            parser.error("Arguments --pipeline_name and --branch_name are required for 'create_new_method_config'")
        # Call the function to create a new method configuration
        method_config_name = api.create_new_method_config(args.branch_name, args.pipeline_name)
        print(method_config_name)
        if method_config_name:
            logging.info(f"Method configuration created with name: {method_config_name}")
        else:
            logging.error("Failed to create method configuration.")
    elif args.action == "delete_method_config":
        if not args.method_config_name:
            parser.error("Argument --method_config_name is required for 'delete_method_config'")
        else:
            # Delete the method configuration
            result = api.delete_method_config(args.method_config_name)
            if result:
                logging.info("Method configuration deleted successfully.")
            else:
                logging.error("Failed to delete method configuration.")
    elif args.action == "cancel_old_submissions":
        if not all([args.pipeline_name, args.branch_name]):
            parser.error("Arguments --pipeline_name and --branch_name are required for 'cancel_old_submissions'")

        # Cancel old submissions
        cancelled_count = api.cancel_old_submissions(
            args.pipeline_name,
            args.branch_name
        )
        print(f"Cancelled {cancelled_count} old submissions")




