import base64
import json
import logging
import os
import requests
import traceback
from time import sleep
from datetime import datetime, timezone
from google.auth.transport.requests import Request
from google.auth import credentials
from google.oauth2 import service_account

# Configuration and environment variables
USER = os.getenv("pdt-tester@warp-pipeline-dev.iam.gserviceaccount.com")
WORKSPACE_NAMESPACE = os.getenv("warp-pipelines")
WORKSPACE_NAME = os.getenv("WARP Tests")
METHOD_NAMESPACE = os.getenv("warp-")
METHOD_NAME = os.getenv("METHOD_NAME")
ENTITY_TYPE = os.getenv("ENTITY_TYPE")
ENTITY_ID = os.getenv("ENTITY_ID")
sa_json_b64 = os.environ.get("SA_JSON_B64")

# Configure logging
LOG_FORMAT = "%(asctime)s %(levelname)-8s %(message)s"
LOG_LEVEL = "INFO"
logging.basicConfig(
    format=LOG_FORMAT,
    level=getattr(logging, LOG_LEVEL),
    datefmt="%Y-%m-%d %H:%M:%S",
)
class FireCloudAPI:
    def __init__(self, sa_json_b64, user, workspace_namespace, workspace_name, method_namespace, method_name, entity_type, entity_id):
        """
        Initialize the FireCloudJobManager with configuration.
        """
        self.user = user
        self.workspace_namespace = workspace_namespace
        self.workspace_name = workspace_name
        self.method_namespace = method_namespace
        self.method_name = method_name
        self.entity_type = entity_type
        self.entity_id = entity_id
        self.credentials = self._load_credentials(sa_json_b64)

        # Configure logging
        LOG_FORMAT = "%(asctime)s %(levelname)-8s %(message)s"
        logging.basicConfig(format=LOG_FORMAT, level=logging.INFO, datefmt="%Y-%m-%d %H:%M:%S")

    def _load_credentials(self, sa_json_b64):
        """
        Load the service account credentials.
        """
        scopes = ["profile", "email", "openid"]
        decoded_sa = json.loads(sa_json_b64)
        sa_credentials = service_account.Credentials.from_service_account_info(decoded_sa, scopes=scopes)
        return sa_credentials.with_subject(self.user)

    def get_user_token(self):
        """
        Obtain and refresh the user access token.
        """
        if not self.credentials.valid or (self.credentials.expiry and (self.credentials.expiry - datetime.now(timezone.utc)).total_seconds() < 60):
            logging.info("Refreshing user access token.")
            self.credentials.refresh(Request())
        return self.credentials.token

    def build_auth_headers(self, token):
        """
        Construct standard authorization headers.
        """
        return {
            "content-type": "application/json",
            "Authorization": f"Bearer {token}",
        }

    def submit_job(self):
        """
        Submit a job to the specified workspace.
        """
        logging.info(f"Submitting job for method {self.method_namespace}/{self.method_name} in workspace {self.workspace_namespace}/{self.workspace_name}.")
        uri = f"https://api.firecloud.org/api/workspaces/{self.workspace_namespace}/{self.workspace_name}/submissions"
        token = self.get_user_token()
        headers = self.build_auth_headers(token)
        body = {
            "deleteIntermediateOutputFiles": False,
            "methodConfigurationNamespace": self.method_namespace,
            "methodConfigurationName": self.method_name,
            "entityType": self.entity_type,
            "entityName": self.entity_id,
            "useCallCache": False,
        }

        response = requests.post(uri, json=body, headers=headers)
        if response.status_code != 201:
            logging.error(f"Failed to submit job. Status code: {response.status_code}. Response: {response.text}")
            raise Exception("Submission failed.")

        submission_id = response.json().get("submissionId")
        logging.info(f"Job submitted successfully. Submission ID: {submission_id}")
        return submission_id

    def poll_submission_status(self, submission_id):
        """
        Poll the status of the submission until completion.
        """
        logging.info(f"Polling status for submission ID: {submission_id}")
        uri = f"https://api.firecloud.org/api/workspaces/{self.workspace_namespace}/{self.workspace_name}/submissions/{submission_id}"
        token = self.get_user_token()
        headers = self.build_auth_headers(token)

        response = requests.get(uri, headers=headers)
        if response.status_code != 200:
            logging.error(f"Error polling submission status. Status code: {response.status_code}. Response: {response.text}")
            raise Exception("Failed to poll submission status.")

        submission_status = response.json().get("status")
        workflows = response.json().get("workflows", [])
        workflow_status = workflows[0]["status"] if workflows else "Unknown"
        logging.info(f"Submission status: {submission_status}. Workflow status: {workflow_status}")
        return submission_status, workflow_status

    def monitor_submission(self, submission_id, timeout_minutes=60):
        """
        Monitor the submission until it completes or times out.
        """
        sleep_seconds = 60
        max_polls = timeout_minutes * 60 // sleep_seconds

        for _ in range(max_polls):
            submission_status, workflow_status = self.poll_submission_status(submission_id)
            if submission_status == "Done":
                return workflow_status
            sleep(sleep_seconds)

        raise TimeoutError("Monitoring submission timed out.")


# Example Usage
if __name__ == "__main__":
    import argparse

    # Parse arguments passed via CLI or a configuration file like YAML
    parser = argparse.ArgumentParser(description="FireCloud Job Manager")
    parser.add_argument("--sa-json-b64", required=True, help="Base64 encoded service account JSON")
    parser.add_argument("--user", required=True, help="Impersonated user email")
    parser.add_argument("--workspace-namespace", required=True, help="Workspace namespace")
    parser.add_argument("--workspace-name", required=True, help="Workspace name")
    parser.add_argument("--method-namespace", required=True, help="Method configuration namespace")
    parser.add_argument("--method-name", required=True, help="Method configuration name")
    parser.add_argument("--entity-type", required=True, help="Entity type for the job")
    parser.add_argument("--entity-id", required=True, help="Entity ID for the job")
    args = parser.parse_args()

    try:
        # Initialize manager
        manager = FireCloudJobAPI(
            args.sa_json_b64,
            args.user,
            args.workspace_namespace,
            args.workspace_name,
            args.method_namespace,
            args.method_name,
            args.entity_type,
            args.entity_id,
        )

        # Submit job
        submission_id = manager.submit_job()

        # Monitor job
        workflow_status = manager.monitor_submission(submission_id)
        if workflow_status == "Succeeded":
            logging.info("Job completed successfully.")
        else:
            logging.error(f"Job failed with workflow status: {workflow_status}")
            exit(1)

    except Exception as e:
        logging.error(f"Error during execution: {str(e)}")
        traceback.print_exc()
        exit(1)
