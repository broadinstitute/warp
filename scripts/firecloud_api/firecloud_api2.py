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
USER = os.getenv("USER")
WORKSPACE_NAMESPACE = os.getenv("WORKSPACE_NAMESPACE")
WORKSPACE_NAME = os.getenv("WORKSPACE_NAME")
METHOD_NAMESPACE = os.getenv("METHOD_NAMESPACE")
METHOD_NAME = os.getenv("METHOD_NAME")
ENTITY_TYPE = os.getenv("ENTITY_TYPE")
ENTITY_ID = os.getenv("ENTITY_ID")
sa_json_b64 = os.environ.get("SA_JSON_B64")
print(sa_json_b64)

# Configure logging
LOG_FORMAT = "%(asctime)s %(levelname)-8s %(message)s"
LOG_LEVEL = "INFO"
logging.basicConfig(
    format=LOG_FORMAT,
    level=getattr(logging, LOG_LEVEL),
    datefmt="%Y-%m-%d %H:%M:%S",
)

def get_user_token(credentials: credentials.Credentials):
    """
    Obtain and refresh the user access token.
    """
    if not credentials.valid or (credentials.expiry and (credentials.expiry - datetime.now(timezone.utc)).total_seconds() < 60):
        logging.info("Refreshing user access token.")
        credentials.refresh(Request())
    return credentials.token

def build_auth_headers(token: str):
    """
    Construct standard authorization headers.
    """
    return {
        "content-type": "application/json",
        "Authorization": f"Bearer {token}",
    }

def submit_job(credentials: credentials.Credentials):
    """
    Submit a job to the specified workspace.
    """
    logging.info(f"Submitting job for method {METHOD_NAMESPACE}/{METHOD_NAME} in workspace {WORKSPACE_NAMESPACE}/{WORKSPACE_NAME}.")
    uri = f"https://api.firecloud.org/api/workspaces/{WORKSPACE_NAMESPACE}/{WORKSPACE_NAME}/submissions"
    token = get_user_token(credentials)
    headers = build_auth_headers(token)
    body = {
        "deleteIntermediateOutputFiles": False,
        "methodConfigurationNamespace": METHOD_NAMESPACE,
        "methodConfigurationName": METHOD_NAME,
        "entityType": ENTITY_TYPE,
        "entityName": ENTITY_ID,
        "useCallCache": False,
    }

    response = requests.post(uri, json=body, headers=headers)
    if response.status_code != 201:
        logging.error(f"Failed to submit job. Status code: {response.status_code}. Response: {response.text}")
        raise Exception("Submission failed.")

    submission_id = response.json().get("submissionId")
    logging.info(f"Job submitted successfully. Submission ID: {submission_id}")
    return submission_id

def poll_submission_status(credentials: credentials.Credentials, submission_id: str):
    """
    Poll the status of the submission until completion.
    """
    logging.info(f"Polling status for submission ID: {submission_id}")
    uri = f"https://api.firecloud.org/api/workspaces/{WORKSPACE_NAMESPACE}/{WORKSPACE_NAME}/submissions/{submission_id}"
    token = get_user_token(credentials)
    headers = build_auth_headers(token)

    response = requests.get(uri, headers=headers)
    if response.status_code != 200:
        logging.error(f"Error polling submission status. Status code: {response.status_code}. Response: {response.text}")
        raise Exception("Failed to poll submission status.")

    submission_status = response.json().get("status")
    workflows = response.json().get("workflows", [])
    workflow_status = workflows[0]["status"] if workflows else "Unknown"
    logging.info(f"Submission status: {submission_status}. Workflow status: {workflow_status}")
    return submission_status, workflow_status

def monitor_submission(credentials: credentials.Credentials, submission_id: str, timeout_minutes=60):
    """
    Monitor the submission until it completes or times out.
    """
    sleep_seconds = 60
    max_polls = timeout_minutes * 60 // sleep_seconds

    for _ in range(max_polls):
        submission_status, workflow_status = poll_submission_status(credentials, submission_id)
        if submission_status == "Done":
            return workflow_status
        sleep(sleep_seconds)

    raise TimeoutError("Monitoring submission timed out.")

def main():
    """
    Main workflow execution function.
    """
    try:
        logging.info("Starting job submission and monitoring process.")
        scopes = ["profile", "email", "openid"]
        #decoded_sa = base64.b64decode(SA_JSON_B64).decode("utf-8")
        decoded_sa = SA_JSON_B64
        sa_credentials = service_account.Credentials.from_service_account_info(
            json.loads(decoded_sa), scopes=scopes
        )
        delegated_credentials = sa_credentials.with_subject(USER)

        # Submit the job and monitor its progress
        submission_id = submit_job(delegated_credentials)
        workflow_status = monitor_submission(delegated_credentials, submission_id)

        if workflow_status == "Succeeded":
            logging.info("Job completed successfully.")
        else:
            logging.error(f"Job failed with workflow status: {workflow_status}")
            exit(1)
    except Exception as e:
        logging.error(f"Error during execution: {str(e)}")
        traceback.print_exc()
        exit(1)

if __name__ == "__main__":
    main()
