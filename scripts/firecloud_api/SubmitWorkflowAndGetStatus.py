import requests
import time
import argparse
import json

def create_submission(token, workspace_namespace, workspace_name, submission_data):
    # API endpoint
    base_url = f'https://api.firecloud.org/api/workspaces/{workspace_namespace}/{workspace_name}/submissions'

    # Headers to make API requests
    headers = {
        'accept': 'application/json',
        'Authorization': f'Bearer {token}',
        'Content-Type': 'application/json'
    }

    # Create the submission
    # send an HTTP POST request to the API endpoint
    response = requests.post(base_url, headers=headers, json=submission_data)
    # convert the response json into a dictionary
    submission_response = response.json()
    # extract the submission ID from the response dictionary
    submission_id = submission_response.get("submissionId")

    if not submission_id:
        print("Failed to create submission.")
    else:
        print(f"Submission created with ID: {submission_id}")
        return submission_id

def poll_submission_status(token, workspace_namespace, workspace_name, submission_id):

    # Status endpoint
    status_url = f'https://api.firecloud.org/api/workspaces/{workspace_namespace}/{workspace_name}/submissions/{submission_id}'

    # Headers to make API requests
    headers = {
        'accept': 'application/json',
        'Authorization': f'Bearer {token}'
    }

    # polling the submission status
    # create an emptu list to store the previous workflow status
    previous_workflow_status = []

    # loop until the submission is done
    while True:
        # send a get request and convert the response json into a dictionary
        status_response = requests.get(status_url, headers=headers)
        status_data = status_response.json()

        # get the submission status
        submission_status = status_data.get("status")
        # get the workflow status of each workflow in the submission
        workflows_status = [workflow.get("status") for workflow in status_data.get("workflows", [])]

        # print the workflow status to stdout if it has changed
        if workflows_status != previous_workflow_status:
            print(f"Workflows Status: {workflows_status}")
            previous_workflow_status = workflows_status

        # Check if the submission has completed
        if submission_status == "Done" and "Failed" in workflows_status:
            print("At least one workflow has failed.")
            break
        elif submission_status == "Done":
            break

        # Wait for 10 seconds before polling again
        time.sleep(10)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Submit and monitor a job.')
    parser.add_argument('--token', required=True, help='API access token')
    parser.add_argument('--workspace-namespace', required=True, help='Workspace namespace')
    parser.add_argument('--workspace-name', required=True, help='Workspace name')
    parser.add_argument('--submission-data-file', required=True, help='Path to the JSON file containing submission data')

    args = parser.parse_args()

    # load submission data from JSON file
    with open(args.submission_data_file, 'r') as file:
        submission_data = json.load(file)

    # create submission and get submission ID
    submission_id = create_submission(args.token, args.workspace_namespace, args.workspace_name, submission_data)

    if submission_id:
        # Poll submission status
        poll_submission_status(args.token, args.workspace_namespace, args.workspace_name, submission_id)
