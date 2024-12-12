import requests
import sys

def fetch_workflow_and_version_ids(token, repository, subclass, version_name):
    url = f"https://dockstore.org/api/workflows/path/workflow/{repository}/published"
    headers = {
        "Authorization": f"Bearer {token}",
        "Accept": "application/json",
    }

    response = requests.get(url, headers=headers)
    response.raise_for_status()
    data = response.json()

    # Extract workflow and version IDs
    workflow_id = data.get("id")
    version_id = next(
        (version["id"] for version in data.get("workflowVersions", [])
         if version["name"] == version_name),
        None
    )
    #increase the version number by 1
    version_id = version_id + 1

    if not workflow_id or not version_id:
        raise ValueError("Workflow ID or Version ID could not be found.")

    return workflow_id, version_id

def fetch_commit_id(token, workflow_id, version_id):
    url = f"https://dockstore.org/api/workflows/{workflow_id}/workflowVersions/{version_id}"
    headers = {
        "Authorization": f"Bearer {token}",
        "Accept": "application/json",
    }

    response = requests.get(url, headers=headers)
    response.raise_for_status()
    data = response.json()

    # Extract commit ID
    commit_id = data.get("commitID")
    if not commit_id:
        raise ValueError("Commit ID could not be found.")

    return commit_id

if __name__ == "__main__":
    if len(sys.argv) != 5:
        print("Usage: python fetch_dockstore_commit.py <token> <repository> <subclass> <version_name>")
        sys.exit(1)

    _, token, repository, subclass, version_name = sys.argv

    try:
        workflow_id, version_id = fetch_workflow_and_version_ids(token, repository, subclass, version_name)
        print(f"Workflow ID: {workflow_id}")
        print(f"Version ID: {version_id}")

        commit_id = fetch_commit_id(token, workflow_id, version_id)
        print(f"Commit ID: {commit_id}")

    except Exception as e:
        print(f"Error: {e}")
