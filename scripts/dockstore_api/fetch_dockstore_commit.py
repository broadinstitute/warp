import requests
import sys

def fetch_commit_id(token, repository, version_name):
    # Fetch the workflow data
    url = f"https://dockstore.org/api/workflows/path/workflow/github.com%2Fbroadinstitute%2Fwarp%2F{repository}/published"
    headers = {
        "Authorization": f"Bearer {token}",
        "Accept": "application/json",
    }

    response = requests.get(url, headers=headers)
    response.raise_for_status()
    data = response.json()

    # Extract workflow ID and version ID
    workflow_id = data.get("id")
    version_id = next(
        (version["id"] for version in data.get("workflowVersions", [])
         if version["name"] == version_name),
        None
    )

    if not workflow_id or not version_id:
        raise ValueError("Workflow ID or Version ID could not be found.")

    # Fetch the specific version details to get the commit ID
    version_url = f"https://dockstore.org/api/workflows/{workflow_id}/workflowVersions/{version_id}"
    version_response = requests.get(version_url, headers=headers)
    version_response.raise_for_status()
    version_data = version_response.json()

    # Extract commit ID
    commit_id = version_data.get("commitID")
    if not commit_id:
        raise ValueError("Commit ID could not be found.")

    return commit_id

if __name__ == "__main__":
    if len(sys.argv) != 4:
        print("Usage: python fetch_dockstore_commit.py <token> <repository> <version_name>")
        sys.exit(1)

    _, token, repository, version_name = sys.argv

    try:
        commit_id = fetch_commit_id(token, repository, version_name)
        print(commit_id)
    except Exception as e:
        print(f"Error: {e}")