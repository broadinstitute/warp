import requests
import sys

def fetch_workflow_id(token, repository, subclass, version_name):
    url = f"https://dockstore.org/api/workflows/path/workflow/{repository}/published?subclass={subclass}&versionName={version_name}"
    print(url)
    params = {
        "subclass": subclass,
        "versionName": version_name,
    }
    headers = {
        "Authorization": f"Bearer {token}",
        "Accept": "application/json",
    }

    response = requests.get(url, params=params, headers=headers)
    response.raise_for_status()
    workflow_data = response.json()
    print(workflow_data.get("id"))

if __name__ == "__main__":
    if len(sys.argv) != 5:
        print("Usage: python fetch_dockstore_id.py <token> <repository> <subclass> <version_name>")
        sys.exit(1)

    _, token, repository, subclass, version_name = sys.argv
    fetch_workflow_id(token, repository, subclass, version_name)
