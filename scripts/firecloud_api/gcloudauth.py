# gcloud_auth_list.py
import subprocess
import logging

# Configure logging
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')

def gcloud_auth_list():
    try:
        result = subprocess.run(
            ["gcloud", "auth", "list", "--format=json"],
            capture_output=True,
            text=True,
            check=True
        )
        logging.info("gcloud auth list output:")
        print(result.stdout)
    except subprocess.CalledProcessError as e:
        logging.error(f"Error executing 'gcloud auth list': {e.stderr}")
    except Exception as e:
        logging.error(f"An unexpected error occurred: {str(e)}")

if __name__ == "__main__":
    gcloud_auth_list()
