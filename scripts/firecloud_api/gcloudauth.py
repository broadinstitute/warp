import subprocess
import logging
import os

def gcloud_auth_list(service_account_key_file=None):
    try:
        # If a service account key file is provided, activate the service account
        if service_account_key_file and os.path.exists(service_account_key_file):
            logging.info(f"Activating service account using key file: {service_account_key_file}")
            # Command to activate the service account
            subprocess.run(
                ["gcloud", "auth", "activate-service-account", "--key-file", service_account_key_file],
                check=True
            )
        else:
            logging.warning("No service account key file provided or the file does not exist.")

        # List authenticated accounts
        logging.info("Listing authenticated accounts:")
        result = subprocess.run(
            ["gcloud", "auth", "list"],
            capture_output=True, text=True, check=True
        )

        # Display the output
        if result.stdout.strip():
            logging.info(f"gcloud auth list output:\n{result.stdout}")
        else:
            logging.warning("No authenticated accounts found.")

    except subprocess.CalledProcessError as e:
        logging.error(f"An error occurred while executing gcloud command: {e}")
        logging.error(f"Error details: {e.stderr}")

# Example usage:
# Replace 'YOUR_SERVICE_ACCOUNT_KEY.json' with the actual path to your service account key file
gcloud_auth_list("path_to_your_service_account_key.json")
