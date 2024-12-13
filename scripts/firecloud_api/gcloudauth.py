import logging
import subprocess
import os
import base64
import tempfile

# Set up logging to print output to both console and file (optional)
logging.basicConfig(level=logging.DEBUG, format='%(asctime)s - %(levelname)s - %(message)s')

def gcloud_auth_list(sa_json_b64):
    try:
        # Decode the Base64 service account key
        if sa_json_b64:
            logging.info("Decoding service account JSON...")
            decoded_json = base64.b64decode(sa_json_b64).decode('utf-8')

            # Create a temporary file to store the decoded JSON
            with tempfile.NamedTemporaryFile(delete=False, mode='w', suffix='.json') as tmp_key_file:
                tmp_key_file.write(decoded_json)
                tmp_key_file_path = tmp_key_file.name

            logging.info(f"Service account key file created at: {tmp_key_file_path}")

            # Activate the service account using the decoded JSON file
            logging.info(f"Activating service account using key file: {tmp_key_file_path}")
            result = subprocess.run(
                ["gcloud", "auth", "activate-service-account", "--key-file", tmp_key_file_path],
                capture_output=True, text=True, check=True
            )

            # Check if activation was successful
            if result.returncode == 0:
                logging.info("Service account activated successfully.")
            else:
                logging.error("Failed to activate service account.")
                logging.error(result.stderr)

            # Now list authenticated accounts
            logging.info("Listing authenticated accounts:")
            result = subprocess.run(
                ["gcloud", "auth", "list"],
                capture_output=True, text=True, check=True
            )

            # Print the output directly to console
            print("gcloud auth list output:")
            print(result.stdout)

            # Log the output
            if result.stdout.strip():
                logging.info(f"gcloud auth list output:\n{result.stdout}")
            else:
                logging.warning("No authenticated accounts found.")

            # Clean up the temporary key file
            os.remove(tmp_key_file_path)
            logging.info("Temporary service account key file removed.")

        else:
            logging.error("No service account key (Base64) provided.")

    except subprocess.CalledProcessError as e:
        logging.error(f"An error occurred while executing gcloud command: {e}")
        logging.error(f"Error details: {e.stderr}")
    except Exception as e:
        logging.error(f"An unexpected error occurred: {e}")

# Example usage (you would pass SA_JSON_B64 from environment):
# gcloud_auth_list(SA_JSON_B64)
