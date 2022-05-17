# Scala Test Generator

Quick script to generate the wrapper WDL for the new scala test framework.

The script should generate about 95% of the wrapper wdl for you, you will only need to provide the input_file/files names for the calls to `GetValidationInputs`.

Note that there may be some naming conflicts if the main WDL has `vault_token_path` as input (this is the case in the Arrays pipelines as Arrays needs a vault token to read from the CloudSQL db). In that case
you will have two inputs of the name `vault token_path`. One will need to be changed to `vault_token_path_arrays` and passed as input to the subworkflow. See [TestArrays.wdl](../scripts/TestArrays.wdl) for an example.

## Quick Start

The script only needs two arguments: the workflow to generate the test wrapper for, and the validation wdl for that workflow.

```bash
pip3 install requirements.txt

python3 gentest.py --workflow Arrays --validation VerifyArrays
```

