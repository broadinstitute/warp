# Workflow Tests
These is a test suite designed to make testing cromwell workflows easier. It contains utilities that makes writing tests simple and provides layers of abstraction to avoid rewriting the same code over and over again.

## Environments
Predefined environments, implemented by exetending the trait `CromwellEnvironment`, allow the workflow tests to be run in multiple environments

### Dev
The dev environment tells the test suite to run tests on the gotc dev cromwell instance. It will make use of dev clio and firecloud-alpha if they're needed. This is the environment that will most likely be used.

### Prod
The prod environment tells the test suite to run tests on the gotc production cromwell instance. It will make use of dev clio and firecloud-alpha if they're needed.

## Tests

There are a smaller number of tests that test the workflow test framework itself. They can be run using sbt:
```$bash
$ sbt test
```

### Dummy  workflow

There is a dummy test that can be run to verify that the test suite is working or that a cromwell instance is correctly up and running. The dummy WDL simply prints a line to the console, and the test checks that the line was printed. This is the default if no test is specified, and it defaults to running in the `Dev` environment.

```$bash
$ java -jar workflow_tests.jar
```

### Cram Data Delivery
```
Command: CramDataDelivery [options]
Test the Cram data-delivery workflow
  -d <value> | --data-type <value>
        Data type to test delivery of [WGS, Exome]
  -r <value> | --requester <value>
        Email address of the requester so that they can access a workspace created with the --leave-workspace option
  -e <value> | --env <value>
        The environment that this should run in [staging|jgdev|pharma5|jgprod|dev|prod]
  --leave-workspace
        Leave the workspace in firecloud for manual verification
```

The cramDataDelivery test will launch a workflow that creates a FireCloud workspace and test that the workspace attributes, participants, sample sets, samples, and tags are correct. It will then clean up after itself, deleting the workspace and moving files back to their original location. If the revert fails, there is a script in `src/main/scripts/cramDataDelivery/cleanup` that can help.
```$bash
$ java -jar workflow_tests.jar CramDataDelivery -e Dev
````

### Germline Single Sample
```
Command: GermlineSingleSample [options]
Test the GermlineSingleSample workflow
  -t <value> | --test <value>
        The type of test to run
  -d <value> | --data-type <value>
        The data type to test [WGS,Exome]
  -b <value> | --branch <value>
        The branch of truth data to test against (Defaults to develop)
  --use-timestamp <value>
        Do not run the workflows. Instead, just use a previous runs timestamp (yyyy-MM-dd-HH-mm-ss)
  -u | --uncached
        Disable call-caching for the main workflow run
  -e <value> | --env <value>
        The environment that this should run in [staging|jgdev|pharma5|jgprod|dev|prod]
  --update-truth
        Update the truth data with the results of this run.
  --papi-version <value>
        The version of Pipelines API to use
```
This test will kick off the GermlineSingleSample test of your choice. It runs the WDL and then checks the resulting data against known truth data.
