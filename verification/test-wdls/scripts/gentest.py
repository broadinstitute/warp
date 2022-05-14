#!/usr/bin/env python

from __future__ import annotations
import argparse
import logging
import sys
import os

logging.basicConfig(level=logging.INFO)


class TestGenerator:
    PIPELINE_ROOT = "../../../pipelines"
    VALIDATION_ROOT = "../../../verification"

    def __init__(self, workflow: str, validation: str) -> None:
        self.imports = None
        self.workflow = None
        self.validation = None
        self.workflow_text = None
        self.workflow_inputs = None
        self.workflow_outputs = None

        # Get the absolute path for the provided workflow
        for path, _, files in os.walk(self.PIPELINE_ROOT):
            for f in files:
                if f == f"{workflow}.wdl":
                    workflow_path = os.path.join(os.path.abspath(path), f)
                    logging.info(f" - Workflow found at path -> {workflow_path}")
                    self.workflow = workflow_path

        # Get the absolute path for the provided validation wdl
        for path, _, files in os.walk(self.VALIDATION_ROOT):
            for f in files:
                if f == f"{validation}.wdl":
                    validation_path = os.path.join(os.path.abspath(path), f)
                    logging.info(
                        f" - Validation wdl found at path -> {validation_path}"
                    )
                    self.validation = validation_path

        if not self.workflow or not self.validation:
            raise FileNotFoundError(
                f"- Error unable to find one or more of the following files: {workflow}.wdl , {validation}.wdl"
            )

    def parse_workflow(self) -> TestGenerator:
        """Parses the wdl specified, strips out newline characters, empty lines and comments
        then inserts each line in into self.workflow_text

        Args:
            self (TestGenerator): The current instance of class TestGenerator

        Returns:
            self (TestGenerator): Returns the current instance of class TestGenerator for method chaining
        """

        with open(self.workflow, "r") as f:
            self.workflow_text = [
                line.strip()
                for line in f.readlines()
                if line.strip() and not line.strip().startswith("#")
            ]

        for line in self.workflow_text:
            print(f"{line}")
        return self


def main(args):
    generator = TestGenerator(args.workflow, args.validation)
    generator.parse_workflow()


def _parse_args(argv):
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--workflow",
        "-w",
        required=True,
        dest="workflow",
        help="The name of the workflow to generate the test WDL for ",
    )
    parser.add_argument(
        "--validation",
        "-v",
        required=True,
        dest="validation",
        help="The name of the validation WDL for this this workflow",
    )
    return parser.parse_args(argv)


if __name__ == "__main__":
    args = _parse_args(sys.argv[1:])
    main(args)
