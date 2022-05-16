#!/usr/bin/env python

from __future__ import annotations
from collections import namedtuple
from jinja2 import Environment, FileSystemLoader
import argparse
import logging
import sys
import os
import re

logging.basicConfig(level=logging.INFO)

# Get the path in import statement excluding ../../
reg_expression = re.compile(r"([/A-Za-z]+.wdl)")

WdlImport = namedtuple("WdlImport", "path, name")
WdlInput = namedtuple("WdlInput", "type, name")
SubWorkflowInput = namedtuple("SubWorkflowInput", "name")


class TestGenerator:
    PIPELINE_ROOT = "../../../pipelines"
    VALIDATION_ROOT = "../../../verification"

    def __init__(self, workflow: str, validation: str) -> None:
        # Names of the main workflow and validation wdl
        self.workflow = workflow
        self.validation = validation

        # Paths for the workflow and validation wdls
        self.workflow_path = None
        self.validation_path = None

        # {
        #   array_file : [],
        #   file : [] ,
        #   optional : []
        # }
        self.workflow_outputs = dict()
        self.workflow_metrics = dict()

        # Hold values to be passed to template
        (
            self.imports,
            self.workflow_text,
            self.workflow_inputs,
            self.raw_inputs,
            self.raw_outputs,
            self.subworkflow_inputs,
        ) = (list(), list(), list(), list(), list(), list())

        # Get the absolute path for the provided workflow
        for path, _, files in os.walk(self.PIPELINE_ROOT):
            for f in files:
                if f == f"{workflow}.wdl":
                    workflow_path = os.path.join(os.path.abspath(path), f)
                    logging.info(f" - Workflow found at path -> {workflow_path}")
                    self.workflow_path = workflow_path

        # Get the absolute path for the provided validation wdl
        for path, _, files in os.walk(self.VALIDATION_ROOT):
            for f in files:
                if f == f"{validation}.wdl":
                    validation_path = os.path.join(os.path.abspath(path), f)
                    logging.info(
                        f" - Validation wdl found at path -> {validation_path}"
                    )
                    self.validation_path = validation_path

        if not self.workflow or not self.validation:
            raise FileNotFoundError(
                f"- Error unable to find one or more of the following files: {workflow}.wdl , {validation}.wdl"
            )

        return

    def parse_workflow(self) -> TestGenerator:
        """Parses the wdl specified, strips out newline characters, empty lines and comments
        then inserts each line in into self.workflow_text

        Args:
            self (TestGenerator): The current instance of class TestGenerator

        Returns:
            self (TestGenerator): Returns the current instance of class TestGenerator for method chaining
        """

        with open(self.workflow_path, "r") as f:
            self.workflow_text = [
                line.strip()
                for line in f.readlines()
                if line.strip() and not line.strip().startswith("#")
            ]

        return self

    def get_imports(self) -> TestGenerator:
        """Parses the wdl specified to grab all of the import statements and makes them relative to
        the verification/test-wdls dir. Then adds them as a list of WdlImports to self.imports

        Args:
            self (TestGenerator): The current instance of class TestGenerator

        Returns:
            self (TestGenerator): Returns the current instance of class TestGenerator for method chaining
        """

        imports = [
            line.split() for line in self.workflow_text if line.startswith("import")
        ]

        for i in imports:
            match = reg_expression.findall(i[1])[0]
            rel_path = f"../..{match}"
            self.imports.append(WdlImport(rel_path, i[-1]))

        # Add the main workflow as an import
        main_workflow_path = f"../..{self.workflow_path.split('/warp')[-1]}"
        self.imports.append(WdlImport(main_workflow_path, self.workflow))

        return self

    def get_inputs(self) -> TestGenerator:
        """Parses the wdl specified to grab all of the workflow inputs, includes defaults. Then adds them as list of WdlInputs
        to self.inputs

        Args:
            self (TestGenerator): The current instance of class TestGenerator

        Returns:
            self (TestGenerator): Returns the current instance of class TestGenerator for method chaining
        """
        start = end = 0
        text = self.workflow_text

        for i in range(len(text)):
            if text[i].startswith("input {"):
                start = i
                while not text[i].startswith("}"):
                    i += 1
                end = i
                break

        # Add raw lines to instance
        self.raw_inputs = text[start + 1 : end]

        # Spit into tuple and include default (Type, Name)
        # e.g. (Boolean, read_fingerprint_from_mercury = false)
        self.workflow_inputs = [
            WdlInput(t, n) for t, n in (line.split(None, 1) for line in self.raw_inputs)
        ]

        return self

    def get_subworkflow_inputs(self) -> TestGenerator:
        """Adds only the names of the workflow inputs as SubWorkflowInputs to self.subworkflow_inputs

        Args:
            self (TestGenerator): The current instance of class TestGenerator

        Returns:
            self (TestGenerator): Returns the current instance of class TestGenerator for method chaining
        """

        # Spit into tuple without default (Name)
        # e.g. (read_fingerprint_from_mercury)
        self.subworkflow_inputs = [
            SubWorkflowInput(l[1]) for l in (line.split() for line in self.raw_inputs)
        ]

        return self

    def get_outputs(self) -> TestGenerator:
        """ Parses the output of the main workflow and seperates it into regular outputs and metrics outputs

        Those outputs are then organized into a dict based on the type: Array[File], File, File?
        """

        start = end = 0
        text = self.workflow_text

        for i in range(len(text)):
            if text[i].startswith("output {"):
                start = i
                while not text[i].startswith("}"):
                    i += 1
                end = i
                break

        # Add raw lines to instance
        self.raw_outputs = text[start + 1 : end]

        for line in self.raw_outputs:
            if "_metrics" in line:
                self._parse_metric(line)
            else:
                self._parse_output(line)
        print(self.workflow_outputs)
        return self

    def _parse_metric(self, line: str) -> None:
        vals = line.split()
        t, n = vals[0], vals[1]

        if t == "Array[File]":
            self.workflow_metrics.update(
                {"array_file": [n] + self.workflow_metrics.get("array_file", [])}
            )
        elif t == "File":
            self.workflow_metrics.update(
                {"file": [n] + self.workflow_metrics.get("file", [])}
            )
        elif t == "File?":
            self.workflow_metrics.update(
                {"optional": [n] + self.workflow_metrics.get("optional", [])}
            )

        return

    def _parse_output(self, line: str) -> None:
        vals = line.split()
        t, n = vals[0], vals[1]
        print(t, n)

        if t == "Array[File]":
            self.workflow_outputs.update(
                {"array_file": [n] + self.workflow_outputs.get("array_file", [])}
            )
        elif t == "File":
            self.workflow_outputs.update(
                {"file": [n] + self.workflow_outputs.get("file", [])}
            )
        elif t == "File?":
            self.workflow_outputs.update(
                {"optional": [n] + self.workflow_outputs.get("optional", [])}
            )
        return

    def get_validation_inputs(self) -> TestGenerator:
        return self


def main(args):
    generator = TestGenerator(workflow=args.workflow, validation=args.validation)
    (
        generator.parse_workflow()
        .get_imports()
        .get_inputs()
        .get_subworkflow_inputs()
        .get_outputs()
        .get_validation_inputs()
    )

    tester_name = f"Test{args.workflow}"
    jinja_env = Environment(loader=FileSystemLoader("templates"))
    template = jinja_env.get_template("TestTemplate.wdl.j2")
    content = template.render(
        workflow=args.workflow,
        imports=generator.imports,
        inputs=generator.workflow_inputs,
        test_workflow=tester_name,
        subworkflow_inputs=generator.subworkflow_inputs,
    )

    with open("foo.wdl", "w") as f:
        f.write(content)


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
