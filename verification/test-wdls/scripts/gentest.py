#!/usr/bin/env python

from __future__ import annotations
from typing import List
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
        self.validation_inputs = dict()

        # Hold values to be passed to template
        (
            self.imports,
            self.workflow_text,
            self.validation_text,
            self.workflow_inputs,
            self.raw_validation_inputs,
            self.raw_inputs,
            self.raw_outputs,
            self.subworkflow_inputs,
        ) = (list(), list(), list(), list(), list(), list(), list(), list())

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
        then inserts each line in into self.workflow_text"""

        with open(self.workflow_path, "r") as f:
            self.workflow_text = [
                line.strip()
                for line in f.readlines()
                if line.strip() and not line.strip().startswith("#")
            ]

        return self

    def parse_validation(self) -> TestGenerator:
        """Parses the valdiation wdl specified, strips out newline characters, empty lines and comments
        then inserts each line in into self.validation_text"""

        with open(self.validation_path, "r") as f:
            self.validation_text = [
                line.strip()
                for line in f.readlines()
                if line.strip() and not line.strip().startswith("#")
            ]

        return self

    def get_imports(self) -> TestGenerator:
        """Add the necessary imports to the test wdl: Copy, Utilities, main Wdl and Validation Wdl. Append them
        to self.imports"""

        main_workflow_path = f"../..{self.workflow_path.split('/warp')[-1]}"
        validation_workflow_path = f"../..{self.validation_path.split('warp')[-1]}"

        self.imports.append(WdlImport(main_workflow_path, self.workflow))
        self.imports.append(WdlImport(validation_workflow_path, self.validation))
        self.imports.append(WdlImport("../../tasks/wdl/Utilities.wdl", "Utilities"))
        self.imports.append(
            WdlImport("../../tasks/wdl/CopyFilesFromCloudToCloud.wdl", "Copy")
        )

        return self

    def _input_indexes(self, text: List[str]) -> (int, int):
        """Get the start and end index of an input{} section"""

        start = end = 0

        for i in range(len(text)):
            if text[i].startswith("input {"):
                start = i
                while not text[i].startswith("}"):
                    i += 1
                end = i
                break

        return start, end

    def _output_indexes(self, text: List[str]) -> (int, int):
        """Get the start and end index of an output{} section"""

        start = end = 0

        for i in range(len(text)):
            if text[i].startswith("output {"):
                start = i
                while not text[i].startswith("}"):
                    i += 1
                end = i
                break

        return start, end

    def get_inputs(self) -> TestGenerator:
        """Parses the wdl specified to grab all of the workflow inputs, includes defaults. Then adds them as list of WdlInputs
        to self.inputs"""

        text = self.workflow_text
        start, end = self._input_indexes(text)

        # Add raw lines to instance
        self.raw_inputs = text[start + 1 : end]

        # Spit into tuple and include default (Type, Name)
        # e.g. (Boolean, read_fingerprint_from_mercury = false)
        self.workflow_inputs = [
            WdlInput(t, n) for t, n in (line.split(None, 1) for line in self.raw_inputs)
        ]

        logging.info(f" - Generating inputs for Test{self.workflow}.wdl...")

        return self

    def get_validation_inputs(self) -> TestGenerator:
        """Parses the validation wdl specified to grab all of the validation wdl inputs
        Also formats them in PascalCase and checks whether they are a File of Array[File]"""

        text = self.validation_text
        start, end = self._input_indexes(text)

        self.raw_validation_inputs = text[start + 1 : end]

        for i in self.raw_validation_inputs:
            t, v = i.split()
            is_file = t == "File"

            # If we aren't at a paired input then skip
            if not v.startswith("truth_") and not v.startswith("test_"):
                continue

            _, input_name = v.split("_", 1)

            # check if we have an entry for this validation pair, if not then add it
            if input_name not in self.validation_inputs:
                # convert red_idat_md5 to RedIdatMd5
                formatted_name = "".join(
                    [text.capitalize() for text in input_name.split("_")]
                )
                self.validation_inputs[input_name] = {
                    "value": input_name,
                    "format_name": formatted_name,
                    "is_file": is_file,
                }

        logging.info(
            f" - Collecting and formatting inputs for {self.validation}.wdl..."
        )

    def get_subworkflow_inputs(self) -> TestGenerator:
        """Adds only the names of the workflow inputs as SubWorkflowInputs to self.subworkflow_inputs"""

        # Spit into tuple without default (Name)
        # e.g. (read_fingerprint_from_mercury)
        self.subworkflow_inputs = [
            SubWorkflowInput(l[1]) for l in (line.split() for line in self.raw_inputs)
        ]

        logging.info(f" - Generating subworkflow inputs for {self.workflow}.wdl...")

        return self

    def get_outputs(self) -> TestGenerator:
        """ Parses the output of the main workflow and seperates it into regular outputs and metrics outputs"""

        text = self.workflow_text
        start, end = self._output_indexes(text)

        # Add raw lines to instance
        self.raw_outputs = text[start + 1 : end]

        for line in self.raw_outputs:
            if "_metrics" in line:
                self._parse_metric(line)
            else:
                self._parse_output(line)
        return self

        logging.info(f" - Parsing metrics and outputs for {self.workflow}.wdl...")

    def _parse_metric(self, line: str) -> None:
        """Create a list for each type of metric output and map it to that type
        {
            array_file : [],
            file: [],
            optional: []
        }
        """

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
        """Create a list for each type of regular output and map it to that type
        {
            array_file : [],
            file: [],
            optional: []
        }
        """

        vals = line.split()
        t, n = vals[0], vals[1]

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


def main(args):
    generator = TestGenerator(workflow=args.workflow, validation=args.validation)
    (
        generator.parse_workflow()
        .parse_validation()
        .get_imports()
        .get_inputs()
        .get_subworkflow_inputs()
        .get_outputs()
        .get_validation_inputs()
    )

    tester_name = f"Test{args.workflow}"
    jinja_env = Environment(loader=FileSystemLoader("templates"))
    template = jinja_env.get_template("TestTemplate.wdl.j2")

    logging.info(f" - Building {tester_name}.wdl...")

    content = template.render(
        workflow=args.workflow,
        imports=generator.imports,
        inputs=generator.workflow_inputs,
        outputs=generator.workflow_outputs,
        metrics=generator.workflow_metrics,
        test_workflow=tester_name,
        validation=generator.validation,
        subworkflow_inputs=generator.subworkflow_inputs,
        validation_inputs=generator.validation_inputs,
    )

    with open(f"{tester_name}.wdl", "w") as f:
        f.write(content)

    logging.info(f" - Successfully generated {tester_name}.wdl!")


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
