#!/usr/bin/env python3
"""
WDL Command Block Extractor to Bash Scripts

Parses WDL pipeline files using miniwdl and extracts all task command blocks,
exporting each as executable bash scripts with proper variable substitution.
"""

import argparse
import json
import logging
import os
import re
import subprocess
import sys
from pathlib import Path
from typing import Dict, List, Set, Tuple, Optional, Any
import tempfile

try:
    import WDL
except ImportError:
    print("Error: miniwdl is not installed. Install it with: pip install miniwdl")
    sys.exit(1)


class WDLToBashExtractor:
    """Extracts WDL task command blocks and converts them to bash scripts."""

    def __init__(self, output_dir: str = "extracted_tasks", log_level: str = "INFO"):
        """
        Initialize the extractor.

        Args:
            output_dir: Directory to save extracted bash scripts
            log_level: Logging level (DEBUG, INFO, WARNING, ERROR)
        """
        self.output_dir = Path(output_dir)
        self.output_dir.mkdir(parents=True, exist_ok=True)
        
        # Set up logging
        log_format = '%(asctime)s - %(name)s - %(levelname)s - %(message)s'
        logging.basicConfig(level=getattr(logging, log_level), format=log_format)
        self.logger = logging.getLogger(__name__)
        
        # Track processed files to avoid infinite loops
        self.processed_files: Set[str] = set()
        self.task_definitions: Dict[str, Dict[str, Any]] = {}
        self.variable_substitutions: Dict[str, List[Tuple[str, str]]] = {}
        self.bash_errors: Dict[str, str] = {}  # Task name -> error message

    def resolve_import_path(self, import_path: str, base_dir: Path) -> Optional[Path]:
        """
        Resolve WDL import path relative to base directory.

        Args:
            import_path: The import path from WDL file
            base_dir: The directory containing the WDL file

        Returns:
            Resolved Path if found, None otherwise
        """
        # Handle both absolute and relative paths
        if import_path.startswith("/"):
            path = Path(import_path)
        else:
            path = (base_dir / import_path).resolve()

        if path.exists() and path.suffix == ".wdl":
            self.logger.debug(f"Resolved import: {import_path} -> {path}")
            return path
        
        self.logger.warning(f"Could not resolve import: {import_path} from {base_dir}")
        return None

    def extract_imports(self, wdl_content: str, base_dir: Path) -> List[Path]:
        """
        Extract and resolve all import statements from WDL content.

        Args:
            wdl_content: The WDL file content
            base_dir: The directory containing the WDL file

        Returns:
            List of resolved import paths
        """
        imports = []
        # Match: import "path/to/file.wdl" as alias
        # Also handles: import "path/to/file.wdl"
        pattern = r'import\s+"([^"]+)"\s*(?:as\s+\w+)?'
        
        for match in re.finditer(pattern, wdl_content):
            import_path = match.group(1)
            resolved_path = self.resolve_import_path(import_path, base_dir)
            if resolved_path:
                imports.append(resolved_path)
        
        return imports

    def extract_command_block(self, command_str: str) -> str:
        """
        Extract the raw command block from task definition.

        Args:
            command_str: The command section as string

        Returns:
            Clean command block text
        """
        # Remove leading/trailing whitespace and command block delimiters
        # WDL uses <<< >>> or <<< ... >>>
        text = command_str.strip()
        
        # Remove the <<< and >>> markers if present
        if text.startswith("<<<"):
            text = text[3:]
        if text.endswith(">>>"):
            text = text[:-3]
        
        return text.strip()

    def extract_input_variables(self, inputs: List[WDL.Decl]) -> Dict[str, Dict[str, str]]:
        """
        Extract input variables and their types.

        Args:
            inputs: List of WDL input declarations

        Returns:
            Dictionary mapping variable name to type info
        """
        variables = {}
        for decl in inputs:
            var_name = decl.name
            var_type = str(decl.type)
            
            # Determine if optional
            is_optional = decl.type.optional
            
            variables[var_name] = {
                "type": var_type,
                "optional": is_optional,
                "wdl_type": str(decl.type.base_type) if hasattr(decl.type, 'base_type') else var_type
            }
            self.logger.debug(f"Found input: {var_name}: {var_type}")
        
        return variables

    def extract_output_variables(self, outputs: List[WDL.Decl]) -> Dict[str, str]:
        """
        Extract output declarations.

        Args:
            outputs: List of WDL output declarations

        Returns:
            Dictionary mapping output name to type
        """
        outputs_dict = {}
        for decl in outputs:
            outputs_dict[decl.name] = str(decl.type)
            self.logger.debug(f"Found output: {decl.name}: {str(decl.type)}")
        
        return outputs_dict

    def substitute_variables(self, command: str, input_vars: Dict[str, str]) -> Tuple[str, List[Tuple[str, str]]]:
        """
        Substitute WDL variables with bash variables.

        Args:
            command: The command block text
            input_vars: Dictionary of input variable names

        Returns:
            Tuple of (modified_command, list_of_substitutions)
        """
        modified_command = command
        substitutions: List[Tuple[str, str]] = []

        # First, protect $${...} sequences by replacing them temporarily
        protected_patterns = {}
        pattern_count = 0
        
        # Find all $${...} patterns and protect them
        for match in re.finditer(r'\$\$\{([^}]+)\}', modified_command):
            placeholder = f"__PROTECTED_PATTERN_{pattern_count}__"
            protected_patterns[placeholder] = match.group(0)
            modified_command = modified_command.replace(match.group(0), placeholder)
            pattern_count += 1

        # Replace both ~{var} and ${var} syntax (WDL uses both)
        for var_name in sorted(input_vars.keys(), key=len, reverse=True):
            # Match ~{var_name} (WDL newer syntax)
            pattern1 = r'~\{' + re.escape(var_name) + r'\}'
            replacement = f'${var_name}'
            
            new_command = re.sub(pattern1, replacement, modified_command)
            if new_command != modified_command:
                substitutions.append((f'~{{{var_name}}}', replacement))
                modified_command = new_command
                self.logger.debug(f"Substituted: ~{{{var_name}}} -> ${var_name}")

            # Match ${var_name} (WDL old syntax or escaped)
            pattern2 = r'\$\{' + re.escape(var_name) + r'\}'
            new_command = re.sub(pattern2, replacement, modified_command)
            if new_command != modified_command:
                substitutions.append((f'${{{var_name}}}', replacement))
                modified_command = new_command
                self.logger.debug(f"Substituted: ${{{var_name}}} -> ${var_name}")

        # Restore protected patterns (convert $${...} to ${...})
        for placeholder, original in protected_patterns.items():
            # Convert $${...} to ${...} for shell variable syntax
            shell_var = original.replace("$$", "$")
            modified_command = modified_command.replace(placeholder, shell_var)
            substitutions.append((original, shell_var))
            self.logger.debug(f"Preserved shell sequence: {original} -> {shell_var}")

        return modified_command, substitutions

    def generate_bash_script(self, task_name: str, command: str, input_vars: Dict[str, Dict[str, str]], 
                           outputs: Dict[str, str], source_file: str) -> str:
        """
        Generate a complete bash script from task definition.

        Args:
            task_name: Name of the task
            command: The command block
            input_vars: Dictionary of input variables with types
            outputs: Dictionary of output variables
            source_file: Source WDL file path

        Returns:
            Complete bash script as string
        """
        # Perform variable substitution
        substituted_command, substitutions = self.substitute_variables(
            command, 
            {name: info["type"] for name, info in input_vars.items()}
        )
        
        # Store substitutions for logging
        self.variable_substitutions[task_name] = substitutions

        # Build the bash script
        script_lines = []
        
        # Shebang
        script_lines.append("#!/bin/bash")
        script_lines.append("")
        
        # Header comment
        script_lines.append(f"# Generated bash script for WDL task: {task_name}")
        script_lines.append(f"# Source: {source_file}")
        script_lines.append(f"# Generated: {__import__('datetime').datetime.now().isoformat()}")
        script_lines.append("")
        
        # Input documentation
        if input_vars:
            script_lines.append("# Input Variables:")
            script_lines.append("# ===============")
            for var_name, var_info in sorted(input_vars.items()):
                optional_marker = " (optional)" if var_info["optional"] else ""
                script_lines.append(f"# {var_name}: {var_info['type']}{optional_marker}")
            script_lines.append("")
        
        # Output documentation
        if outputs:
            script_lines.append("# Output Variables:")
            script_lines.append("# =================")
            for out_name, out_type in sorted(outputs.items()):
                script_lines.append(f"# {out_name}: {out_type}")
            script_lines.append("")
        
        # Variable declarations
        script_lines.append("# Input Variable Declarations")
        script_lines.append("# ============================")
        for var_name, var_info in sorted(input_vars.items()):
            if var_info["optional"]:
                script_lines.append(var_name + '="${' + var_name + ':-}"  # ' + var_info["type"])
            else:
                script_lines.append(var_name + '="${TODO_' + var_name + '}"  # ' + var_info["type"] + " - REQUIRED")
        script_lines.append("")
        
        # Error handling
        script_lines.append("# Exit on error")
        script_lines.append("set -euo pipefail")
        script_lines.append("")
        
        # Command block
        script_lines.append("# Command Block")
        script_lines.append("# =============")
        script_lines.append(substituted_command)
        script_lines.append("")
        
        return "\n".join(script_lines)

    def validate_bash_syntax(self, script: str, task_name: str) -> Tuple[bool, str]:
        """
        Validate bash script syntax using bash -n.

        Args:
            script: The bash script content
            task_name: Task name (for logging)

        Returns:
            Tuple of (is_valid, error_message)
        """
        try:
            result = subprocess.run(
                ["bash", "-n"],
                input=script.encode(),
                capture_output=True,
                timeout=5
            )
            if result.returncode == 0:
                self.logger.debug(f"Bash syntax validation passed for {task_name}")
                return True, ""
            else:
                error_msg = result.stderr.decode().strip()
                self.logger.warning(f"Bash syntax error in {task_name}: {error_msg}")
                # Track error for test reporting
                self.bash_errors[task_name] = error_msg
                return False, error_msg
        except subprocess.TimeoutExpired:
            error_msg = f"Bash validation timed out for {task_name}"
            self.logger.warning(error_msg)
            self.bash_errors[task_name] = error_msg
            return False, error_msg
        except Exception as e:
            error_msg = f"Error validating bash syntax: {e}"
            self.logger.warning(error_msg)
            self.bash_errors[task_name] = error_msg
            return False, error_msg

    def save_bash_script(self, task_name: str, script: str, validate: bool = True) -> bool:
        """
        Save bash script to file.

        Args:
            task_name: Name of the task
            script: Script content
            validate: Whether to validate syntax before saving

        Returns:
            True if successful, False otherwise
        """
        # Validate if requested
        if validate:
            is_valid, error_msg = self.validate_bash_syntax(script, task_name)
            if not is_valid:
                self.logger.error(f"Skipping {task_name} due to bash syntax error: {error_msg}")
                return False

        # Save script
        script_path = self.output_dir / f"{task_name}.sh"
        try:
            with open(script_path, 'w') as f:
                f.write(script)
            # Make executable
            os.chmod(script_path, 0o755)
            self.logger.info(f"Saved bash script: {script_path}")
            return True
        except Exception as e:
            self.logger.error(f"Error saving script for {task_name}: {e}")
            return False

    def generate_input_json_template(self, task_name: str, input_vars: Dict[str, Dict[str, str]]) -> str:
        """
        Generate JSON template for task inputs.

        Args:
            task_name: Name of the task
            input_vars: Dictionary of input variables

        Returns:
            JSON string with input template
        """
        template = {}
        for var_name, var_info in input_vars.items():
            var_type = var_info["wdl_type"].lower()
            
            # Determine appropriate default based on type
            if "file" in var_type:
                default = f"TODO: path/to/{var_name}"
            elif "string" in var_type:
                default = f"TODO: {var_name}_value"
            elif "int" in var_type:
                default = "TODO: <integer>"
            elif "float" in var_type:
                default = "TODO: <float>"
            elif "boolean" in var_type:
                default = True
            elif "array" in var_type:
                default = ["TODO: array_element"]
            else:
                default = f"TODO: {var_name}"
            
            template[var_name] = default

        return json.dumps(template, indent=2)

    def save_input_template(self, task_name: str, input_vars: Dict[str, Dict[str, str]]) -> bool:
        """
        Save input JSON template.

        Args:
            task_name: Name of the task
            input_vars: Dictionary of input variables

        Returns:
            True if successful, False otherwise
        """
        if not input_vars:
            return True

        json_content = self.generate_input_json_template(task_name, input_vars)
        json_path = self.output_dir / f"{task_name}_inputs.json"
        
        try:
            with open(json_path, 'w') as f:
                f.write(json_content)
            self.logger.info(f"Saved input template: {json_path}")
            return True
        except Exception as e:
            self.logger.error(f"Error saving input template for {task_name}: {e}")
            return False

    def process_task(self, task: WDL.Task, source_file: str, save_files: bool = True) -> bool:
        """
        Process a single task and extract command block.

        Args:
            task: WDL Task object
            source_file: Source WDL file path
            save_files: Whether to save extracted scripts (False for test mode)

        Returns:
            True if successful, False otherwise
        """
        task_name = task.name
        
        self.logger.info(f"Processing task: {task_name} from {source_file}")
        
        # Extract command
        if not hasattr(task, 'command') or task.command is None:
            self.logger.warning(f"Task {task_name} has no command block, skipping")
            return False

        try:
            # Task.command is a WDL.Expr.TaskCommand object
            # Convert to string which reconstructs from parts
            command = self.extract_command_block(str(task.command))
        except Exception as e:
            self.logger.error(f"Error extracting command from {task_name}: {e}")
            return False

        # Extract inputs
        input_vars = self.extract_input_variables(task.inputs)

        # Extract outputs
        outputs = self.extract_output_variables(task.outputs)

        # Store task definition
        self.task_definitions[task_name] = {
            "source": source_file,
            "inputs": input_vars,
            "outputs": outputs,
            "command": command
        }

        # Generate bash script
        bash_script = self.generate_bash_script(task_name, command, input_vars, outputs, source_file)

        # Always validate to catch errors
        is_valid, error_msg = self.validate_bash_syntax(bash_script, task_name)
        
        # Only save if requested and valid (or if save_files is True and validation passed)
        if not save_files:
            # Test mode: just validate, don't save
            return is_valid

        # Extraction mode: save files
        script_saved = self.save_bash_script(task_name, bash_script, validate=False)
        
        if script_saved and input_vars:
            self.save_input_template(task_name, input_vars)

        return script_saved

    def process_wdl_file(self, wdl_path: Path, save_files: bool = True) -> int:
        """
        Process a WDL file and its imports recursively.

        Args:
            wdl_path: Path to WDL file
            save_files: Whether to save extracted scripts (False for test mode)

        Returns:
            Number of tasks processed
        """
        wdl_path = wdl_path.resolve()
        
        # Avoid infinite loops
        if str(wdl_path) in self.processed_files:
            self.logger.debug(f"Already processed: {wdl_path}")
            return 0

        self.processed_files.add(str(wdl_path))
        
        if not wdl_path.exists():
            self.logger.error(f"File not found: {wdl_path}")
            return 0

        self.logger.info(f"Processing WDL file: {wdl_path}")

        # Read file content
        try:
            with open(wdl_path, 'r') as f:
                content = f.read()
        except Exception as e:
            self.logger.error(f"Error reading file {wdl_path}: {e}")
            return 0

        # Parse with miniwdl
        try:
            doc = WDL.parse_document(content)
        except Exception as e:
            self.logger.error(f"Error parsing WDL file {wdl_path}: {e}")
            return 0

        # Process imports first (depth-first)
        base_dir = wdl_path.parent
        imports = self.extract_imports(content, base_dir)
        
        tasks_processed = 0
        for import_path in imports:
            tasks_processed += self.process_wdl_file(import_path, save_files=save_files)

        # Process tasks in this file
        for task in doc.tasks:
            if self.process_task(task, str(wdl_path), save_files=save_files):
                tasks_processed += 1

        return tasks_processed

    def generate_summary_report(self) -> str:
        """
        Generate a summary report of all processing.

        Returns:
            Summary report as string
        """
        report_lines = []
        report_lines.append("=" * 70)
        report_lines.append("WDL TO BASH EXTRACTION SUMMARY")
        report_lines.append("=" * 70)
        report_lines.append("")
        
        report_lines.append(f"Tasks Processed: {len(self.task_definitions)}")
        report_lines.append(f"Output Directory: {self.output_dir.absolute()}")
        report_lines.append(f"Files Processed: {len(self.processed_files)}")
        report_lines.append("")
        
        if self.task_definitions:
            report_lines.append("Tasks Extracted:")
            report_lines.append("-" * 70)
            for task_name, task_info in sorted(self.task_definitions.items()):
                report_lines.append(f"  {task_name}:")
                report_lines.append(f"    Source: {task_info['source']}")
                report_lines.append(f"    Inputs: {len(task_info['inputs'])}")
                report_lines.append(f"    Outputs: {len(task_info['outputs'])}")
                
                if task_name in self.variable_substitutions:
                    subs = self.variable_substitutions[task_name]
                    report_lines.append(f"    Variable Substitutions: {len(subs)}")
                    for old, new in subs[:3]:  # Show first 3
                        report_lines.append(f"      {old} -> {new}")
                    if len(subs) > 3:
                        report_lines.append(f"      ... and {len(subs) - 3} more")
            report_lines.append("")
        
        report_lines.append("=" * 70)
        return "\n".join(report_lines)

    def extract_all(self, wdl_path: Path, validate: bool = True, 
                   generate_inputs: bool = True) -> bool:
        """
        Main entry point to extract all tasks from a WDL file.

        Args:
            wdl_path: Path to primary WDL file
            validate: Whether to validate bash syntax
            generate_inputs: Whether to generate input JSON templates

        Returns:
            True if successful, False otherwise
        """
        try:
            tasks_processed = self.process_wdl_file(wdl_path)
            
            # Print summary
            report = self.generate_summary_report()
            print(report)
            
            if tasks_processed == 0:
                self.logger.error("No tasks were extracted")
                return False
            
            self.logger.info(f"Successfully extracted {tasks_processed} tasks")
            return True
            
        except Exception as e:
            self.logger.error(f"Fatal error during extraction: {e}")
            return False

    def test_wdl_for_bash_errors(self, wdl_path: Path) -> bool:
        """
        Test mode: Scan all tasks in a WDL file for bash syntax errors.
        Does not extract or save scripts, only validates and reports errors.

        Args:
            wdl_path: Path to primary WDL file

        Returns:
            True if no errors found, False if errors detected
        """
        try:
            self.logger.info(f"Starting bash error detection for: {wdl_path}")
            tasks_processed = self.process_wdl_file(wdl_path, save_files=False)
            
            # Print test results
            report = self.generate_test_report()
            print(report)
            
            if tasks_processed == 0:
                self.logger.error("No tasks found in WDL file")
                return False
            
            # Return success only if no bash errors were found
            return len(self.bash_errors) == 0
            
        except Exception as e:
            self.logger.error(f"Fatal error during bash error detection: {e}")
            return False

    def generate_test_report(self) -> str:
        """
        Generate a test report showing bash error detection results.

        Returns:
            Test report as string
        """
        report_lines = []
        report_lines.append("=" * 70)
        report_lines.append("WDL BASH ERROR DETECTION REPORT")
        report_lines.append("=" * 70)
        report_lines.append("")
        
        total_tasks = len(self.task_definitions)
        error_count = len(self.bash_errors)
        valid_count = total_tasks - error_count
        
        report_lines.append(f"Total Tasks Found: {total_tasks}")
        report_lines.append(f"Valid Tasks: {valid_count}")
        report_lines.append(f"Tasks with Bash Errors: {error_count}")
        report_lines.append(f"Files Processed: {len(self.processed_files)}")
        report_lines.append("")
        
        if self.bash_errors:
            report_lines.append("TASKS WITH BASH SYNTAX ERRORS:")
            report_lines.append("-" * 70)
            for task_name, error_msg in sorted(self.bash_errors.items()):
                report_lines.append(f"\n❌ {task_name}")
                report_lines.append(f"   Error: {error_msg}")
                if task_name in self.task_definitions:
                    source = self.task_definitions[task_name].get("source", "unknown")
                    report_lines.append(f"   Source: {source}")
            report_lines.append("")
        else:
            report_lines.append("✅ All tasks passed bash syntax validation!")
            report_lines.append("")
        
        # Success/failure summary
        report_lines.append("=" * 70)
        if error_count == 0:
            report_lines.append("RESULT: PASS - No bash syntax errors detected")
        else:
            report_lines.append(f"RESULT: FAIL - {error_count} task(s) have bash syntax errors")
        report_lines.append("=" * 70)
        
        return "\n".join(report_lines)


def main():
    """Command-line interface."""
    parser = argparse.ArgumentParser(
        description="Extract WDL task command blocks as executable bash scripts or test for bash errors",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Extract all tasks from a pipeline
  python3 wdl_to_bash_extractor.py pipeline.wdl
  
  # Extract with custom output directory
  python3 wdl_to_bash_extractor.py pipeline.wdl -o ./my_scripts
  
  # Extract with debug logging
  python3 wdl_to_bash_extractor.py pipeline.wdl -l DEBUG
  
  # TEST MODE: Scan for bash syntax errors (no extraction)
  python3 wdl_to_bash_extractor.py pipeline.wdl --test
  
  # Extract without input JSON templates
  python3 wdl_to_bash_extractor.py pipeline.wdl --no-input-templates
        """
    )
    
    parser.add_argument(
        "wdl_file",
        help="Path to WDL file to process"
    )
    parser.add_argument(
        "-o", "--output",
        default="extracted_tasks",
        help="Output directory for bash scripts (default: extracted_tasks)"
    )
    parser.add_argument(
        "-l", "--log-level",
        default="INFO",
        choices=["DEBUG", "INFO", "WARNING", "ERROR"],
        help="Logging level (default: INFO)"
    )
    parser.add_argument(
        "--test",
        action="store_true",
        help="Test mode: Scan for bash syntax errors without extracting scripts"
    )
    parser.add_argument(
        "--no-validation",
        action="store_true",
        help="Skip bash syntax validation (ignored in test mode)"
    )
    parser.add_argument(
        "--no-input-templates",
        action="store_true",
        help="Skip generating input JSON templates"
    )
    
    args = parser.parse_args()
    
    # Validate input file
    wdl_path = Path(args.wdl_file)
    if not wdl_path.exists():
        print(f"Error: File not found: {wdl_path}", file=sys.stderr)
        sys.exit(1)
    
    # Create extractor
    extractor = WDLToBashExtractor(
        output_dir=args.output,
        log_level=args.log_level
    )
    
    # Run appropriate mode
    if args.test:
        # Test mode: scan for errors
        success = extractor.test_wdl_for_bash_errors(wdl_path)
    else:
        # Extract mode: extract scripts
        success = extractor.extract_all(
            wdl_path,
            validate=not args.no_validation,
            generate_inputs=not args.no_input_templates
        )
    
    sys.exit(0 if success else 1)


if __name__ == "__main__":
    main()
