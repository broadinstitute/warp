#!/usr/bin/env python3
"""
Test suite for WDL to Bash Extractor
Tests various WDL parsing scenarios and script generation
"""

import tempfile
import subprocess
import os
from pathlib import Path


def test_simple_task():
    """Test extraction of a simple task with basic variables"""
    wdl_content = '''
    task SimpleTask {
      input {
        String name
        File input_file
        Int count
        Boolean flag
      }
      command <<<
        echo "Processing ${name}"
        wc -l ${input_file}
        for i in $(seq 1 ${count}); do
          echo "Iteration $i"
        done
      >>>
      output {
        String result = read_string(stdout())
      }
    }
    '''
    
    print("Test 1: Simple Task Extraction")
    print("=" * 50)
    
    with tempfile.TemporaryDirectory() as tmpdir:
        wdl_file = Path(tmpdir) / "simple.wdl"
        wdl_file.write_text(wdl_content)
        
        # Run extractor
        result = subprocess.run(
            ["python3", "wdl_to_bash_extractor.py", str(wdl_file), 
             "-o", f"{tmpdir}/output"],
            capture_output=True,
            text=True
        )
        
        if result.returncode == 0:
            script = Path(f"{tmpdir}/output/SimpleTask.sh")
            if script.exists():
                print("✓ Script generated successfully")
                content = script.read_text()
                
                # Check for key elements
                checks = [
                    ("Shebang", "#!/bin/bash" in content),
                    ("Variables declared", "$name" in content and "$input_file" in content),
                    ("Command block", "Processing" in content and "Iteration" in content),
                    ("Set error handling", "set -euo pipefail" in content),
                ]
                
                for check_name, passed in checks:
                    print(f"  {'✓' if passed else '✗'} {check_name}")
            else:
                print("✗ Script not created")
        else:
            print(f"✗ Extraction failed: {result.stderr}")
    
    print()


def test_optional_inputs():
    """Test handling of optional input types"""
    wdl_content = '''
    task OptionalTask {
      input {
        String required_input
        String? optional_string
        File? optional_file
        Int? optional_int
      }
      command <<<
        set -e
        if [ -n "${optional_string}" ]; then
          echo "Got optional: ${optional_string}"
        fi
        if [ -f "${optional_file}" ]; then
          ls -la ${optional_file}
        fi
      >>>
      output {
        String status = "complete"
      }
    }
    '''
    
    print("Test 2: Optional Input Types")
    print("=" * 50)
    
    with tempfile.TemporaryDirectory() as tmpdir:
        wdl_file = Path(tmpdir) / "optional.wdl"
        wdl_file.write_text(wdl_content)
        
        result = subprocess.run(
            ["python3", "wdl_to_bash_extractor.py", str(wdl_file),
             "-o", f"{tmpdir}/output"],
            capture_output=True,
            text=True
        )
        
        if result.returncode == 0:
            script = Path(f"{tmpdir}/output/OptionalTask.sh")
            if script.exists():
                print("✓ Script generated successfully")
                content = script.read_text()
                
                checks = [
                    ("Required variable marked", "TODO_required_input" in content),
                    ("Optional variable with default", "optional_string:-" in content),
                    ("Optional file with default", "optional_file:-" in content),
                ]
                
                for check_name, passed in checks:
                    print(f"  {'✓' if passed else '✗'} {check_name}")
            else:
                print("✗ Script not created")
        else:
            print(f"✗ Extraction failed: {result.stderr}")
    
    print()


def test_variable_substitution():
    """Test WDL to bash variable substitution patterns"""
    wdl_content = '''
    task SubstitutionTask {
      input {
        String input_name
        File data_file
      }
      command <<<
        # Test different variable patterns
        echo ${input_name}
        cat ${data_file}
        # Shell variable (should be preserved)
        echo $${USER}
        echo "Home: $${HOME}"
      >>>
      output {
        String log = stdout()
      }
    }
    '''
    
    print("Test 3: Variable Substitution")
    print("=" * 50)
    
    with tempfile.TemporaryDirectory() as tmpdir:
        wdl_file = Path(tmpdir) / "substitute.wdl"
        wdl_file.write_text(wdl_content)
        
        result = subprocess.run(
            ["python3", "wdl_to_bash_extractor.py", str(wdl_file),
             "-o", f"{tmpdir}/output", "-l", "DEBUG"],
            capture_output=True,
            text=True
        )
        
        if result.returncode == 0:
            script = Path(f"{tmpdir}/output/SubstitutionTask.sh")
            if script.exists():
                print("✓ Script generated successfully")
                content = script.read_text()
                
                # Check substitutions
                checks = [
                    ("input_name substituted", "$input_name" in content and "${input_name}" not in content),
                    ("data_file substituted", "$data_file" in content),
                    ("Shell variables preserved", "${USER}" in content and "${HOME}" in content),
                ]
                
                for check_name, passed in checks:
                    print(f"  {'✓' if passed else '✗'} {check_name}")
                
                # Show variable substitutions from log
                if "Substituted" in result.stderr:
                    print("\n  Variable substitutions logged:")
                    for line in result.stderr.split('\n'):
                        if "Substituted" in line and "->" in line:
                            print(f"    {line.split('Substituted:')[1]}")
            else:
                print("✗ Script not created")
        else:
            print(f"✗ Extraction failed: {result.stderr}")
    
    print()


def test_bash_validation():
    """Test bash syntax validation"""
    wdl_content = '''
    task ValidTask {
      input {
        String filename
      }
      command <<<
        set -euo pipefail
        
        for f in *.txt; do
          echo "Processing $f"
          grep "pattern" "$f" > "output_$f"
        done
        
        ls -la output_*
      >>>
      output {
        Array[File] outputs = glob("output_*")
      }
    }
    '''
    
    print("Test 4: Bash Syntax Validation")
    print("=" * 50)
    
    with tempfile.TemporaryDirectory() as tmpdir:
        wdl_file = Path(tmpdir) / "valid.wdl"
        wdl_file.write_text(wdl_content)
        
        result = subprocess.run(
            ["python3", "wdl_to_bash_extractor.py", str(wdl_file),
             "-o", f"{tmpdir}/output"],
            capture_output=True,
            text=True
        )
        
        if result.returncode == 0:
            script = Path(f"{tmpdir}/output/ValidTask.sh")
            if script.exists():
                print("✓ Script generated successfully")
                
                # Manually validate with bash -n
                validation = subprocess.run(
                    ["bash", "-n", str(script)],
                    capture_output=True,
                    text=True
                )
                
                if validation.returncode == 0:
                    print("✓ Bash syntax validation passed")
                else:
                    print(f"✗ Bash syntax error: {validation.stderr}")
            else:
                print("✗ Script not created")
        else:
            print(f"✗ Extraction failed: {result.stderr}")
    
    print()


def test_input_json_generation():
    """Test generation of input JSON templates"""
    wdl_content = '''
    task JSONTask {
      input {
        String sample_name
        File bam_file
        Int threads = 4
        Boolean run_qc
        Float? contamination_threshold
        Array[String] analysis_types
      }
      command {
        echo "Analyzing ${sample_name}"
      }
      output {
        String result = "done"
      }
    }
    '''
    
    print("Test 5: Input JSON Template Generation")
    print("=" * 50)
    
    with tempfile.TemporaryDirectory() as tmpdir:
        wdl_file = Path(tmpdir) / "json_test.wdl"
        wdl_file.write_text(wdl_content)
        
        result = subprocess.run(
            ["python3", "wdl_to_bash_extractor.py", str(wdl_file),
             "-o", f"{tmpdir}/output"],
            capture_output=True,
            text=True
        )
        
        if result.returncode == 0:
            json_file = Path(f"{tmpdir}/output/JSONTask_inputs.json")
            if json_file.exists():
                print("✓ JSON template generated successfully")
                
                import json
                try:
                    data = json.loads(json_file.read_text())
                    print("\n  Generated template:")
                    print("  {")
                    for key, value in data.items():
                        print(f'    "{key}": {json.dumps(value)},')
                    print("  }")
                    
                    checks = [
                        ("String fields", "sample_name" in data or "sample name" in data),
                        ("File fields", "bam_file" in data),
                        ("Array fields", "analysis_types" in data),
                    ]
                    
                    print("\n  Content checks:")
                    for check_name, passed in checks:
                        print(f"    {'✓' if passed else '✗'} {check_name}")
                        
                except json.JSONDecodeError as e:
                    print(f"✗ Invalid JSON generated: {e}")
            else:
                print("✗ JSON template not created")
        else:
            print(f"✗ Extraction failed: {result.stderr}")
    
    print()


def run_all_tests():
    """Run all tests"""
    print("\n")
    print("=" * 60)
    print("WDL to Bash Extractor - Test Suite")
    print("=" * 60)
    print()
    
    tests = [
        ("Simple Task", test_simple_task),
        ("Optional Inputs", test_optional_inputs),
        ("Variable Substitution", test_variable_substitution),
        ("Bash Validation", test_bash_validation),
        ("JSON Generation", test_input_json_generation),
    ]
    
    passed = 0
    failed = 0
    
    for test_name, test_func in tests:
        try:
            test_func()
            passed += 1
        except Exception as e:
            print(f"✗ {test_name} failed with exception: {e}\n")
            failed += 1
    
    print("=" * 60)
    print(f"Test Summary: {passed} passed, {failed} failed")
    print("=" * 60)


if __name__ == "__main__":
    import sys
    
    # Change to script directory
    os.chdir(Path(__file__).parent)
    
    run_all_tests()
