"""
WDL parser module for generating Mermaid diagrams from WDL workflows.

Uses miniwdl's WDL.Tree to parse WDL files and extract task dependencies,
then generates Mermaid diagram representations of the workflow DAG.
"""

import sys
from pathlib import Path
from typing import Dict, Set, List, Tuple, Optional
import json

try:
    import WDL
except ImportError:
    WDL = None


def load_wdl(wdl_path: str) -> Optional[WDL.Workflow]:
    """
    Load and parse a WDL file using miniwdl's WDL.Tree.
    
    Args:
        wdl_path: Path to the WDL file
    
    Returns:
        WDL workflow object or None if parsing fails
    
    Raises:
        ImportError: If miniwdl is not installed
    """
    if WDL is None:
        raise ImportError(
            "miniwdl is not installed. Install it with: pip install miniwdl"
        )
    
    try:
        wdl_file = Path(wdl_path).resolve()
        
        # Read the WDL file content
        with open(wdl_file, 'r') as f:
            wdl_content = f.read()
        
        # Parse the WDL content (not the file path)
        doc = WDL.parse_document(wdl_content, uri=str(wdl_file))
        
        # Try to typecheck, but don't fail if imports or dependencies are missing
        try:
            doc.typecheck()
        except (WDL.Error.MultipleValidationErrors, WDL.Error.EvalError, Exception) as e:
            # If imports or task dependencies are not found, we can still try to analyze
            # the workflow structure. Just skip typecheck.
            pass
        
        # Get the workflow from the document
        workflow = doc.workflow
        
        if workflow is None:
            print(f"Error: No workflow found in {wdl_path}", file=sys.stderr)
            return None
        
        return workflow
    
    except WDL.Error.SyntaxError as e:
        print(f"Syntax error in WDL file: {e}", file=sys.stderr)
        return None
    except WDL.Error.ValidationError as e:
        print(f"Validation error in WDL file: {e}", file=sys.stderr)
        return None
    except FileNotFoundError:
        print(f"Error: WDL file not found: {wdl_path}", file=sys.stderr)
        return None
    except Exception as e:
        print(f"Error parsing WDL file: {type(e).__name__}: {e}", file=sys.stderr)
        import traceback
        traceback.print_exc()
        return None


def extract_calls_and_dependencies(workflow: WDL.Workflow) -> Tuple[Dict[str, str], List[Tuple[str, str]]]:
    """
    Extract all calls and their dependencies from a workflow.
    
    Args:
        workflow: WDL workflow object
    
    Returns:
        Tuple of (call_map, dependencies) where:
        - call_map: Dict mapping call names to task names
        - dependencies: List of (source_call, target_call) tuples representing edges
    """
    call_map = {}  # Maps call name to task name
    calls_by_id = {}  # Maps workflow_node_id to call name
    dependencies = []  # List of (from_call, to_call) tuples
    
    def extract_calls_recursive(body_element):
        """Recursively extract calls from workflow body elements."""
        if isinstance(body_element, WDL.Call):
            call_name = body_element.name
            
            # Handle cases where callee might be None (e.g., when imports are missing)
            if body_element.callee is not None:
                task_name = body_element.callee.name
            else:
                # Try to get the task name from the callee_id if available
                if hasattr(body_element, 'callee_id'):
                    callee_id = body_element.callee_id
                    # callee_id might be a list like ['Tasks', 'TaskName'] or a string
                    if isinstance(callee_id, list):
                        task_name = callee_id[-1]  # Get the last element
                    else:
                        task_name = str(callee_id)
                else:
                    task_name = call_name
            
            call_map[call_name] = task_name
            calls_by_id[body_element.workflow_node_id] = call_name
        
        elif isinstance(body_element, (WDL.Scatter, WDL.Conditional)):
            # Recursively process scatter and conditional bodies
            if hasattr(body_element, 'body'):
                if isinstance(body_element.body, list):
                    for elem in body_element.body:
                        extract_calls_recursive(elem)
                else:
                    # body might be a single element
                    extract_calls_recursive(body_element.body)
    
    # First pass: extract all calls recursively
    if hasattr(workflow, 'body'):
        for element in workflow.body:
            extract_calls_recursive(element)
    
    # Second pass: extract dependencies using miniwdl's node_id system
    if hasattr(workflow, 'body'):
        def collect_calls_recursive(body_element, all_calls=None):
            """Collect all call objects recursively."""
            if all_calls is None:
                all_calls = []
            
            if isinstance(body_element, WDL.Call):
                all_calls.append(body_element)
            elif isinstance(body_element, (WDL.Scatter, WDL.Conditional)):
                if hasattr(body_element, 'body'):
                    if isinstance(body_element.body, list):
                        for elem in body_element.body:
                            collect_calls_recursive(elem, all_calls)
                    else:
                        collect_calls_recursive(body_element.body, all_calls)
            
            return all_calls
        
        # Collect all calls
        all_calls = []
        for element in workflow.body:
            all_calls = collect_calls_recursive(element, all_calls)
        
        # Extract dependencies from each call
        for call_obj in all_calls:
            call_name = call_obj.name
            
            # Get dependencies for this call
            for dep_id in call_obj.workflow_node_dependencies:
                # Filter for call dependencies (not declarations)
                if dep_id.startswith('call-'):
                    # Extract the call name from the dependency id
                    dep_call_name = dep_id.split('-', 1)[1]
                    
                    # Add edge from dependency to this call
                    if dep_call_name in call_map and dep_call_name != call_name:
                        dep_pair = (dep_call_name, call_name)
                        if dep_pair not in dependencies:
                            dependencies.append(dep_pair)
    
    return call_map, dependencies


def generate_mermaid_diagram(
    call_map: Dict[str, str],
    dependencies: List[Tuple[str, str]],
    wdl_path: str
) -> str:
    """
    Generate a Mermaid diagram from calls and dependencies.
    
    Args:
        call_map: Dict mapping call aliases to task names
        dependencies: List of (source_call, target_call) tuples
        wdl_path: Path to the WDL file (for reference)
    
    Returns:
        Mermaid diagram as a string
    """
    diagram_lines = []
    
    # Add metadata
    wdl_name = Path(wdl_path).name
    num_calls = len(call_map)
    
    diagram_lines.append("```mermaid")
    diagram_lines.append("graph LR")
    
    # Add node definitions
    for call_alias, task_name in sorted(call_map.items()):
        if call_alias != task_name:
            # Show alias and task name for different aliases
            node_label = f"{call_alias}\\n({task_name})"
        else:
            # Just show the call name if it's the same as task name
            node_label = call_alias
        
        diagram_lines.append(f"    {call_alias}[\"{node_label}\"]")
    
    # Add edges (dependencies)
    if dependencies:
        # Remove duplicates and sort for consistent output
        unique_deps = sorted(set(dependencies))
        for source, target in unique_deps:
            diagram_lines.append(f"    {source} --> {target}")
    
    diagram_lines.append("```")
    
    return "\n".join(diagram_lines)


def generate_summary(
    wdl_path: str,
    workflow: WDL.Workflow,
    call_map: Dict[str, str]
) -> str:
    """
    Generate a summary of the workflow structure.
    
    Args:
        wdl_path: Path to the WDL file
        workflow: WDL workflow object
        call_map: Dict mapping call aliases to task names
    
    Returns:
        Summary text
    """
    summary_lines = []
    wdl_name = Path(wdl_path).name
    
    summary_lines.append(f"WDL: {wdl_name}")
    
    # Count unique tasks referenced
    unique_tasks = len(set(call_map.values()))
    summary_lines.append(f"Tasks: {unique_tasks}")
    summary_lines.append(f"Calls: {len(call_map)}")
    summary_lines.append("")
    
    # List available tasks
    if call_map:
        summary_lines.append("Available tasks:")
        for task_name in sorted(set(call_map.values())):
            summary_lines.append(f"  - {task_name}")
    
    return "\n".join(summary_lines)


def generate_diagram(wdl_path: str, output_file: Optional[str] = None) -> Optional[str]:
    """
    Generate a Mermaid diagram from a WDL workflow.
    
    Args:
        wdl_path: Path to the WDL file
        output_file: Optional output file path for the diagram
    
    Returns:
        The generated diagram as a string, or None if generation failed
    
    Example:
        >>> diagram = generate_diagram("/path/to/workflow.wdl")
        >>> print(diagram)
        
        >>> generate_diagram("/path/to/workflow.wdl", "output.mmd")
    """
    # Check if miniwdl is available
    if WDL is None:
        print("Error: miniwdl is not installed. Install it with: pip install miniwdl", file=sys.stderr)
        return None
    
    # Load and parse the WDL file
    workflow = load_wdl(wdl_path)
    if workflow is None:
        return None
    
    # Extract calls and dependencies
    call_map, dependencies = extract_calls_and_dependencies(workflow)
    
    if not call_map:
        print(f"Warning: No calls found in workflow {wdl_path}", file=sys.stderr)
    
    # Generate the diagram
    diagram = generate_mermaid_diagram(call_map, dependencies, wdl_path)
    
    # Generate summary
    summary = generate_summary(wdl_path, workflow, call_map)
    
    # Combine summary and diagram
    output = f"{summary}\n\n{diagram}"
    
    # Write to file if specified
    if output_file:
        try:
            output_path = Path(output_file)
            output_path.parent.mkdir(parents=True, exist_ok=True)
            output_path.write_text(output)
            print(f"Diagram saved to: {output_file}")
        except Exception as e:
            print(f"Error writing output file: {e}", file=sys.stderr)
            return None
    
    return output


if __name__ == "__main__":
    if len(sys.argv) < 2:
        print("Usage: python wdl_parser.py <wdl_file> [output_file]", file=sys.stderr)
        sys.exit(1)
    
    wdl_path = sys.argv[1]
    output_file = sys.argv[2] if len(sys.argv) > 2 else None
    
    diagram = generate_diagram(wdl_path, output_file)
    if diagram and output_file is None:
        print(diagram)
    
    sys.exit(0 if diagram else 1)
