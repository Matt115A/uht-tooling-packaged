"""External tool validation utilities."""
import shutil
import subprocess
from typing import Dict, List, Optional, Tuple


class ToolNotFoundError(Exception):
    """Raised when a required external tool is not found."""

    pass


def check_tool_available(tool_name: str) -> Tuple[bool, Optional[str]]:
    """
    Check if a tool is available on PATH.

    Args:
        tool_name: Name of the executable to check.

    Returns:
        Tuple of (available, version_or_error).
        If available is True, version_or_error contains the version string.
        If available is False, version_or_error contains the error message.
    """
    path = shutil.which(tool_name)
    if not path:
        return False, f"'{tool_name}' not found on PATH"

    # Try to get version
    try:
        result = subprocess.run(
            [tool_name, "--version"],
            capture_output=True,
            text=True,
            timeout=5,
        )
        version = (result.stdout or result.stderr).strip().split("\n")[0]
        return True, version
    except subprocess.TimeoutExpired:
        return True, "version unknown (timeout)"
    except Exception:
        return True, "version unknown"


def validate_tools(tools: List[str], raise_on_missing: bool = True) -> Dict[str, dict]:
    """
    Validate multiple tools, optionally raising ToolNotFoundError.

    Args:
        tools: List of tool names to validate.
        raise_on_missing: If True, raise ToolNotFoundError for missing tools.

    Returns:
        Dictionary mapping tool names to their status:
        {
            "tool_name": {
                "available": bool,
                "version": str or None,
                "error": str or None
            }
        }

    Raises:
        ToolNotFoundError: If raise_on_missing is True and any tool is missing.
    """
    results: Dict[str, dict] = {}
    missing: List[str] = []

    for tool in tools:
        available, info = check_tool_available(tool)
        if available:
            results[tool] = {
                "available": True,
                "version": info,
                "error": None,
            }
        else:
            results[tool] = {
                "available": False,
                "version": None,
                "error": info,
            }
            missing.append(tool)

    if raise_on_missing and missing:
        missing_str = ", ".join(missing)
        raise ToolNotFoundError(
            f"Missing required external tool(s): {missing_str}. "
            f"Install via conda: conda install -c bioconda {' '.join(missing)}"
        )

    return results


# Tool requirements per workflow
WORKFLOW_TOOLS: Dict[str, List[str]] = {
    "mutation_caller": ["mafft"],
    "umi_hunter": ["mafft"],
    "ep_library_profile": ["minimap2", "NanoFilt"],
}


def validate_workflow_tools(workflow: str, raise_on_missing: bool = True) -> Dict[str, dict]:
    """
    Validate tools required for a specific workflow.

    Args:
        workflow: Name of the workflow (e.g., "mutation_caller", "umi_hunter").
        raise_on_missing: If True, raise ToolNotFoundError for missing tools.

    Returns:
        Dictionary mapping tool names to their status.

    Raises:
        ValueError: If workflow is not recognized.
        ToolNotFoundError: If raise_on_missing is True and any tool is missing.
    """
    if workflow not in WORKFLOW_TOOLS:
        # No external tools required for this workflow
        return {}

    tools = WORKFLOW_TOOLS[workflow]
    return validate_tools(tools, raise_on_missing=raise_on_missing)


def get_tool_requirements_message(workflow: str) -> str:
    """
    Get a human-readable message about tool requirements for a workflow.

    Args:
        workflow: Name of the workflow.

    Returns:
        A message describing required tools, or empty string if none required.
    """
    if workflow not in WORKFLOW_TOOLS:
        return ""

    tools = WORKFLOW_TOOLS[workflow]
    return (
        f"This workflow requires: {', '.join(tools)}. "
        f"Install via: conda install -c bioconda {' '.join(tools)}"
    )
