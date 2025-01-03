import os
import re
from collections import defaultdict

def find_python_files(root_dir):
    """Find all Python files in the project, excluding venv and __pycache__ directories."""
    python_files = []
    for dirpath, _, filenames in os.walk(root_dir):
        if "venv" in dirpath or "__pycache__" in dirpath:
            continue
        for file in filenames:
            if file.endswith(".py"):
                python_files.append(os.path.relpath(os.path.join(dirpath, file), root_dir))
    return python_files


def extract_imports(file_path):
    """Extract imports from a Python file."""
    imports = set()
    with open(file_path, "r", encoding="utf-8") as file:
        for line in file:
            # Match imports
            match = re.match(r"^\s*(?:from\s+([\w\.]+)|import\s+([\w\.]+))", line)
            if match:
                module = match.group(1) or match.group(2)
                imports.add(module.split(".")[0])  # Take only the top-level module
    return imports


def map_dependencies(root_dir, python_files):
    """Map dependencies between Python files."""
    dependencies = defaultdict(set)
    all_modules = {os.path.splitext(file)[0].replace("/", ".") for file in python_files}

    for file in python_files:
        file_path = os.path.join(root_dir, file)
        imports = extract_imports(file_path)
        file_module = os.path.splitext(file)[0].replace("/", ".")
        for module in imports:
            if module in all_modules:
                dependencies[file_module].add(module)
    return dependencies


def find_unused_files(dependencies, all_files):
    """Find Python files that are not imported by any other file."""
    imported_files = {dep for deps in dependencies.values() for dep in deps}
    unused_files = [file for file in all_files if os.path.splitext(file)[0].replace("/", ".") not in imported_files]
    return unused_files


if __name__ == "__main__":
    root_dir = "."  # Replace with your project root directory
    python_files = find_python_files(root_dir)
    dependencies = map_dependencies(root_dir, python_files)

    # Save dependencies to file
    with open("dependencies.txt", "w") as dep_file:
        for file, deps in dependencies.items():
            dep_file.write(f"{file} imports:\n")
            for dep in deps:
                dep_file.write(f"  - {dep}\n")
            dep_file.write("\n")

    # Find and save unused files
    unused_files = find_unused_files(dependencies, python_files)
    with open("unused_files.txt", "w") as unused_file:
        unused_file.writelines(f"{file}\n" for file in unused_files)

    print("Dependency mapping saved to 'dependencies.txt'.")
    print("Unused Python files saved to 'unused_files.txt'.")
