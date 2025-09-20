#!/usr/bin/env python3
import re
from collections import defaultdict
import json
import sys

def parse_bcheck_output(filename):
    """Parse bcheck output and organize by source file"""

    file_to_functions = defaultdict(list)
    current_function = None
    current_errors = []

    with open(filename, 'r') as f:
        lines = f.readlines()

    i = 0
    while i < len(lines):
        line = lines[i].strip()

        # Skip empty lines and known noise
        if not line or line.startswith('ERROR: too many states'):
            i += 1
            continue

        # New function
        if line.startswith('Function '):
            # Save previous function if exists
            if current_function and current_errors:
                # Extract source file from error lines
                source_file = None
                for error in current_errors:
                    match = re.search(r'/([^/]+\.cpp):\d+$', error)
                    if match:
                        source_file = match.group(1)
                        break

                if source_file:
                    file_to_functions[source_file].append({
                        'function': current_function,
                        'errors': current_errors
                    })

            current_function = line[9:]  # Remove "Function "
            current_errors = []

        # Error line
        elif line.startswith('['):
            current_errors.append(line)

        i += 1

    # Don't forget the last function
    if current_function and current_errors:
        source_file = None
        for error in current_errors:
            match = re.search(r'/([^/]+\.cpp):\d+$', error)
            if match:
                source_file = match.group(1)
                break
        if source_file:
            file_to_functions[source_file].append({
                'function': current_function,
                'errors': current_errors
            })

    return file_to_functions

def filter_relevant_functions(file_to_functions):
    """Filter out Rcpp and other library functions"""
    filtered = defaultdict(list)

    for source_file, functions in file_to_functions.items():
        # Skip Rcpp headers
        if 'Rcpp/include' in source_file:
            continue

        for func_info in functions:
            func_name = func_info['function']
            # Skip Rcpp internal functions
            if func_name.startswith('Rcpp::'):
                continue
            # Skip RcppExports functions (auto-generated)
            if '_gflow_' in func_name:
                continue

            filtered[source_file].append(func_info)

    return filtered

def generate_reports(file_to_functions, output_dir='rchk_reports'):
    """Generate individual reports for each source file"""
    import os
    os.makedirs(output_dir, exist_ok=True)

    # Summary report
    summary = []

    for source_file, functions in sorted(file_to_functions.items()):
        if not functions:
            continue

        # Create report file
        base_name = source_file.replace('.cpp', '')
        report_file = f"{output_dir}/{base_name}_rchk_errors.txt"

        with open(report_file, 'w') as f:
            f.write(f"RCHK Error Report for {source_file}\n")
            f.write("=" * 60 + "\n\n")

            for func_info in functions:
                f.write(f"Function {func_info['function']}\n")
                for error in func_info['errors']:
                    f.write(f"  {error}\n")
                f.write("\n")

        # Add to summary
        func_names = [f['function'] for f in functions]
        summary.append(f"{source_file}: {', '.join(func_names)}")

    # Write summary
    with open(f"{output_dir}/summary.txt", 'w') as f:
        f.write("RCHK Error Summary\n")
        f.write("=" * 60 + "\n\n")
        for line in summary:
            f.write(f"{line}\n")

    # Also create JSON for programmatic access
    with open(f"{output_dir}/errors.json", 'w') as f:
        json.dump(dict(file_to_functions), f, indent=2)

if __name__ == "__main__":
    if len(sys.argv) != 2:
        print("Usage: python3 process_rchk.py <bcheck_file>")
        sys.exit(1)

    bcheck_file = sys.argv[1]
    file_to_functions = parse_bcheck_output(bcheck_file)
    filtered = filter_relevant_functions(file_to_functions)
    generate_reports(filtered)

    print(f"Generated reports in rchk_reports/")
    print(f"Found errors in {len(filtered)} source files")
