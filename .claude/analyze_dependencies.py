#!/usr/bin/env python3

import os
import re
import csv
from pathlib import Path

def extract_function_definitions(content):
    """Extract function names defined in R file"""
    # Match patterns like: function_name <- function(...) or function_name = function(...)
    pattern = r'^([a-zA-Z][a-zA-Z0-9._]*)\s*(?:<-|=)\s*function\s*\('
    matches = re.findall(pattern, content, re.MULTILINE)
    return set(matches)

def extract_function_calls(content):
    """Extract function calls from R file"""
    # Remove comments first
    lines = content.split('\n')
    cleaned_lines = []
    for line in lines:
        # Remove comments (everything after #)
        if '#' in line:
            line = line[:line.index('#')]
        cleaned_lines.append(line)
    content = '\n'.join(cleaned_lines)
    
    # Remove string literals to avoid false positives
    content = re.sub(r'"[^"]*"', '""', content)
    content = re.sub(r"'[^']*'", "''", content)
    
    # Match function calls - word followed by (
    pattern = r'\b([a-zA-Z][a-zA-Z0-9._]*)\s*\('
    matches = re.findall(pattern, content)
    
    # Filter out R keywords and common functions
    r_keywords = {
        'if', 'else', 'for', 'while', 'repeat', 'function', 'return', 'next', 'break',
        'switch', 'tryCatch', 'try', 'stop', 'warning', 'message', 'cat', 'print',
        'invisible', 'suppressWarnings', 'suppressMessages', 'library', 'require',
        'source', 'list', 'c', 'seq', 'rep', 'length', 'nrow', 'ncol', 'dim',
        'matrix', 'array', 'data.frame', 'as.matrix', 'as.data.frame', 'as.numeric',
        'as.character', 'as.logical', 'as.integer', 'is.null', 'is.na', 'is.numeric',
        'is.character', 'is.logical', 'any', 'all', 'which', 'sum', 'mean', 'min', 'max',
        'paste', 'paste0', 'sprintf', 'grep', 'gsub', 'sub', 'str', 'eval', 'parse',
        'deparse', 'substitute', 'quote', 'expression', 'call', 'do.call', 'lapply',
        'sapply', 'apply', 'tapply', 'mapply', 'vapply', 'replicate', 'unique',
        'sort', 'order', 'rank', 'table', 'unlist', 'names', 'colnames', 'rownames',
        'setNames', 'match', 'merge', 'rbind', 'cbind', 'split', 'aggregate',
        'with', 'within', 'transform', 'subset', 'head', 'tail', 'sample', 'set.seed',
        'runif', 'rnorm', 'rpois', 'rbinom', 'dnorm', 'pnorm', 'qnorm', 'lm', 'glm',
        'predict', 'coef', 'summary', 'plot', 'lines', 'points', 'abline', 'par',
        'dev.off', 'pdf', 'png', 'jpeg', 'svg', 'read.csv', 'write.csv', 'read.table',
        'write.table', 'readLines', 'writeLines', 'file', 'dir', 'getwd', 'setwd',
        'Sys.time', 'proc.time', 'system.time', 'gc', 'rm', 'ls', 'exists', 'get',
        'assign', 'attach', 'detach', 'search', 'options', 'getOption', 'setOption',
        'class', 'typeof', 'mode', 'attributes', 'attr', 'structure', 'unclass',
        'inherits', 'methods', 'UseMethod', 'NextMethod', 'standardGeneric', 'floor',
        'ceiling', 'round', 'abs', 'sqrt', 'exp', 'log', 'log10', 'sin', 'cos', 'tan',
        'asin', 'acos', 'atan', 'sinh', 'cosh', 'tanh', 'var', 'sd', 'cor', 'cov',
        'median', 'quantile', 'range', 'diff', 'cumsum', 'cumprod', 'cummin', 'cummax',
        'duplicated', 'rev', 'append', 'intersect', 'union', 'setdiff', 'setequal',
        '%in%', 'match.arg', 'missing', 'on.exit', 'sys.call', 'sys.frame', 'parent.frame',
        't', 'crossprod', 'tcrossprod', 'solve', 'chol', 'qr', 'svd', 'eigen', 'det',
        'diag', 'upper.tri', 'lower.tri', 'I', 'identity', 'identical', 'all.equal',
        'isTRUE', 'isFALSE', 'xor', 'bitwAnd', 'bitwOr', 'bitwXor', 'bitwNot', 'intToBits',
        'packBits', 'strsplit', 'nchar', 'tolower', 'toupper', 'chartr', 'abbreviate',
        'make.names', 'format', 'formatC', 'prettyNum', 'strftime', 'strptime', 'as.Date',
        'as.POSIXct', 'as.POSIXlt', 'weekdays', 'months', 'quarters', 'julian', 'date',
        'Sys.Date', 'factor', 'levels', 'nlevels', 'reorder', 'relevel', 'cut', 'pretty',
        'findInterval', 'approx', 'approxfun', 'spline', 'splinefun', 'loess', 'lowess',
        'smooth', 'fft', 'mvfft', 'filter', 'convolve', 'spectrum', 'acf', 'pacf', 'ccf',
        'arima', 'Box.test', 'PP.test', 'adf.test', 't.test', 'wilcox.test', 'ks.test',
        'shapiro.test', 'chisq.test', 'fisher.test', 'binom.test', 'prop.test', 'mcnemar.test',
        'dist', 'hclust', 'kmeans', 'prcomp', 'princomp', 'cancor', 'lda', 'qda', 'naive_bayes',
        'double', 'integer', 'complex', 'character', 'logical', 'raw', 'NULL', 'NA', 'NaN',
        'Inf', 'TRUE', 'FALSE', 'T', 'F', 'pi', 'letters', 'LETTERS', 'month.abb', 'month.name',
        'body', 'formals', 'environment', 'environmentName', 'env.profile', 'sys.parent',
        'sys.nframe', 'sys.calls', 'sys.frames', 'sys.parents', 'sys.on.exit', 'sys.status',
        'globalenv', 'baseenv', 'emptyenv', 'parent.env', 'is.na<-', 'length<-', 'levels<-',
        'names<-', 'dimnames<-', 'dim<-', 'class<-', 'attr<-', 'attributes<-', 'comment<-',
        'row.names<-', 'rownames<-', 'colnames<-', 'body<-', 'environment<-', 'formals<-',
        'args', 'commandArgs', 'do.call', 'match.call', 'match.fun', 'browser', 'debug',
        'undebug', 'isdebugged', 'debugonce', 'trace', 'untrace', 'traceback', 'locator'
    }
    
    # Filter out R keywords and duplicates
    function_calls = set(m for m in matches if m not in r_keywords)
    return function_calls

def analyze_file(filepath, all_msr2_functions, all_gflow_functions):
    """Analyze dependencies for a single R file"""
    with open(filepath, 'r', encoding='utf-8', errors='ignore') as f:
        content = f.read()
    
    defined_functions = extract_function_definitions(content)
    called_functions = extract_function_calls(content)
    
    # Functions from this file
    total_functions = len(defined_functions)
    
    # Find dependent functions (calls to other msr2 functions not in gflow)
    dependent_calls = called_functions & all_msr2_functions
    dependent_calls = dependent_calls - all_gflow_functions
    dependent_calls = dependent_calls - defined_functions  # Remove self-references
    
    dependent_functions = len(dependent_calls)
    
    return total_functions, dependent_functions, defined_functions, dependent_calls

def main():
    msr2_dir = Path('/Users/pgajer/current_projects/msr2/R')
    gflow_dir = Path('/Users/pgajer/current_projects/gflow/R')
    output_file = Path('/Users/pgajer/current_projects/gflow/.claude/msr2_dependencies.csv')
    
    # Read list of files only in msr2
    with open('/tmp/msr2_only_files.txt', 'r') as f:
        msr2_only_files = [line.strip() for line in f if line.strip()]
    
    print(f"Found {len(msr2_only_files)} files in msr2/R but not in gflow/R")
    
    # First pass: collect all function definitions
    print("\nFirst pass: Extracting all function definitions...")
    
    all_msr2_functions = set()
    all_gflow_functions = set()
    file_functions = {}
    
    # Get functions from msr2 files not in gflow
    for i, filename in enumerate(msr2_only_files):
        if (i + 1) % 10 == 0:
            print(f"  Processing msr2 file {i+1}/{len(msr2_only_files)}")
        
        filepath = msr2_dir / filename
        if filepath.exists() and filepath.is_file():
            with open(filepath, 'r', encoding='utf-8', errors='ignore') as f:
                content = f.read()
            functions = extract_function_definitions(content)
            file_functions[filename] = functions
            all_msr2_functions.update(functions)
    
    print(f"  Found {len(all_msr2_functions)} functions in msr2-only files")
    
    # Get functions from gflow files
    gflow_files = list(gflow_dir.glob('*.R'))
    for i, filepath in enumerate(gflow_files):
        if (i + 1) % 10 == 0:
            print(f"  Processing gflow file {i+1}/{len(gflow_files)}")
        
        with open(filepath, 'r', encoding='utf-8', errors='ignore') as f:
            content = f.read()
        functions = extract_function_definitions(content)
        all_gflow_functions.update(functions)
    
    print(f"  Found {len(all_gflow_functions)} functions in gflow files")
    
    # Second pass: analyze dependencies
    print("\nSecond pass: Analyzing dependencies...")
    
    results = []
    for i, filename in enumerate(msr2_only_files):
        if (i + 1) % 10 == 0:
            print(f"  Analyzing file {i+1}/{len(msr2_only_files)}")
        
        filepath = msr2_dir / filename
        if filepath.exists() and filepath.is_file():
            total_functions, dependent_functions, defined_funcs, dependent_calls = analyze_file(
                filepath, all_msr2_functions, all_gflow_functions
            )
            results.append({
                'filename': filename,
                'total_functions': total_functions,
                'dependent_functions': dependent_functions,
                'function_names': ', '.join(sorted(defined_funcs)) if defined_funcs else '',
                'dependent_calls': ', '.join(sorted(dependent_calls)) if dependent_calls else ''
            })
    
    # Sort results
    results.sort(key=lambda x: (x['dependent_functions'], x['total_functions']))
    
    # Save to CSV
    print(f"\nSaving results to {output_file}")
    output_file.parent.mkdir(parents=True, exist_ok=True)
    
    with open(output_file, 'w', newline='') as csvfile:
        fieldnames = ['filename', 'total_functions', 'dependent_functions', 'function_names', 'dependent_calls']
        writer = csv.DictWriter(csvfile, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(results)
    
    print(f"Analysis complete! Results saved to {output_file}")
    
    # Print summary
    print("\nSummary:")
    print(f"  Files with 0 dependencies: {sum(1 for r in results if r['dependent_functions'] == 0)}")
    print(f"  Files with 1-5 dependencies: {sum(1 for r in results if 1 <= r['dependent_functions'] <= 5)}")
    print(f"  Files with 6-10 dependencies: {sum(1 for r in results if 6 <= r['dependent_functions'] <= 10)}")
    print(f"  Files with >10 dependencies: {sum(1 for r in results if r['dependent_functions'] > 10)}")

if __name__ == '__main__':
    main()