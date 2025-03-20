#!/usr/bin/env python3
"""
Check environment setup for mdCATH processing.
"""

import os
import sys
import subprocess
import shutil
from colorama import init, Fore, Style

# Initialize colorama
init()

def print_status(component, status, message=None):
    """Print status with color coding."""
    if status == "OK":
        status_str = f"{Fore.GREEN}[OK]{Style.RESET_ALL}"
    elif status == "WARNING":
        status_str = f"{Fore.YELLOW}[WARNING]{Style.RESET_ALL}"
    else:
        status_str = f"{Fore.RED}[FAILED]{Style.RESET_ALL}"
    
    print(f"{status_str} {component}")
    if message:
        print(f"     {message}")

def check_python_version():
    """Check Python version."""
    version = sys.version_info
    if version.major == 3 and version.minor >= 8:
        print_status("Python version", "OK", f"Python {version.major}.{version.minor}.{version.micro}")
    else:
        print_status("Python version", "FAILED", 
                     f"Python {version.major}.{version.minor}.{version.micro} (Python 3.8+ recommended)")

def check_package(package_name):
    """Check if a Python package is installed."""
    try:
        __import__(package_name)
        print_status(f"Package: {package_name}", "OK")
        return True
    except ImportError:
        print_status(f"Package: {package_name}", "FAILED", "Not installed")
        return False

def check_executable(name, command=None):
    """Check if an executable is available."""
    if command is None:
        command = name
    
    path = shutil.which(command)
    if path:
        print_status(f"Executable: {name}", "OK", f"Found at: {path}")
        return True
    else:
        print_status(f"Executable: {name}", "FAILED", "Not found in PATH")
        return False

def check_directory(path, label=None):
    """Check if a directory exists."""
    if label is None:
        label = path
    
    if os.path.exists(path):
        if os.path.isdir(path):
            print_status(f"Directory: {label}", "OK", f"Found at: {path}")
            return True
        else:
            print_status(f"Directory: {label}", "FAILED", f"Path exists but is not a directory: {path}")
            return False
    else:
        print_status(f"Directory: {label}", "FAILED", f"Not found: {path}")
        return False

def check_file(path, label=None):
    """Check if a file exists."""
    if label is None:
        label = path
    
    if os.path.exists(path):
        if os.path.isfile(path):
            print_status(f"File: {label}", "OK", f"Found at: {path}")
            return True
        else:
            print_status(f"File: {label}", "FAILED", f"Path exists but is not a file: {path}")
            return False
    else:
        print_status(f"File: {label}", "FAILED", f"Not found: {path}")
        return False

def check_h5_files(directory, count=5):
    """Check if H5 files exist in the directory."""
    import glob
    h5_files = glob.glob(os.path.join(directory, "*.h5"))
    
    if h5_files:
        print_status("H5 files", "OK", f"Found {len(h5_files)} files in {directory}")
        print(f"    First {min(count, len(h5_files))} files:")
        for i, f in enumerate(h5_files[:count]):
            print(f"     - {os.path.basename(f)}")
        return True
    else:
        print_status("H5 files", "FAILED", f"No .h5 files found in {directory}")
        return False

def check_msms_executables(msms_dir):
    """Check MSMS executables."""
    executables = [
        "msms.x86_64Linux2.2.6.1",
        "pdb_to_xyzr"
    ]
    
    all_exist = True
    for exe in executables:
        path = os.path.join(msms_dir, exe)
        if os.path.exists(path):
            if os.access(path, os.X_OK):
                print_status(f"MSMS executable: {exe}", "OK")
            else:
                print_status(f"MSMS executable: {exe}", "WARNING", "File exists but is not executable")
                print(f"     Run: chmod +x {path}")
                all_exist = False
        else:
            print_status(f"MSMS executable: {exe}", "FAILED", f"Not found: {path}")
            all_exist = False
    
    return all_exist

def main():
    """Main function."""
    print(f"{Fore.CYAN}=== mdCATH Processing Environment Check ==={Style.RESET_ALL}")
    print()
    
    # Check Python version
    check_python_version()
    print()
    
    # Check required Python packages
    print(f"{Fore.CYAN}Checking Python packages:{Style.RESET_ALL}")
    required_packages = [
        "h5py", "numpy", "pandas", "matplotlib", "seaborn", 
        "yaml", "biopython", "tqdm", "colorama"
    ]
    
    for package in required_packages:
        check_package(package)
    print()
    
    # Check executables
    print(f"{Fore.CYAN}Checking executables:{Style.RESET_ALL}")
    check_executable("make-frame-dataset")
    print()
    
    # Check MSMS executables
    print(f"{Fore.CYAN}Checking MSMS executables:{Style.RESET_ALL}")
    msms_dir = "./msms_executables"
    check_msms_executables(msms_dir)
    print()
    
    # Check configuration
    print(f"{Fore.CYAN}Checking configuration:{Style.RESET_ALL}")
    check_file("./src/mdcath/config/default_config.yaml", "Default config")
    print()
    
    # Check mdCATH dataset
    print(f"{Fore.CYAN}Checking mdCATH dataset:{Style.RESET_ALL}")
    mdcath_dir = "/mnt/datasets/MD_CATH/data"
    if check_directory(mdcath_dir, "mdCATH data directory"):
        check_h5_files(mdcath_dir)
    print()
    
    # Check output directories
    print(f"{Fore.CYAN}Checking output directories:{Style.RESET_ALL}")
    output_dir = "./outputs"
    if not check_directory(output_dir, "Output directory"):
        print(f"     Will be created during execution")
    print()
    
    print(f"{Fore.CYAN}=== Environment Check Complete ==={Style.RESET_ALL}")

if __name__ == "__main__":
    main()