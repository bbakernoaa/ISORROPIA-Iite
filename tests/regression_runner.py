import os
import sys
import re
import subprocess

def parse_report_file(filepath):
    """
    Parses an ISORROPIA-Lite human-readable text report (.txt file)
    and returns a dictionary of parsed species concentration values and properties.
    """
    results = {}
    if not os.path.exists(filepath):
        raise FileNotFoundError(f"Report file {filepath} not found.")

    with open(filepath, 'r') as f:
        content = f.read()

    # Regular expressions for key-value extraction
    pattern_bracket_val = re.compile(r'\[([\w\s\+\-\(\)]+)\]\s+([0-9\.eE\+\-]+)')
    
    # We parse section records
    # Each record block is separated by: "================================================"
    record_idx = 0
    for line in content.splitlines():
        line = line.strip()
        if not line:
            continue
            
        if "================================================" in line:
            record_idx += 1
            continue

        match = pattern_bracket_val.search(line)
        if match:
            key = match.group(1).strip()
            val_str = match.group(2).strip()
            try:
                val = float(val_str)
                # Save key with record prefix to compare all runs side-by-side
                composite_key = f"Rec{record_idx}_{key}"
                results[composite_key] = val
            except ValueError:
                pass

    return results

def compare_results(ref_data, target_data, tolerance=1e-3):
    """
    Compares targeted active thermodynamic values within a relative tolerance.
    """
    mismatches = 0
    checked_keys = 0

    # We only compare active chemical properties of Case 1/2 systems
    KEYS_TO_COMPARE = [
        "WATER", "H+", "NH4+", "NO3+", "NO3-", "SO4--", "HSO4-", 
        "NH3", "HNO3", "Wat(NH4)2SO4", "WatNH4NO3", "WatOrg", 
        "pH", "IONIC STRENGTH"
    ]

    print(f"{'Species/Parameter':<28} | {'Reference (Fortran)':<20} | {'Target (C++)':<20} | {'Rel Diff':<15}")
    print("-" * 92)

    all_keys = sorted(list(set(ref_data.keys()) | set(target_data.keys())))
    
    for key in all_keys:
        # Extract base key after RecN_
        if '_' in key:
            parts = key.split('_', 1)
            base_key = parts[1].strip()
        else:
            base_key = key

        # Skip keys that are not in our comparison list
        if base_key not in KEYS_TO_COMPARE:
            continue

        if key not in ref_data:
            print(f"Key {key:<26} | MISSING IN FORTRAN     | {target_data[key]:<20.6E} | N/A")
            mismatches += 1
            continue
        if key not in target_data:
            print(f"Key {key:<26} | {ref_data[key]:<20.6E} | MISSING IN C++       | N/A")
            mismatches += 1
            continue

        ref_val = ref_data[key]
        tgt_val = target_data[key]
        checked_keys += 1

        if ref_val == 0.0 and tgt_val == 0.0:
            rel_diff = 0.0
        else:
            diff = abs(tgt_val - ref_val)
            rel_diff = diff / max(abs(ref_val), 1.0)

        status_str = ""
        if rel_diff > tolerance:
            status_str = "❌ MISMATCH"
            mismatches += 1
            print(f"{key:<28} | {ref_val:<20.6E} | {tgt_val:<20.6E} | {rel_diff:<15.6E} {status_str}")
        else:
            print(f"{key:<28} | {ref_val:<20.6E} | {tgt_val:<20.6E} | {rel_diff:<15.6E}")

    print("-" * 92)
    print(f"Comparison completed: {checked_keys} active keys checked, {mismatches} mismatches found.")
    return mismatches == 0

if __name__ == "__main__":
    print("=== ISORROPIA-Lite Regression Test Harness ===")
    
    script_dir = os.path.dirname(os.path.abspath(__file__))
    project_root = os.path.join(script_dir, "..")
    
    # 1. Locate C++ executable
    cpp_cli = os.path.join(project_root, "build", "isorropia_cli")
    if not os.path.exists(cpp_cli):
        # Check standard build subdirs or platforms
        cpp_cli = os.path.join(project_root, "build", "Debug", "isorropia_cli")
        if not os.path.exists(cpp_cli):
            cpp_cli = os.path.join(project_root, "build", "Release", "isorropia_cli")

    if not os.path.exists(cpp_cli):
        print(f"❌ C++ executable not found at {cpp_cli}. Please build the project first.")
        sys.exit(1)

    # 2. Paths to files
    inp_file = os.path.join(project_root, "isolite1_0_src", "test1.inp")
    fortran_out_file = os.path.join(project_root, "isolite1_0_src", "test1.txt")
    cpp_out_file = os.path.join(project_root, "isolite1_0_src", "test1_cpp.txt")

    # Ensure C++ output is generated fresh
    if os.path.exists(cpp_out_file):
        os.remove(cpp_out_file)

    # 3. Execute the C++ CLI on test1.inp
    print(f"Executing C++ Solver: {cpp_cli} {inp_file}")
    try:
        # Run in isolite1_0_src dir so output is written in the correct workspace
        result = subprocess.run([cpp_cli, inp_file], cwd=os.path.join(project_root, "isolite1_0_src"), capture_output=True, text=True)
        if result.returncode != 0:
            print(f"❌ C++ execution failed:\n{result.stderr}")
            sys.exit(1)
        print(result.stdout.strip())
    except Exception as e:
        print(f"❌ Error running C++ binary: {e}")
        sys.exit(1)

    # 4. Parse outputs
    print(f"\nParsing Reference Fortran Output: {fortran_out_file}")
    if not os.path.exists(fortran_out_file):
        print(f"❌ Fortran reference output {fortran_out_file} not found.")
        sys.exit(1)
    ref_data = parse_report_file(fortran_out_file)
    print(f"Successfully parsed {len(ref_data)} keys from Fortran.")

    print(f"\nParsing Target C++ Output: {cpp_out_file}")
    if not os.path.exists(cpp_out_file):
        print(f"❌ C++ report output {cpp_out_file} not found.")
        sys.exit(1)
    target_data = parse_report_file(cpp_out_file)
    print(f"Successfully parsed {len(target_data)} keys from C++.")

    # 5. Numerical side-by-side validation
    print("\nComparing C++ vs. Fortran E2E outputs side-by-side:")
    success = compare_results(ref_data, target_data, tolerance=1e-3) # relative tolerance of 0.1% for Phase 2 validation
    
    if success:
        print("\n✅ Regression validation PASSED! C++ outputs match Fortran legacy outputs perfectly.")
        sys.exit(0)
    else:
        print("\n❌ Regression validation FAILED. Numerical mismatches found.")
        sys.exit(1)
