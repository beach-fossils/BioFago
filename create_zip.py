import os
import zipfile
import sys
from pathlib import Path
from datetime import datetime

def should_include_file(filename):
    """Check if file should be included in the ZIP"""
    exclude_patterns = [
        '.DS_Store',
        '__MACOSX',
        '._'
    ]
    return not any(pattern in str(filename) for pattern in exclude_patterns)

def zipdir(path, ziph):
    path = Path(path)
    
    # Collect and sort all files first
    files_to_add = []
    for root, _, files in os.walk(path):
        root_path = Path(root)
        # Skip __MACOSX directories
        if '__MACOSX' in str(root_path):
            continue
            
        for file in sorted(files):  # Sort files for consistent ordering
            if not should_include_file(file):
                continue
                
            file_path = root_path / file
            archive_path = Path(root).relative_to(path.parent) / file
            files_to_add.append((file_path, archive_path))
    
    # Set consistent timestamp for all files
    fixed_date = (2024, 1, 1, 0, 0, 0)  # January 1, 2024, 00:00:00
    
    # Add files to ZIP with consistent timestamp
    for file_path, archive_path in sorted(files_to_add):  # Sort again for consistency
        info = zipfile.ZipInfo(str(archive_path), fixed_date)
        info.external_attr = 0o644 << 16  # Consistent permissions
        
        with open(file_path, 'rb') as f:
            ziph.writestr(info, f.read(), zipfile.ZIP_DEFLATED)

def create_zip(source_dir, output_file):
    source_dir = Path(source_dir)
    output_file = Path(output_file)

    if not source_dir.exists():
        print(f"Error: Source directory '{source_dir}' does not exist.")
        return False

    try:
        with zipfile.ZipFile(output_file, 'w', compression=zipfile.ZIP_DEFLATED, compresslevel=9) as zipf:
            zipdir(source_dir, zipf)
            
        print(f"Zip file created: {output_file}")
        print(f"Zip file size: {output_file.stat().st_size} bytes")
        return True
    except Exception as e:
        print(f"Error creating zip file: {str(e)}")
        return False

if __name__ == "__main__":
    if len(sys.argv) != 3:
        print("Usage: python create_zip.py <source_directory> <output_zip_file>")
        sys.exit(1)

    source_dir = sys.argv[1]
    output_file = sys.argv[2]
    success = create_zip(source_dir, output_file)
    sys.exit(0 if success else 1)
