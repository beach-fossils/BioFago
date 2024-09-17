import os
import zipfile
import sys
from pathlib import Path


def zipdir(path, ziph):
    path = Path(path)
    for root, _, files in os.walk(path):
        for file in files:
            file_path = Path(root) / file
            archive_path = file_path.relative_to(path)
            ziph.write(file_path, archive_path)


def create_zip(source_dir, output_file):
    source_dir = Path(source_dir)
    output_file = Path(output_file)

    if not source_dir.exists():
        print(f"Error: Source directory '{source_dir}' does not exist.")
        return False

    try:
        with zipfile.ZipFile(output_file, 'w', zipfile.ZIP_DEFLATED) as zipf:
            zipdir(source_dir, zipf)
        print(f"Zip file created successfully: {output_file}")
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