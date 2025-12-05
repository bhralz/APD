import os
import shutil
from pathlib import Path
from typing import List
import tkinter as tk
from tkinter import filedialog
from tqdm import tqdm


class FolderMerger:
    def __init__(self):
        # Hide the root tkinter window
        self.root = tk.Tk()
        self.root.withdraw()

        self.source1: Path | None = None
        self.source2: Path | None = None
        self.destination: Path | None = None

    def select_source1(self):
        folder = filedialog.askdirectory(title="Select First Source Folder")
        if folder:
            self.source1 = Path(folder)
            print(f"Source 1: {self.source1}")

    def select_source2(self):
        folder = filedialog.askdirectory(title="Select Second Source Folder")
        if folder:
            self.source2 = Path(folder)
            print(f"Source 2: {self.source2}")

    def select_destination(self):
        folder = filedialog.askdirectory(title="Select Destination Folder")
        if folder:
            self.destination = Path(folder)
            print(f"Destination: {self.destination}")

    def ask_paths(self):
        print("Please select the folders:")
        self.select_source1()
        if not self.source1:
            print("No first source selected. Exiting.")
            return False

        self.select_source2()
        if not self.source2:
            print("No second source selected. Exiting.")
            return False

        self.select_destination()
        if not self.destination:
            print("No destination selected. Exiting.")
            return False

        return True

    @staticmethod
    def is_barcode_folder(name: str) -> bool:
        """
        Customize this function to match your barcode folder pattern.
        Common patterns: starts with digits, contains only digits + letters, etc.
        """
        return name and (name[0].isdigit() or len(name) >= 8)  # Example: starts with digit or long code

    def get_barcode_subfolders(self, root: Path) -> List[Path]:
        if not root.exists():
            return []
        return [
            p for p in root.iterdir()
            if p.is_dir() and self.is_barcode_folder(p.name)
        ]

    def merge_folders(self):
        if not all([self.source1, self.source2, self.destination]):
            print("Paths not fully set!")
            return

        # Ensure destination exists
        self.destination.mkdir(parents=True, exist_ok=True)

        # Collect all unique barcode folders from both sources
        folders1 = {f.name: f for f in self.get_barcode_subfolders(self.source1)}
        folders2 = {f.name: f for f in self.get_barcode_subfolders(self.source2)}

        all_barcode_names = set(folders1.keys()) | set(folders2.keys())
        print(f"Found {len(all_barcode_names)} unique barcode folders to merge.")

        total_files = 0
        for name in all_barcode_names:
            src1 = folders1.get(name)
            src2 = folders2.get(name)
            dest = self.destination / name
            dest.mkdir(exist_ok=True)

            files_to_copy = []
            seen_files = set()

            # From source 1
            if src1:
                for file in src1.iterdir():
                    if file.is_file():
                        files_to_copy.append((file, dest / file.name, src1))
                        seen_files.add(file.name)

            # From source 2 (only if not already in source1, or newer)
            if src2:
                for file in src2.iterdir():
                    if file.is_file():
                        dest_path = dest / file.name
                        if file.name not in seen_files:
                            files_to_copy.append((file, dest_path, src2))
                        else:
                            # Conflict: prefer newer file
                            existing_src = src1 / file.name if src1 else None
                            if existing_src and existing_src.exists():
                                if file.stat().st_mtime > existing_src.stat().st_mtime:
                                    files_to_copy = [f for f in files_to_copy if f[1].name != file.name]
                                    files_to_copy.append((file, dest_path, src2))
                            else:
                                # If no src1 version, take src2
                                files_to_copy.append((file, dest_path, src2))

            total_files += len(files_to_copy)

        print(f"Starting to copy {total_files} files...")

        # Copy with progress bar
        with tqdm(total=total_files, desc="Merging Files", unit="file") as pbar:
            for src_file, dest_file, source_root in (tc for tc in self.get_files_to_copy_generator(all_barcode_names)):
                try:
                    shutil.copy2(src_file, dest_file)
                    pbar.set_postfix({"current": src_file.name[:30]})
                except Exception as e:
                    print(f"\nError copying {src_file}: {e}")
                pbar.update(1)

        print(f"\nMerge completed! All files are in: {self.destination}")

    def get_files_to_copy_generator(self, all_barcode_names):
        """Generator version for memory efficiency with large datasets"""
        folders1 = {f.name: f for f in self.get_barcode_subfolders(self.source1)}
        folders2 = {f.name: f for f in self.get_barcode_subfolders(self.source2)}

        for name in all_barcode_names:
            src1 = folders1.get(name)
            src2 = folders2.get(name)
            dest = self.destination / name
            dest.mkdir(exist_ok=True)

            seen_files = set()

            # Priority: source1 first
            if src1:
                for file in src1.iterdir():
                    if file.is_file():
                        yield file, dest / file.name, src1
                        seen_files.add(file.name)

            if src2:
                for file in src2.iterdir():
                    if file.is_file() and file.name not in seen_files:
                        # Optional: check timestamp to prefer newer
                        dest_path = dest / file.name
                        existing = src1 / file.name if src1 else None
                        if existing and existing.exists():
                            if file.stat().st_mtime <= existing.stat().st_mtime:
                                continue  # skip older duplicate
                        yield file, dest_path, src2


def main():
    print("Barcode Folder Merger - Union with Intelligent Conflict Resolution")
    print("=" * 70)

    merger = FolderMerger()

    if not merger.ask_paths():
        return

    print("\nConfiguration Summary:")
    print(f"   Source 1     : {merger.source1}")
    print(f"   Source 2     : {merger.source2}")
    print(f"   Destination  : {merger.destination}")

    response = input("\nStart merging? (y/n): ").strip().lower()
    if response != 'y':
        print("Cancelled.")
        return

    merger.merge_folders()


if __name__ == "__main__":
    # Required: pip install tqdm
    main()