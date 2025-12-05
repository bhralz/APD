import os
import tkinter as tk
from tkinter import filedialog, messagebox, ttk
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np

class PathogenReportAnalyzer:
    def __init__(self):
        self.root = tk.Tk()
        self.root.title("Pathogen Report Analyzer - Distilled Report")
        self.root.geometry("1100x750")
        self.root.configure(bg="#f8f9fa")

        self.folder_path = None
        self.file_list = []
        self.selected_file = None
        self.df = None
        self.filtered_df = None

        self.setup_ui()
        self.style_ui()

    def style_ui(self):
        style = ttk.Style()
        style.theme_use('clam')
        style.configure("TButton", padding=8, font=('Helvetica', 11, 'bold'))
        style.configure("Treeview", rowheight=25, font=('Helvetica', 10))
        style.configure("Treeview.Heading", font=('Helvetica', 11, 'bold'), background="#4a90e2", foreground="white")

    def setup_ui(self):
        # Header
        header = tk.Label(self.root, text="Pathogen Excel Report Analyzer", 
                          font=("Helvetica", 18, "bold"), bg="#f8f9fa", fg="#2c3e50")
        header.pack(pady=20)

        # Select Folder
        tk.Button(self.root, text="Select Folder with XLSX Files", command=self.select_folder,
                  bg="#3498db", fg="white", font=("Helvetica", 13, "bold"), height=2, width=30).pack(pady=10)

        # File List
        list_frame = tk.Frame(self.root)
        list_frame.pack(pady=10, fill=tk.X, padx=30)
        tk.Label(list_frame, text="Excel Files Found:", font=("Helvetica", 12, "bold")).pack(anchor="w")

        self.file_listbox = tk.Listbox(list_frame, height=6, font=("Courier", 11), selectbackground="#a6dcef")
        self.file_listbox_scroll = ttk.Scrollbar(list_frame, orient="vertical", command=self.file_listbox.yview)

        # FIXED: was self.filebox_scroll → now correct name
        self.file_listbox.config(yscrollcommand=self.file_listbox_scroll.set)

        self.file_listbox.pack(side=tk.LEFT, fill=tk.X, expand=True, pady=5)
        self.file_listbox_scroll.pack(side=tk.RIGHT, fill=tk.Y)

        self.file_listbox.bind("<<ListboxSelect>>", self.on_file_select)

        # Data Table
        table_frame = tk.LabelFrame(self.root, text=" Distilled Report Data Preview ", font=("Helvetica", 12, "bold"), padx=10, pady=10)
        table_frame.pack(fill=tk.BOTH, expand=True, padx=30, pady=10)

        self.tree = ttk.Treeview(table_frame, show="headings")
        vsb = ttk.Scrollbar(table_frame, orient="vertical", command=self.tree.yview)
        hsb = ttk.Scrollbar(table_frame, orient="horizontal", command=self.tree.xview)
        self.tree.configure(yscrollcommand=vsb.set, xscrollcommand=hsb.set)

        self.tree.pack(side=tk.LEFT, fill=tk.BOTH, expand=True)
        vsb.pack(side=tk.RIGHT, fill=tk.Y)
        hsb.pack(side=tk.BOTTOM, fill=tk.X)

        # Threshold Input
        input_frame = tk.Frame(self.root, bg="#f8f9fa")
        input_frame.pack(pady=20)

        tk.Label(input_frame, text="Minimum Number of Hits:", font=("Helvetica", 13, "bold"), bg="#f8f9fa").pack(side=tk.LEFT, padx=10)
        self.threshold_entry = tk.Entry(input_frame, font=("Helvetica", 13), width=10, justify="center", relief="sunken", bd=3)
        self.threshold_entry.insert(0, "5")
        self.threshold_entry.pack(side=tk.LEFT, padx=10)

        tk.Button(input_frame, text="Generate Beautiful Plots", command=self.apply_threshold,
                  bg="#e74c3c", fg="white", font=("Helvetica", 14, "bold"), width=20, height=2).pack(side=tk.LEFT, padx=20)

    # ... [rest of the methods unchanged: select_folder, on_file_select, load_excel_file, display_data, apply_threshold, generate_plots] ...

    def select_folder(self):
        self.folder_path = filedialog.askdirectory(title="Select Folder with Excel Reports")
        if not self.folder_path:
            return
        self.file_list = [f for f in os.listdir(self.folder_path) if f.lower().endswith(('.xlsx', '.xls'))]
        self.file_listbox.delete(0, tk.END)
        for file in self.file_list:
            self.file_listbox.insert(tk.END, file)
        messagebox.showinfo("Success", f"Found {len(self.file_list)} Excel file(s)")

    def on_file_select(self, event):
        selection = self.file_listbox.curselection()
        if not selection:
            return
        self.selected_file = self.file_list[selection[0]]
        self.load_excel_file()

    def load_excel_file(self):
        file_path = os.path.join(self.folder_path, self.selected_file)
        try:
            xl = pd.ExcelFile(file_path)
            if "Distilled Report" not in xl.sheet_names:
                messagebox.showerror("Error", "Sheet 'Distilled Report' not found!")
                return
            self.df = pd.read_excel(file_path, sheet_name="Distilled Report")

            required = ["full pathogen name", "Number of Hits", "Pathogen Type",
                       "coverage ratio", "similarity value", "quality of read"]
            missing = [c for c in required if c not in self.df.columns]
            if missing:
                messagebox.showerror("Missing Columns", f"Missing: {', '.join(missing)}")
                return

            self.df = self.df[required].copy()
            self.df["Number of Hits"] = pd.to_numeric(self.df["Number of Hits"], errors='coerce').fillna(0).astype(int)

            self.display_data()
            messagebox.showinfo("Loaded", f"Loaded {len(self.df)} rows from {self.selected_file}")

        except Exception as e:
            messagebox.showerror("Error", str(e))

    def display_data(self):
        self.tree.delete(*self.tree.get_children())
        self.tree["columns"] = list(self.df.columns)
        for col in self.df.columns:
            self.tree.heading(col, text=col)
            self.tree.column(col, width=180, anchor="center")
        for _, row in self.df.iterrows():
            values = [str(v) if pd.notna(v) else "" for v in row]
            self.tree.insert("", "end", values=values)

    def apply_threshold(self):
        if self.df is None:
            messagebox.showwarning("No File", "Select a file first!")
            return
        try:
            threshold = int(self.threshold_entry.get())
            if threshold < 0: raise ValueError
        except:
            messagebox.showerror("Invalid", "Enter a valid number")
            return

        self.filtered_df = self.df[self.df["Number of Hits"] > threshold]
        if self.filtered_df.empty:
            messagebox.showinfo("No Data", f"No pathogens with >{threshold} hits")
            return
        self.generate_plots(threshold)

    def generate_plots(self, threshold):
        df = self.filtered_df
        top = df.groupby("full pathogen name")["Number of Hits"].sum().sort_values(ascending=False).head(15)
        colors = plt.cm.tab20(np.linspace(0, 1, len(top))) if len(top)>10 else plt.cm.Set3(np.linspace(0, 1, len(top)))

        plt.style.use('default')
        fig = plt.figure(figsize=(18, 11), dpi=120)
        fig.suptitle(f"Pathogen Analysis – Hits > {threshold} – {self.selected_file}", fontsize=20, fontweight='bold', y=0.98)

        # Pie
        plt.subplot(2,2,1)
        plt.pie(top.values, labels=top.index, autopct='%1.1f%%', colors=colors, textprops={'fontsize':10})
        plt.title("Share by Total Hits", fontweight='bold', fontsize=14)

        # Bar
        plt.subplot(2,2,2)
        bars = plt.barh(range(len(top)), top.values, color=colors)
        plt.yticks(range(len(top)), top.index)
        plt.gca().invert_yaxis()
        plt.title("Top Pathogens by Hits", fontweight='bold')
        for i, v in enumerate(top.values):
            plt.text(v, i, f" {v:,}", va='center', fontweight='bold')

        # Stacked metrics
        plt.subplot(2,2,3)
        metrics = df.groupby("full pathogen name")[["coverage ratio", "similarity value", "quality of read"]].mean().loc[top.index]
        scaled = metrics.div(metrics.max(axis=0), axis=1) * 100
        bottom = np.zeros(len(scaled))
        for i, col in enumerate(scaled.columns):
            plt.bar(scaled.index, scaled[col], bottom=bottom, label=col.replace(' ', '\n'), color=['#e74c3c','#3498db','#2ecc71'][i])
            bottom += scaled[col].values
        plt.xticks(rotation=50, ha='right')
        plt.legend()
        plt.title("Avg Coverage / Similarity / Quality", fontweight='bold')

        # Type pie
        plt.subplot(2,2,4)
        types = df["Pathogen Type"].value_counts()
        plt.pie(types.values, labels=types.index, autopct='%1.1f%%', colors=plt.cm.Pastel1(range(len(types))))
        plt.title("Pathogen Type Diversity", fontweight='bold')

        plt.tight_layout()
        plt.show()

    def run(self):
        self.root.mainloop()

if __name__ == "__main__":
    app = PathogenReportAnalyzer()
    app.run()