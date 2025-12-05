#!/usr/bin/env python3
"""
Phylogenetic Tree Analysis Pipeline - HPC Compatible Version
=============================================================
Designed for headless HPC environments without X11 display.
Avoids matplotlib C++ library conflicts by using pure Python/SVG output.

Author: Generated for HPC bioinformatics analysis
"""

import os
import glob
import subprocess
import tempfile
import xml.etree.ElementTree as ET
from pathlib import Path
from typing import List, Dict, Optional, Tuple
from dataclasses import dataclass, field
from abc import ABC, abstractmethod
import logging
import math

# Bioinformatics imports (pure Python, no C dependencies)
from Bio import SeqIO, AlignIO, Phylo
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Phylo.TreeConstruction import DistanceCalculator, DistanceTreeConstructor

# Configure logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s - %(message)s'
)
logger = logging.getLogger(__name__)


# =============================================================================
# DATA CLASSES
# =============================================================================

@dataclass
class SequenceData:
    """Data class to hold sequence information."""
    id: str
    sequence: str
    source_file: str
    barcode: Optional[str] = None
    description: str = ""
    
    @property
    def length(self) -> int:
        return len(self.sequence)
    
    def to_seqrecord(self) -> SeqRecord:
        """Convert to BioPython SeqRecord."""
        return SeqRecord(
            Seq(self.sequence),
            id=self.id,
            description=self.description
        )


@dataclass
class PhyloTreeResult:
    """Data class to hold phylogenetic analysis results."""
    tree: Phylo.BaseTree.Tree
    alignment_file: str
    tree_file: str
    method: str
    num_sequences: int
    alignment_length: int
    metadata: Dict = field(default_factory=dict)


# =============================================================================
# SEQUENCE LOADERS
# =============================================================================

class SequenceLoader(ABC):
    """Abstract base class for sequence loading strategies."""
    
    @abstractmethod
    def load(self, path: str) -> List[SequenceData]:
        pass


class FastaLoader(SequenceLoader):
    """Loader for FASTA format sequences."""
    
    def load(self, path: str) -> List[SequenceData]:
        sequences = []
        try:
            for record in SeqIO.parse(path, "fasta"):
                seq_data = SequenceData(
                    id=record.id,
                    sequence=str(record.seq),
                    source_file=path,
                    description=record.description
                )
                sequences.append(seq_data)
            logger.info(f"Loaded {len(sequences)} sequences from {path}")
        except Exception as e:
            logger.error(f"Error loading {path}: {e}")
        return sequences


class BarcodeDirectoryLoader(SequenceLoader):
    """Loader for barcode directory structure (barcode01, barcode02, etc.)."""
    
    def __init__(self, fasta_pattern: str = "*.fasta"):
        self.fasta_pattern = fasta_pattern
        self.fasta_loader = FastaLoader()
    
    def load(self, base_path: str) -> List[SequenceData]:
        sequences = []
        base_path = Path(base_path)
        
        # Find all barcode directories
        barcode_dirs = sorted(glob.glob(str(base_path / "barcode*")))
        
        if not barcode_dirs:
            # Try loading directly if no barcode subdirs
            fasta_files = glob.glob(str(base_path / self.fasta_pattern))
            for fasta_file in fasta_files:
                sequences.extend(self.fasta_loader.load(fasta_file))
            return sequences
        
        for barcode_dir in barcode_dirs:
            barcode_name = os.path.basename(barcode_dir)
            fasta_files = glob.glob(os.path.join(barcode_dir, self.fasta_pattern))
            
            for fasta_file in fasta_files:
                file_sequences = self.fasta_loader.load(fasta_file)
                for seq in file_sequences:
                    seq.barcode = barcode_name
                    if barcode_name not in seq.id:
                        seq.id = f"{barcode_name}_{seq.id}"
                sequences.extend(file_sequences)
        
        logger.info(f"Total sequences loaded from barcodes: {len(sequences)}")
        return sequences


class ReferenceLoader(SequenceLoader):
    """Loader for reference genome sequences."""
    
    def __init__(self, fasta_pattern: str = "*.fasta"):
        self.fasta_pattern = fasta_pattern
        self.fasta_loader = FastaLoader()
    
    def load(self, directory: str) -> List[SequenceData]:
        sequences = []
        fasta_files = glob.glob(os.path.join(directory, self.fasta_pattern))
        
        for fasta_file in fasta_files:
            file_sequences = self.fasta_loader.load(fasta_file)
            for seq in file_sequences:
                seq.description = f"[REF] {seq.description}"
            sequences.extend(file_sequences)
        
        logger.info(f"Loaded {len(sequences)} reference sequences")
        return sequences


# =============================================================================
# MULTIPLE SEQUENCE ALIGNMENT
# =============================================================================

class MultipleSequenceAligner(ABC):
    """Abstract base class for MSA methods."""
    
    @abstractmethod
    def align(self, sequences: List[SequenceData], output_file: str) -> str:
        pass


class MafftAligner(MultipleSequenceAligner):
    """MAFFT-based multiple sequence aligner."""
    
    def __init__(self, 
                 algorithm: str = "auto",
                 max_iterate: int = 1000,
                 thread: int = -1):
        self.algorithm = algorithm
        self.max_iterate = max_iterate
        self.thread = thread
    
    def align(self, sequences: List[SequenceData], output_file: str) -> str:
        # Create temporary input file
        with tempfile.NamedTemporaryFile(mode='w', suffix='.fasta', delete=False) as tmp:
            for seq in sequences:
                SeqIO.write(seq.to_seqrecord(), tmp, "fasta")
            input_file = tmp.name
        
        try:
            # Build MAFFT command
            cmd = ["mafft"]
            
            if self.algorithm == "linsi":
                cmd.extend(["--localpair", "--maxiterate", str(self.max_iterate)])
            elif self.algorithm == "ginsi":
                cmd.extend(["--globalpair", "--maxiterate", str(self.max_iterate)])
            elif self.algorithm == "einsi":
                cmd.extend(["--genafpair", "--maxiterate", str(self.max_iterate)])
            elif self.algorithm == "auto":
                cmd.append("--auto")
            
            if self.thread > 0:
                cmd.extend(["--thread", str(self.thread)])
            
            cmd.append(input_file)
            
            logger.info(f"Running MAFFT alignment with {len(sequences)} sequences...")
            result = subprocess.run(
                cmd, 
                capture_output=True, 
                text=True,
                check=True
            )
            
            with open(output_file, 'w') as f:
                f.write(result.stdout)
            
            logger.info(f"Alignment saved to {output_file}")
            return output_file
            
        except subprocess.CalledProcessError as e:
            logger.error(f"MAFFT failed: {e.stderr}")
            raise
        except FileNotFoundError:
            logger.warning("MAFFT not found, trying muscle or clustalo...")
            return self._try_alternative_aligners(input_file, output_file, sequences)
        finally:
            if os.path.exists(input_file):
                os.unlink(input_file)
    
    def _try_alternative_aligners(self, input_file: str, output_file: str, 
                                   sequences: List[SequenceData]) -> str:
        """Try alternative aligners if MAFFT unavailable."""
        # Try MUSCLE
        try:
            result = subprocess.run(
                ["muscle", "-align", input_file, "-output", output_file],
                capture_output=True, text=True, check=True
            )
            logger.info("Used MUSCLE for alignment")
            return output_file
        except:
            pass
        
        # Try Clustal Omega
        try:
            result = subprocess.run(
                ["clustalo", "-i", input_file, "-o", output_file, "--force"],
                capture_output=True, text=True, check=True
            )
            logger.info("Used Clustal Omega for alignment")
            return output_file
        except:
            pass
        
        # Last resort: copy input (sequences must be same length)
        logger.warning("No aligner found! Using input sequences as-is")
        import shutil
        shutil.copy(input_file, output_file)
        return output_file


# =============================================================================
# TREE BUILDERS
# =============================================================================

class TreeBuilder(ABC):
    """Abstract base class for phylogenetic tree construction."""
    
    @abstractmethod
    def build(self, alignment_file: str) -> Phylo.BaseTree.Tree:
        pass


class NeighborJoiningTreeBuilder(TreeBuilder):
    """Neighbor-Joining tree construction."""
    
    def __init__(self, distance_model: str = "identity"):
        self.distance_model = distance_model
    
    def build(self, alignment_file: str) -> Phylo.BaseTree.Tree:
        logger.info("Building Neighbor-Joining tree...")
        
        alignment = AlignIO.read(alignment_file, "fasta")
        calculator = DistanceCalculator(self.distance_model)
        distance_matrix = calculator.get_distance(alignment)
        
        constructor = DistanceTreeConstructor()
        tree = constructor.nj(distance_matrix)
        
        logger.info("NJ tree construction complete")
        return tree


class UPGMATreeBuilder(TreeBuilder):
    """UPGMA tree construction."""
    
    def __init__(self, distance_model: str = "identity"):
        self.distance_model = distance_model
    
    def build(self, alignment_file: str) -> Phylo.BaseTree.Tree:
        logger.info("Building UPGMA tree...")
        
        alignment = AlignIO.read(alignment_file, "fasta")
        calculator = DistanceCalculator(self.distance_model)
        distance_matrix = calculator.get_distance(alignment)
        
        constructor = DistanceTreeConstructor()
        tree = constructor.upgma(distance_matrix)
        
        logger.info("UPGMA tree construction complete")
        return tree


# =============================================================================
# PURE SVG TREE VISUALIZER (NO MATPLOTLIB!)
# =============================================================================

class PureSVGTreeVisualizer:
    """
    Pure Python SVG tree visualizer - NO matplotlib required!
    Works in headless HPC environments.
    """
    
    def __init__(self, 
                 width: int = 800,
                 height: int = 600,
                 margin: int = 50,
                 font_size: int = 12,
                 line_width: int = 2):
        self.width = width
        self.height = height
        self.margin = margin
        self.font_size = font_size
        self.line_width = line_width
        
        # Colors
        self.bg_color = "#ffffff"
        self.line_color = "#4a5568"
        self.node_color = "#4299e1"
        self.text_color = "#2d3748"
        self.ref_color = "#38a169"
    
    def visualize(self, tree: Phylo.BaseTree.Tree, output_file: str, 
                  title: str = "Phylogenetic Tree") -> str:
        """Generate SVG visualization of the tree."""
        
        # Get all terminal nodes (tips)
        terminals = tree.get_terminals()
        num_tips = len(terminals)
        
        # Calculate layout
        plot_width = self.width - 2 * self.margin - 150  # Leave space for labels
        plot_height = self.height - 2 * self.margin - 40  # Leave space for title
        
        # Get tree depth for x-scaling
        max_depth = max(tree.distance(tree.root, t) for t in terminals)
        if max_depth == 0:
            max_depth = 1
        
        # Assign y-positions to terminals
        y_positions = {}
        for i, term in enumerate(terminals):
            y_positions[term] = self.margin + 40 + (i * plot_height / (num_tips - 1 if num_tips > 1 else 1))
        
        # Calculate internal node positions (average of children)
        def get_y_position(clade):
            if clade in y_positions:
                return y_positions[clade]
            children_y = [get_y_position(child) for child in clade.clades]
            y_pos = sum(children_y) / len(children_y)
            y_positions[clade] = y_pos
            return y_pos
        
        # Calculate y for all nodes
        get_y_position(tree.root)
        
        # Calculate x-positions based on distance from root
        x_positions = {}
        def get_x_position(clade, parent_dist=0):
            dist = parent_dist + (clade.branch_length or 0)
            x_pos = self.margin + (dist / max_depth) * plot_width
            x_positions[clade] = x_pos
            for child in clade.clades:
                get_x_position(child, dist)
        
        x_positions[tree.root] = self.margin
        for child in tree.root.clades:
            get_x_position(child, 0)
        
        # Start building SVG
        svg_elements = []
        
        # SVG header
        svg_elements.append(f'''<svg xmlns="http://www.w3.org/2000/svg" 
            width="{self.width}" height="{self.height}" 
            viewBox="0 0 {self.width} {self.height}">''')
        
        # Background
        svg_elements.append(f'<rect width="100%" height="100%" fill="{self.bg_color}"/>')
        
        # Title
        svg_elements.append(f'''<text x="{self.width/2}" y="30" 
            text-anchor="middle" font-family="Arial, sans-serif" 
            font-size="18" font-weight="bold" fill="{self.text_color}">{title}</text>''')
        
        # Draw edges
        def draw_edges(clade):
            if clade.clades:
                x = x_positions[clade]
                y = y_positions[clade]
                
                for child in clade.clades:
                    child_x = x_positions[child]
                    child_y = y_positions[child]
                    
                    # Draw horizontal then vertical line (cladogram style)
                    svg_elements.append(f'''<path d="M {x} {y} L {x} {child_y} L {child_x} {child_y}" 
                        fill="none" stroke="{self.line_color}" stroke-width="{self.line_width}"/>''')
                    
                    draw_edges(child)
        
        draw_edges(tree.root)
        
        # Draw nodes and labels
        for clade in tree.find_clades():
            x = x_positions[clade]
            y = y_positions[clade]
            
            # Node circle
            color = self.node_color
            if clade.name and "REF" in str(clade.name):
                color = self.ref_color
            
            svg_elements.append(f'''<circle cx="{x}" cy="{y}" r="4" 
                fill="{color}" stroke="{color}" stroke-width="1"/>''')
            
            # Label for terminals
            if clade.is_terminal() and clade.name:
                label = clade.name
                svg_elements.append(f'''<text x="{x + 10}" y="{y + 4}" 
                    font-family="Arial, sans-serif" font-size="{self.font_size}" 
                    fill="{self.text_color}">{label}</text>''')
        
        # Close SVG
        svg_elements.append('</svg>')
        
        # Write to file
        svg_content = '\n'.join(svg_elements)
        with open(output_file, 'w') as f:
            f.write(svg_content)
        
        logger.info(f"SVG tree saved to {output_file}")
        return output_file
    
    def visualize_circular(self, tree: Phylo.BaseTree.Tree, output_file: str,
                           title: str = "Phylogenetic Tree (Circular)") -> str:
        """Generate circular SVG visualization."""
        
        terminals = tree.get_terminals()
        num_tips = len(terminals)
        
        cx = self.width / 2
        cy = self.height / 2
        max_radius = min(cx, cy) - self.margin - 80
        
        # Get tree depth
        max_depth = max(tree.distance(tree.root, t) for t in terminals)
        if max_depth == 0:
            max_depth = 1
        
        # Assign angles to terminals
        angles = {}
        for i, term in enumerate(terminals):
            angles[term] = (2 * math.pi * i / num_tips) - math.pi / 2
        
        # Calculate internal node angles (average of children)
        def get_angle(clade):
            if clade in angles:
                return angles[clade]
            children_angles = [get_angle(child) for child in clade.clades]
            angle = sum(children_angles) / len(children_angles)
            angles[clade] = angle
            return angle
        
        get_angle(tree.root)
        
        # Calculate radii based on distance
        radii = {tree.root: 0}
        def get_radius(clade, parent_dist=0):
            dist = parent_dist + (clade.branch_length or 0)
            radii[clade] = (dist / max_depth) * max_radius
            for child in clade.clades:
                get_radius(child, dist)
        
        for child in tree.root.clades:
            get_radius(child, 0)
        
        # Convert to cartesian
        def polar_to_cart(r, theta):
            return cx + r * math.cos(theta), cy + r * math.sin(theta)
        
        # Build SVG
        svg_elements = []
        svg_elements.append(f'''<svg xmlns="http://www.w3.org/2000/svg" 
            width="{self.width}" height="{self.height}">''')
        svg_elements.append(f'<rect width="100%" height="100%" fill="{self.bg_color}"/>')
        svg_elements.append(f'''<text x="{cx}" y="25" text-anchor="middle" 
            font-family="Arial" font-size="16" font-weight="bold" 
            fill="{self.text_color}">{title}</text>''')
        
        # Draw edges
        def draw_circular_edges(clade):
            if clade.clades:
                r = radii[clade]
                theta = angles[clade]
                x, y = polar_to_cart(r, theta)
                
                for child in clade.clades:
                    child_r = radii[child]
                    child_theta = angles[child]
                    child_x, child_y = polar_to_cart(child_r, child_theta)
                    
                    # Arc from parent to child angle at parent radius, then line to child
                    mid_x, mid_y = polar_to_cart(r, child_theta)
                    svg_elements.append(f'''<path d="M {x} {y} L {mid_x} {mid_y} L {child_x} {child_y}" 
                        fill="none" stroke="{self.line_color}" stroke-width="{self.line_width}"/>''')
                    
                    draw_circular_edges(child)
        
        draw_circular_edges(tree.root)
        
        # Draw nodes and labels
        for clade in tree.find_clades():
            r = radii[clade]
            theta = angles[clade]
            x, y = polar_to_cart(r, theta)
            
            color = self.ref_color if (clade.name and "REF" in str(clade.name)) else self.node_color
            svg_elements.append(f'<circle cx="{x}" cy="{y}" r="3" fill="{color}"/>')
            
            if clade.is_terminal() and clade.name:
                # Rotate label based on position
                label_r = r + 10
                lx, ly = polar_to_cart(label_r, theta)
                rotation = math.degrees(theta)
                if rotation > 90 or rotation < -90:
                    rotation += 180
                    anchor = "end"
                else:
                    anchor = "start"
                
                svg_elements.append(f'''<text x="{lx}" y="{ly}" 
                    transform="rotate({rotation}, {lx}, {ly})"
                    text-anchor="{anchor}" font-family="Arial" 
                    font-size="10" fill="{self.text_color}">{clade.name}</text>''')
        
        svg_elements.append('</svg>')
        
        with open(output_file, 'w') as f:
            f.write('\n'.join(svg_elements))
        
        logger.info(f"Circular SVG tree saved to {output_file}")
        return output_file


# =============================================================================
# ASCII TREE VISUALIZER (Terminal output)
# =============================================================================

class ASCIITreeVisualizer:
    """ASCII tree visualization for terminal output."""
    
    def visualize(self, tree: Phylo.BaseTree.Tree, output_file: Optional[str] = None) -> str:
        """Generate ASCII tree representation."""
        from io import StringIO
        
        output = StringIO()
        Phylo.draw_ascii(tree, file=output)
        ascii_tree = output.getvalue()
        
        if output_file:
            with open(output_file, 'w') as f:
                f.write(ascii_tree)
            logger.info(f"ASCII tree saved to {output_file}")
        
        return ascii_tree


# =============================================================================
# HTML REPORT GENERATOR
# =============================================================================

class HTMLReportGenerator:
    """Generate HTML reports with embedded SVG."""
    
    def __init__(self, title: str = "Phylogenetic Analysis Report"):
        self.title = title
    
    def generate(self, output_file: str, tree_result: PhyloTreeResult,
                 sequences: List[SequenceData], svg_file: Optional[str] = None) -> str:
        
        # Read SVG if provided
        svg_content = ""
        if svg_file and os.path.exists(svg_file):
            with open(svg_file, 'r') as f:
                svg_content = f.read()
        
        # Count references
        num_refs = sum(1 for s in sequences if "[REF]" in s.description)
        num_samples = len(sequences) - num_refs
        
        # Build sequence table
        seq_rows = ""
        for seq in sequences:
            is_ref = "[REF]" in seq.description
            badge = '<span class="badge-ref">Reference</span>' if is_ref else '<span class="badge-sample">Sample</span>'
            seq_rows += f"""
            <tr>
                <td>{seq.id}</td>
                <td>{seq.length:,}</td>
                <td>{os.path.basename(seq.source_file)}</td>
                <td>{badge}</td>
            </tr>"""
        
        html = f"""<!DOCTYPE html>
<html lang="en">
<head>
    <meta charset="UTF-8">
    <title>{self.title}</title>
    <style>
        * {{ margin: 0; padding: 0; box-sizing: border-box; }}
        body {{ 
            font-family: -apple-system, BlinkMacSystemFont, 'Segoe UI', Roboto, sans-serif;
            background: #f5f5f5; color: #333; line-height: 1.6;
        }}
        .header {{
            background: linear-gradient(135deg, #2c5282 0%, #1a365d 100%);
            color: white; padding: 2rem; text-align: center;
        }}
        .header h1 {{ font-weight: 300; font-size: 2rem; }}
        .container {{ max-width: 1200px; margin: 0 auto; padding: 2rem; }}
        .stats {{ display: flex; gap: 1rem; margin: 1.5rem 0; flex-wrap: wrap; }}
        .stat {{
            background: white; padding: 1.5rem; border-radius: 8px;
            text-align: center; flex: 1; min-width: 150px;
            box-shadow: 0 2px 4px rgba(0,0,0,0.1);
        }}
        .stat-value {{ font-size: 2rem; color: #2c5282; font-weight: bold; }}
        .stat-label {{ color: #666; font-size: 0.9rem; }}
        .section {{
            background: white; border-radius: 8px; padding: 1.5rem;
            margin: 1.5rem 0; box-shadow: 0 2px 4px rgba(0,0,0,0.1);
        }}
        .section h2 {{
            color: #2c5282; margin-bottom: 1rem;
            padding-bottom: 0.5rem; border-bottom: 2px solid #4299e1;
        }}
        .tree-container {{
            background: #fafafa; border-radius: 4px;
            padding: 1rem; overflow-x: auto; text-align: center;
        }}
        .tree-container svg {{ max-width: 100%; height: auto; }}
        table {{ width: 100%; border-collapse: collapse; }}
        th, td {{ padding: 0.75rem; text-align: left; border-bottom: 1px solid #eee; }}
        th {{ background: #f8f9fa; font-weight: 600; }}
        .badge-ref {{
            background: #c6f6d5; color: #22543d;
            padding: 0.25rem 0.5rem; border-radius: 4px; font-size: 0.8rem;
        }}
        .badge-sample {{
            background: #bee3f8; color: #2a4365;
            padding: 0.25rem 0.5rem; border-radius: 4px; font-size: 0.8rem;
        }}
        .footer {{ text-align: center; padding: 2rem; color: #666; }}
    </style>
</head>
<body>
    <div class="header">
        <h1>🧬 {self.title}</h1>
        <p>HPC Phylogenetic Analysis Pipeline</p>
    </div>
    
    <div class="container">
        <div class="stats">
            <div class="stat">
                <div class="stat-value">{num_samples}</div>
                <div class="stat-label">Sample Sequences</div>
            </div>
            <div class="stat">
                <div class="stat-value">{num_refs}</div>
                <div class="stat-label">References</div>
            </div>
            <div class="stat">
                <div class="stat-value">{tree_result.alignment_length:,}</div>
                <div class="stat-label">Alignment Length</div>
            </div>
            <div class="stat">
                <div class="stat-value">{tree_result.method}</div>
                <div class="stat-label">Tree Method</div>
            </div>
        </div>
        
        <div class="section">
            <h2>🌳 Phylogenetic Tree</h2>
            <div class="tree-container">
                {svg_content if svg_content else '<p>Tree visualization not available</p>'}
            </div>
        </div>
        
        <div class="section">
            <h2>📊 Sequence Summary</h2>
            <table>
                <thead>
                    <tr><th>Sequence ID</th><th>Length</th><th>Source</th><th>Type</th></tr>
                </thead>
                <tbody>{seq_rows}</tbody>
            </table>
        </div>
        
        <div class="section">
            <h2>📁 Output Files</h2>
            <ul>
                <li><strong>Alignment:</strong> {os.path.basename(tree_result.alignment_file)}</li>
                <li><strong>Tree (Newick):</strong> {os.path.basename(tree_result.tree_file)}</li>
            </ul>
        </div>
    </div>
    
    <div class="footer">
        <p>Generated by HPC Phylogenetic Analysis Pipeline</p>
    </div>
</body>
</html>"""
        
        with open(output_file, 'w') as f:
            f.write(html)
        
        logger.info(f"HTML report saved to {output_file}")
        return output_file


# =============================================================================
# MAIN PIPELINE
# =============================================================================

class PhylogeneticAnalysisPipeline:
    """Main pipeline orchestrating the entire phylogenetic analysis workflow."""
    
    def __init__(self,
                 sample_dir: str,
                 reference_dir: str,
                 output_dir: str,
                 aligner: Optional[MultipleSequenceAligner] = None,
                 tree_builder: Optional[TreeBuilder] = None):
        
        self.sample_dir = Path(sample_dir)
        self.reference_dir = Path(reference_dir)
        self.output_dir = Path(output_dir)
        
        self.output_dir.mkdir(parents=True, exist_ok=True)
        
        self.sample_loader = BarcodeDirectoryLoader()
        self.reference_loader = ReferenceLoader()
        self.aligner = aligner or MafftAligner()
        self.tree_builder = tree_builder or NeighborJoiningTreeBuilder()
        self.visualizer = PureSVGTreeVisualizer()
        
        self.sequences: List[SequenceData] = []
        self.result: Optional[PhyloTreeResult] = None
    
    def load_sequences(self) -> List[SequenceData]:
        """Load all sequences from samples and references."""
        logger.info("Loading sequences...")
        
        sample_sequences = self.sample_loader.load(str(self.sample_dir))
        reference_sequences = self.reference_loader.load(str(self.reference_dir))
        
        self.sequences = sample_sequences + reference_sequences
        
        logger.info(f"Total sequences: {len(self.sequences)} "
                   f"(Samples: {len(sample_sequences)}, References: {len(reference_sequences)})")
        
        return self.sequences
    
    def run(self, title: str = "Phylogenetic Analysis") -> PhyloTreeResult:
        """Execute the complete analysis pipeline."""
        
        print("\n" + "="*60)
        print("🧬 PHYLOGENETIC ANALYSIS PIPELINE (HPC Compatible)")
        print("="*60 + "\n")
        
        # Step 1: Load sequences
        self.load_sequences()
        
        if len(self.sequences) < 3:
            raise ValueError(f"Insufficient sequences ({len(self.sequences)}). Need at least 3.")
        
        # Step 2: Alignment
        alignment_file = str(self.output_dir / "alignment.fasta")
        self.aligner.align(self.sequences, alignment_file)
        
        # Get alignment length
        try:
            alignment = AlignIO.read(alignment_file, "fasta")
            alignment_length = alignment.get_alignment_length()
        except:
            alignment_length = len(self.sequences[0].sequence)
        
        # Step 3: Build tree
        tree = self.tree_builder.build(alignment_file)
        
        # Step 4: Save tree
        tree_file = str(self.output_dir / "tree.nwk")
        Phylo.write(tree, tree_file, "newick")
        
        # Step 5: Visualize (Pure SVG - no matplotlib!)
        svg_file = str(self.output_dir / "tree.svg")
        self.visualizer.visualize(tree, svg_file, title=title)
        
        # Circular version too
        svg_circular = str(self.output_dir / "tree_circular.svg")
        self.visualizer.visualize_circular(tree, svg_circular, title=f"{title} (Circular)")
        
        # ASCII for terminal viewing
        ascii_viz = ASCIITreeVisualizer()
        ascii_file = str(self.output_dir / "tree_ascii.txt")
        ascii_tree = ascii_viz.visualize(tree, ascii_file)
        
        # Print ASCII tree to terminal
        print("\n" + "="*60)
        print("ASCII TREE VISUALIZATION:")
        print("="*60)
        print(ascii_tree)
        
        # Create result
        self.result = PhyloTreeResult(
            tree=tree,
            alignment_file=alignment_file,
            tree_file=tree_file,
            method=self.tree_builder.__class__.__name__.replace("TreeBuilder", ""),
            num_sequences=len(self.sequences),
            alignment_length=alignment_length,
            metadata={'svg_file': svg_file, 'svg_circular': svg_circular}
        )
        
        # Step 6: Generate HTML report
        report_gen = HTMLReportGenerator(title=title)
        report_file = str(self.output_dir / "report.html")
        report_gen.generate(report_file, self.result, self.sequences, svg_file)
        
        print("\n" + "="*60)
        print("✅ PIPELINE COMPLETE!")
        print(f"📁 Results saved to: {self.output_dir}")
        print(f"   - tree.nwk (Newick format)")
        print(f"   - tree.svg (Vector graphics)")
        print(f"   - tree_circular.svg (Circular layout)")
        print(f"   - tree_ascii.txt (Terminal viewable)")
        print(f"   - report.html (Full report)")
        print("="*60 + "\n")
        
        return self.result


# =============================================================================
# FACTORY FUNCTION
# =============================================================================

def create_pipeline(sample_dir: str,
                   reference_dir: str,
                   output_dir: str,
                   alignment_method: str = "mafft",
                   tree_method: str = "nj",
                   threads: int = 4) -> PhylogeneticAnalysisPipeline:
    """
    Factory function to create a configured pipeline.
    
    Args:
        sample_dir: Path to barcode directories
        reference_dir: Path to reference genomes
        output_dir: Output directory
        alignment_method: 'mafft', 'mafft_accurate', or 'mafft_fast'
        tree_method: 'nj' (Neighbor-Joining) or 'upgma'
        threads: Number of threads for alignment
    """
    # Select aligner
    if alignment_method == 'mafft_accurate':
        aligner = MafftAligner(algorithm='linsi', thread=threads)
    elif alignment_method == 'mafft_fast':
        aligner = MafftAligner(algorithm='fftns', thread=threads)
    else:
        aligner = MafftAligner(algorithm='auto', thread=threads)
    
    # Select tree builder
    if tree_method == 'upgma':
        tree_builder = UPGMATreeBuilder()
    else:
        tree_builder = NeighborJoiningTreeBuilder()
    
    return PhylogeneticAnalysisPipeline(
        sample_dir=sample_dir,
        reference_dir=reference_dir,
        output_dir=output_dir,
        aligner=aligner,
        tree_builder=tree_builder
    )


# =============================================================================
# MAIN
# =============================================================================

def main():
    """Main function - update paths for your data."""
    
    # =========================================
    # UPDATE THESE PATHS FOR YOUR HPC DATA
    # =========================================
    SAMPLE_DIR = "../Data/Denovo/reformed"
    REFERENCE_DIR = "../Data/complete_genomes_fasta_format"
    OUTPUT_DIR = "./phylo_results"
    
    # Check if paths exist
    if not os.path.exists(SAMPLE_DIR):
        print(f"⚠️  Sample directory not found: {SAMPLE_DIR}")
        print("   Please update SAMPLE_DIR in the script.")
        print("\n   Running demo with mock data instead...\n")
        run_demo()
        return
    
    if not os.path.exists(REFERENCE_DIR):
        print(f"⚠️  Reference directory not found: {REFERENCE_DIR}")
        print("   Please update REFERENCE_DIR in the script.")
        return
    
    # Create and run pipeline
    pipeline = create_pipeline(
        sample_dir=SAMPLE_DIR,
        reference_dir=REFERENCE_DIR,
        output_dir=OUTPUT_DIR,
        alignment_method="mafft",
        tree_method="nj",
        threads=4  # Adjust based on your HPC allocation
    )
    
    result = pipeline.run(title="Pathogen Phylogenetic Analysis")
    print(f"\n📊 Analysis Summary:")
    print(f"   Sequences: {result.num_sequences}")
    print(f"   Alignment length: {result.alignment_length}")
    print(f"   Method: {result.method}")


def run_demo():
    """Run demonstration with mock data."""
    
    demo_dir = Path("./demo_data")
    demo_dir.mkdir(exist_ok=True)
    
    # Create mock barcode directories with SAME LENGTH sequences
    # (important when no aligner is available)
    base_seq = "ATGCGATCGATCGATCGATCGATCGATCGATCGACGTACGTACGTACGTACGT"
    
    sample_seqs = [
        ("sample_1", "ATGCGATCGATCGATCGATCGATCGATCGATCGACGTACGTACGTACGTACGT"),
        ("sample_2", "ATGCGATCGATCGTTCGATCGATCGATCGATCGACGTACGTACGTACGTACGT"),
        ("sample_3", "ATGCGATCGAACGATCGATCGATCGATCGATCGACGTACGTACGTACGTACGT"),
        ("sample_4", "ATGCGATCGATCGATCGATCGTTCGATCGATCGACGTACGTACGTACGTACGT"),
    ]
    
    for i, (name, seq) in enumerate(sample_seqs, 1):
        barcode_dir = demo_dir / f"barcode0{i}"
        barcode_dir.mkdir(exist_ok=True)
        record = SeqRecord(Seq(seq), id=name, description="")
        SeqIO.write([record], str(barcode_dir / "sequences.fasta"), "fasta")
    
    # Create mock references (same length)
    ref_dir = demo_dir / "references"
    ref_dir.mkdir(exist_ok=True)
    
    refs = [
        ("PathogenA", "ATGCGATCGATCGATCGATCGATCGATCGATCGACGTACGTACGTACGTACGT"),
        ("PathogenB", "TTGCGATCGATCGATCGATCGATCGATCGATCGACGTACGTACGTACGTACGT"),
        ("PathogenC", "ATGCGATCGATCGATCGATCGATCGATCGATCGTTTTACGTACGTACGTACGT"),
    ]
    for name, seq in refs:
        record = SeqRecord(Seq(seq), id=name, description="")
        SeqIO.write([record], str(ref_dir / f"{name}.fasta"), "fasta")
    
    # Run pipeline
    pipeline = create_pipeline(
        sample_dir=str(demo_dir),
        reference_dir=str(ref_dir),
        output_dir="./demo_output",
        threads=1
    )
    
    pipeline.run(title="Demo Phylogenetic Analysis")


if __name__ == "__main__":
    main()